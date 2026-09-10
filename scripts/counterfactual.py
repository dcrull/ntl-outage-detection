"""Small, offline native-cell radiance models. Arrays are calendar-day x cell.

M0 uncertainty is an approximate Student-t reference, not an exact median posterior.
M1 is a normal/inverse-gamma local-level filter, conditional on empirical-Bayes
scale priors and a fixed process/observation variance ratio. Neither uses labels.
"""
from dataclasses import asdict, dataclass
import warnings

import numpy as np
from scipy.stats import t as student_t


@dataclass(frozen=True)
class Settings:
    min_observations: int = 10
    scale_floor: float = 0.1  # nW cm^-2 sr^-1; provisional, not an outage threshold
    scale_prior_weight: float = 20.0
    process_ratio: float = 0.001  # daily state variance / observation variance

    def __post_init__(self):
        if self.min_observations < 3:
            raise ValueError('At least three observations are required')
        if not np.isfinite(self.scale_floor) or self.scale_floor <= 0:
            raise ValueError('Scale floor must be finite and positive')
        if not np.isfinite(self.scale_prior_weight) or self.scale_prior_weight <= 2:
            raise ValueError('Scale prior weight must exceed two')
        if not np.isfinite(self.process_ratio) or not 0 <= self.process_ratio <= 0.01:
            raise ValueError('This regularized implementation requires process_ratio in [0, .01]')


def fit_models(training, settings=Settings(), *, scale_pool=None):
    """Fit solely to supplied training data, including all missing calendar days.

    Scale pooling uses quartiles of training median brightness among eligible
    cells. M1's prior is pooled within a brightness group, not a cell's own MAD.
    Hyperparameter uncertainty is not integrated; predictive calibration is tested
    on held-out blocks. Zero is a valid observation; infinity is invalid input.
    For spatial chunks, scale_pool=(edges, scales) supplies the three brightness
    quartiles and four MAD scales computed over the full training-only AOI.
    """
    y = np.asarray(training, dtype=float)
    if y.ndim != 2 or not len(y) or np.isinf(y).any():
        raise ValueError('Expected nonempty date x cell data with NaN missingness')
    count = np.isfinite(y).sum(axis=0)
    supported = count >= settings.min_observations
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered')
        center = np.nanmedian(y, axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(y - center), axis=0)
    if scale_pool is not None:
        edges, scales = (np.asarray(v, dtype=float) for v in scale_pool)
        if (edges.shape != (3,) or scales.shape != (4,)
                or not np.isfinite(edges).all() or np.any(np.diff(edges) < 0)
                or not np.isfinite(scales).all() or np.any(scales < settings.scale_floor)):
            raise ValueError('Invalid full-AOI scale pool')
        groups = np.searchsorted(edges, center, side='right')
        pooled = scales[groups]
        group_scales = scales.tolist()
    elif supported.any():
        edges = np.quantile(center[supported], [.25, .5, .75])
        groups = np.searchsorted(edges, center, side='right')
        global_scale = max(float(np.median(mad[supported])), settings.scale_floor)
    else:
        edges = np.full(3, np.nan)
        groups = np.zeros(y.shape[1], dtype=int)
        global_scale = settings.scale_floor
    if scale_pool is None:
        pooled = np.full(y.shape[1], global_scale)
        group_scales = []
        for group in range(4):
            use = supported & (groups == group)
            scale = max(float(np.median(mad[use])), settings.scale_floor) if use.any() else global_scale
            pooled[groups == group] = scale
            group_scales.append(scale)
    weight = count / (count + settings.scale_prior_weight)
    m0_variance = np.maximum(
        weight * mad**2 + (1 - weight) * pooled**2, settings.scale_floor**2)

    # Diffuse location prior, conditioned on the first valid observation: m=y, C=1.
    # Observation variance theta has InvGamma(nu0/2, nu0*S0/2) prior.
    # C is state variance divided by theta. Q=1+C is innovation variance / theta.
    level = np.full(y.shape[1], np.nan)
    c = np.zeros(y.shape[1])
    df = np.full(y.shape[1], settings.scale_prior_weight)
    scale_sum = df * pooled**2
    for observations in y:
        initialized = np.isfinite(level)
        c[initialized] += settings.process_ratio
        valid = np.isfinite(observations)
        update = initialized & valid
        innovation = observations[update] - level[update]
        q = 1 + c[update]
        gain = c[update] / q
        level[update] += gain * innovation
        scale_sum[update] += innovation**2 / q
        df[update] += 1
        c[update] /= q
        first = valid & ~initialized
        level[first] = observations[first]
        c[first] = 1
    return {
        'settings': asdict(settings), 'count': count, 'supported': supported,
        'baseline_median': center, 'baseline_mad_sigma': mad,
        'pool_group': groups, 'pool_edges': edges, 'pool_scales': np.array(group_scales),
        'pooled_sigma': pooled,
        'M0': {'center': center, 'variance': m0_variance,
               'df': np.maximum(count - 1, 3)},
        'M1': {'center': level, 'variance': scale_sum / df, 'df': df, 'state_factor': c},
    }


def forecast(fit, method, horizons):
    """Open-loop predictions. No event/holdout observations enter this function."""
    h = np.asarray(horizons, dtype=float)
    if h.ndim != 1 or not len(h) or not np.isfinite(h).all() or np.any(h < 1):
        raise ValueError('Horizons must be positive elapsed calendar days')
    params = fit[method]
    n = np.maximum(fit['count'], 1)
    if method == 'M0':
        # Asymptotic normal-sample median uncertainty: Var(median) ~ pi*sigma²/(2n).
        factor = np.broadcast_to(1 + np.pi / (2 * n), (len(h), len(n)))
    elif method == 'M1':
        factor = 1 + params['state_factor'][None, :] + h[:, None] * fit['settings']['process_ratio']
    else:
        raise ValueError(f'Unknown method: {method}')
    center = np.broadcast_to(params['center'], factor.shape).copy()
    scale = np.sqrt(params['variance'][None, :] * factor)
    df = np.broadcast_to(params['df'], factor.shape).copy()
    supported = fit['supported'][None, :]
    center = np.where(supported, center, np.nan)
    scale = np.where(supported, scale, np.nan)
    df = np.where(supported, df, np.nan)
    # Degrees of freedom repeat across horizons and cells. Evaluate each unique
    # value once; this is the same Student-t quantile as the dense calculation.
    unique_df, inverse = np.unique(df[0], return_inverse=True)
    lo = student_t.ppf(.05, unique_df)[inverse][None, :]
    hi = student_t.ppf(.95, unique_df)[inverse][None, :]
    return {
        'center': center, 'scale': scale, 'df': df,
        'sd': scale * np.sqrt(df / (df - 2)),
        'q05': center + scale * lo,
        'q95': center + scale * hi,
    }


def evidence(prediction, observed):
    """Missing observations produce no residual, score, or tail probability."""
    y = np.asarray(observed, dtype=float)
    residual = y - prediction['center']
    standardized = residual / prediction['sd']
    tail = student_t.cdf(residual / prediction['scale'], prediction['df'])
    return {'residual': residual, 'standardized_residual': standardized, 'lower_tail': tail}


def validation_folds(baseline_days):
    """Nonoverlapping targets within each design; 7d and 14d designs are separate."""
    folds = []
    for width in (7, 14):
        for end in range(28, baseline_days - width + 1, width):
            folds.append({'design': f'{width}d', 'train_days': end, 'test_days': width})
    return folds


def metrics(prediction, observed, mask=None):
    """Descriptive held-out scores; no independence-based standard errors."""
    y = np.asarray(observed)
    use = np.isfinite(y) & np.isfinite(prediction['center'])
    if mask is not None:
        use &= mask
    result = {'n': int(use.sum())}
    names = ('mae', 'rmse', 'bias', 'coverage90', 'below_q05', 'above_q95',
             'mean_interval_width', 'interval_score90', 'mean_negative_log_density')
    if not use.any():
        return {**result, **{key: None for key in names}}
    actual = y[use]
    center, scale, df, lo, hi = (prediction[key][use] for key in ('center', 'scale', 'df', 'q05', 'q95'))
    error = actual - center
    width = hi - lo
    score = width + 20 * np.maximum(lo - actual, 0) + 20 * np.maximum(actual - hi, 0)
    result.update(mae=float(np.mean(np.abs(error))), rmse=float(np.sqrt(np.mean(error**2))),
                  bias=float(np.mean(error)), coverage90=float(np.mean((actual >= lo) & (actual <= hi))),
                  below_q05=float(np.mean(actual < lo)), above_q95=float(np.mean(actual > hi)),
                  mean_interval_width=float(np.mean(width)), interval_score90=float(np.mean(score)),
                  mean_negative_log_density=float(np.mean(
                      -student_t.logpdf(error / scale, df) + np.log(scale))))
    return result
