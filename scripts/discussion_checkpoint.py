"""Plotting and missing-aware clustering for the colleague checkpoint notebook."""
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap, BoundaryNorm, FuncNorm, SymLogNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from pyproj import Transformer
import rasterio
from rasterio.transform import Affine

UNITS = 'nW cm$^{-2}$ sr$^{-1}$'
COLORS = ['#0072B2', '#009E73', '#E69F00', '#CC79A7', '#D55E00', '#56B4E9']
TO_UTM = Transformer.from_crs(4326, 32616, always_xy=True)


def file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def median(values, axis=None):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered')
        return np.nanmedian(values, axis=axis)


def projected_box(west, south, east, north):
    x, y = TO_UTM.transform([west, west, east, east], [south, north, south, north])
    return (min(x), min(y), max(x), max(y))


@dataclass
class Scene:
    key: str
    title: str
    observed: np.ndarray
    baseline: np.ndarray
    residuals: dict
    mask: np.ndarray
    transform: Affine
    crs: int
    dates: pd.DatetimeIndex
    event_dates: pd.DatetimeIndex
    _regions: dict = field(default_factory=dict)
    _series: dict = field(default_factory=dict)

    def region(self, bounds=None):
        cache_key = None if bounds is None else tuple(bounds)
        if cache_key in self._regions:
            return self._regions[cache_key]
        h, w = self.mask.shape
        if self.crs == 32616:
            if bounds is None:
                r0, r1, c0, c1 = 0, h, 0, w
            else:
                c0 = max(0, int(np.floor((bounds[0] - self.transform.c) / self.transform.a)))
                c1 = min(w, int(np.ceil((bounds[2] - self.transform.c) / self.transform.a)))
                r0 = max(0, int(np.floor((self.transform.f - bounds[3]) / -self.transform.e)))
                r1 = min(h, int(np.ceil((self.transform.f - bounds[1]) / -self.transform.e)))
            xe = self.transform.c + np.arange(c0, c1 + 1) * self.transform.a
            ye = self.transform.f + np.arange(r0, r1 + 1) * self.transform.e
            region_mask = self.mask[r0:r1, c0:c1].copy()
            if bounds is not None:
                region_mask &= ((xe[:-1] + self.transform.a / 2 >= bounds[0]) &
                                (xe[:-1] + self.transform.a / 2 <= bounds[2]))[None, :]
                region_mask &= ((ye[:-1] + self.transform.e / 2 >= bounds[1]) &
                                (ye[:-1] + self.transform.e / 2 <= bounds[3]))[:, None]
            geometry = (xe, ye)
        else:
            xe = self.transform.c + np.arange(w + 1) * self.transform.a
            ye = self.transform.f + np.arange(h + 1) * self.transform.e
            xedges, yedges = TO_UTM.transform(*np.meshgrid(xe, ye))
            xc, yc = TO_UTM.transform(*np.meshgrid(xe[:-1] + self.transform.a / 2,
                                                  ye[:-1] + self.transform.e / 2))
            selected = np.ones((h, w), bool) if bounds is None else (
                (xc >= bounds[0]) & (xc <= bounds[2]) & (yc >= bounds[1]) & (yc <= bounds[3]))
            rr, cc = np.where(selected)
            if not len(rr):
                raise ValueError('Zoom does not intersect the native scene')
            r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
            region_mask = self.mask[r0:r1, c0:c1] & selected[r0:r1, c0:c1]
            geometry = (xedges[r0:r1 + 1, c0:c1 + 1], yedges[r0:r1 + 1, c0:c1 + 1])
        result = (slice(r0, r1), slice(c0, c1), region_mask, geometry)
        self._regions[cache_key] = result
        return result

    def values(self, kind, bounds=None, day=None, method='M0'):
        rr, cc, mask, _ = self.region(bounds)
        if kind == 'baseline':
            values = self.baseline[rr, cc]
        elif kind == 'observed':
            values = self.observed[self.dates.get_loc(pd.Timestamp(day)), rr, cc]
        elif kind == 'residual':
            values = self.residuals[method][self.event_dates.get_loc(pd.Timestamp(day)), rr, cc]
        else:
            raise ValueError(kind)
        return np.where(mask, values, np.nan)


def load_scenes(root):
    root = Path(root)
    base = root / 'data/published/reallocation_check'
    native_manifest = json.loads((base / 'manifest.json').read_text())
    if native_manifest['status'] != 'complete':
        raise ValueError('Native comparison incomplete')
    used = {}

    def verify(path, expected):
        actual = file_hash(path)
        if actual != expected:
            raise ValueError(f'Input checksum mismatch: {path}')
        used[str(path.relative_to(root))] = actual

    verify(base / 'comparison_arrays.npz', native_manifest['output_sha256']['comparison_arrays.npz'])
    with np.load(base / 'comparison_arrays.npz', allow_pickle=False) as a:
        mask = a['county'].astype(bool)
        dates = pd.DatetimeIndex(pd.to_datetime(a['dates']))
        cut = int(a['cutoff'])
        transform = Affine(*a['transform'])

        def restore(flat):
            result = np.full((len(flat), *mask.shape), np.nan, dtype='float32')
            result[:, mask] = flat
            return result

        scenes = {}
        for key, source, title in [('a2', 'direct', 'VIIRS A2 · matched native support'),
                                    ('allocated500', 'allocated', 'Reallocation · native ~500 m')]:
            observed = restore(np.where(a['common'], a[source], np.nan))
            residuals = {m: restore(a[f'{source}_{m}_residual']) for m in ('M0', 'M1')}
            scenes[key] = Scene(key, title, observed, median(observed[:cut], axis=0), residuals,
                                mask, transform, 4326, dates, dates[cut:])
    fine = base / 'fine'
    fine_manifest = json.loads((fine / 'manifest.json').read_text())
    models = root / 'data/published/counterfactual_reallocation_10m'
    model_manifest = json.loads((models / 'manifest.json').read_text())
    if fine_manifest['status'] != 'complete' or model_manifest['status'] != 'complete':
        raise ValueError('Fine input or model run incomplete')
    if fine_manifest['identity']['comparison_manifest_sha256'] != file_hash(base / 'manifest.json'):
        raise ValueError('Fine/native provenance mismatch')
    if model_manifest['identity']['source_manifest_sha256'] != file_hash(fine / 'manifest.json'):
        raise ValueError('Fine model/input provenance mismatch')
    for name, sha in fine_manifest['output_sha256'].items():
        verify(fine / name, sha)
    fine_residuals = {}
    for method in ('M0', 'M1'):
        name = f'{method}/residual.npy'
        verify(models / name, model_manifest['output_sha256'][name])
        fine_residuals[method] = np.load(models / name, mmap_mode='r')
    verify(models / 'baseline_median.tif', model_manifest['output_sha256']['baseline_median.tif'])
    with rasterio.open(models / 'baseline_median.tif') as src:
        baseline, fine_transform = src.read(1), src.transform
    scenes['allocated10'] = Scene('allocated10', 'Reallocation · 10 m proxy',
        np.load(fine / 'reallocation_10m.npy', mmap_mode='r'), baseline, fine_residuals,
        np.load(fine / 'native_aoi_mask.npy', mmap_mode='r'), fine_transform, 32616,
        pd.DatetimeIndex(pd.to_datetime(fine_manifest['dates'])),
        pd.DatetimeIndex(pd.to_datetime(model_manifest['event_dates'])))
    assert all(s.dates.equals(dates) and s.event_dates.equals(dates[cut:]) for s in scenes.values())
    return scenes, used


def paint(ax, scene, array, bounds=None, *, title='', cmap='magma', vmin=0, vmax=20,
          norm=None, max_side=800, missing_label=True):
    _, _, _, geometry = scene.region(bounds)
    colors = plt.get_cmap(cmap).copy() if isinstance(cmap, str) else cmap.copy()
    colors.set_bad('#dddddd')
    kw = dict(cmap=colors, norm=norm) if norm is not None else dict(cmap=colors, vmin=vmin, vmax=vmax)
    if scene.crs == 32616:
        xe, ye = geometry
        stride = max(1, int(np.ceil(max(array.shape) / max_side)))
        image = ax.imshow(array[::stride, ::stride], extent=(xe[0]/1000, xe[-1]/1000, ye[-1]/1000, ye[0]/1000),
                          origin='upper', interpolation='nearest', **kw)
    else:
        x, y = geometry
        image = ax.pcolormesh(x/1000, y/1000, array, shading='flat', rasterized=True, **kw)
    if bounds is not None:
        ax.set_xlim(bounds[0]/1000, bounds[2]/1000)
        ax.set_ylim(bounds[1]/1000, bounds[3]/1000)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    if missing_label and not np.isfinite(array).any():
        ax.text(.5, .5, 'No valid observations', ha='center', va='center', transform=ax.transAxes,
                fontsize=9, color='#444444', bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
    return image


def image_normalizations(scenes, zoom, image_dates, *, close=None, close_dates=(),
                         radiance_vmax=None, residual_limit=None, linthresh=.5):
    """Shared full-data map ranges, log1p radiance and symmetric-log residuals.

    Scan actual displayed values (before display sampling), ignoring NaN. Manual
    bounds remain optional; auto bounds include every displayed finite value.
    Color normalization never transforms the input arrays or statistical inputs.
    """
    radiance_max, residual_max = 0., 0.

    def maximum(values, signed=False):
        finite = values[np.isfinite(values)]
        if not finite.size:
            return 0.
        if not signed and np.any(finite < 0):
            raise ValueError('log1p radiance maps expect nonnegative radiance')
        return float(np.max(np.abs(finite))) if signed else float(np.max(finite))

    for scene in scenes.values():
        radiance_max = max(radiance_max, maximum(scene.values('baseline')))
        for day in image_dates:
            radiance_max = max(radiance_max, maximum(scene.values('observed', zoom, day)))
            for method in ('M0', 'M1'):
                residual_max = max(residual_max, maximum(scene.values('residual', zoom, day, method), True))
        if close is not None:
            for day in close_dates:
                for method in ('M0', 'M1'):
                    residual_max = max(residual_max, maximum(scene.values('residual', close, day, method), True))
    top = max(radiance_max, 1.) if radiance_vmax is None else float(radiance_vmax)
    limit = max(residual_max, linthresh) if residual_limit is None else float(residual_limit)
    if not (np.isfinite(top) and top > 0 and np.isfinite(limit) and limit > 0
            and np.isfinite(linthresh) and linthresh > 0):
        raise ValueError('Color limits and the symmetric-log linear threshold must be positive and finite')
    return (FuncNorm((np.log1p, np.expm1), vmin=0, vmax=top, clip=False),
            SymLogNorm(linthresh=linthresh, linscale=1., base=10, vmin=-limit, vmax=limit, clip=False))


def map_colorbar(fig, image, *, label, **kwargs):
    """Keep colorbar labels in original radiance units on nonlinear scales."""
    if isinstance(image.norm, FuncNorm):
        high = image.norm.vmax
        ticks = [0.] + [10.**p for p in range(int(np.floor(np.log10(high))) + 1)]
        if high < 1:
            ticks = [0., high / 2, high]
        kwargs['ticks'] = ticks
        label += ' · log1p colors'
    elif isinstance(image.norm, SymLogNorm):
        # Omit sub-threshold decade ticks: on short horizontal bars they crowd
        # the zero label even though that central band is linear.
        start = int(np.ceil(np.log10(image.norm.linthresh)))
        stop = int(np.floor(np.log10(image.norm.vmax)))
        positive = 10. ** np.arange(start, stop + 1)
        kwargs['ticks'] = np.r_[-positive[::-1], 0., positive]
        label += ' · symmetric-log colors'
    return fig.colorbar(image, label=label, **kwargs)


def reference_figure(scene, zoom, vmax, *, norm=None):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), layout='constrained')
    im = paint(axes[0], scene, scene.values('baseline'), title='Dyer analysis footprint', vmax=vmax, norm=norm)
    paint(axes[1], scene, scene.values('baseline', zoom), zoom, title='Built-up areas of interest: Dyersburg and Newbern', vmax=vmax, norm=norm)
    axes[0].add_patch(Rectangle((zoom[0]/1000, zoom[1]/1000), (zoom[2]-zoom[0])/1000,
                               (zoom[3]-zoom[1])/1000, fill=False, edgecolor='#00c8ff', lw=1.5))
    for name, lon, lat in [('Dyersburg', -89.386, 36.035), ('Newbern', -89.262, 36.113)]:
        x, y = TO_UTM.transform(lon, lat)
        axes[1].annotate(name, (x/1000, y/1000), xytext=(3, 9), textcoords='offset points',
                         color='white', fontsize=9, bbox=dict(facecolor='black', alpha=.55, edgecolor='none'))
    map_colorbar(fig, im, ax=axes, shrink=.75, label=UNITS)
    fig.suptitle(f'{scene.title}\nPer-cell baseline median · February 5–April 1', fontsize=13)
    return fig


def observed_montage(scene, zoom, dates, vmax, *, norm=None):
    ncols = 5
    fig, axes = plt.subplots(int(np.ceil(len(dates)/ncols)), ncols, figsize=(15, 7.2), layout='constrained', squeeze=False)
    for ax, day in zip(axes.flat, dates):
        im = paint(ax, scene, scene.values('observed', zoom, day), zoom, title=day.strftime('%b %d'), vmax=vmax, norm=norm)
    for ax in list(axes.flat)[len(dates):]:
        ax.set_visible(False)
    map_colorbar(fig, im, ax=axes.ravel().tolist(), shrink=.65, label=UNITS)
    fig.suptitle(f'{scene.title} · observed radiance\nAll 10 calendar nights, April 2–11', fontsize=13)
    return fig


def comparison_panel(scene, method, zoom, dates, vmax, residual_limit, *, radiance_norm=None, residual_norm=None):
    fig, axes = plt.subplots(len(dates), 3, figsize=(11, 2.55*len(dates)), layout='constrained')
    baseline = scene.values('baseline', zoom)
    for row, day in enumerate(dates):
        im = paint(axes[row, 0], scene, baseline, zoom,
                   title=f'{day:%b %d} · Baseline median', vmax=vmax, norm=radiance_norm)
        paint(axes[row, 1], scene, scene.values('observed', zoom, day), zoom,
              title='Observed', vmax=vmax, norm=radiance_norm)
        residual_im = paint(axes[row, 2], scene, scene.values('residual', zoom, day, method), zoom,
              title=f'{method} residual', cmap='RdBu_r', vmin=-residual_limit, vmax=residual_limit, norm=residual_norm)
    map_colorbar(fig, im, ax=axes[:, :2].ravel().tolist(), orientation='horizontal', shrink=.65,
                 fraction=.015, label=f'Radiance ({UNITS})')
    map_colorbar(fig, residual_im, ax=axes[:, 2].ravel().tolist(), orientation='horizontal',
                 fraction=.015, label=f'Observed − {method} prediction ({UNITS})')
    fig.suptitle(f'{scene.title} · {method} · built-up-area zoom', fontsize=14)
    return fig


def radiance_series(scene, dates, bounds=None):
    key = (None if bounds is None else tuple(bounds), tuple(pd.DatetimeIndex(dates).astype(str)))
    if key in scene._series:
        return scene._series[key]
    baseline_days = scene.dates[scene.dates < scene.event_dates[0]]
    baseline_nightly = np.array([median(scene.values('observed', bounds, day)) for day in baseline_days])
    table = pd.DataFrame({'date': dates,
        'median_radiance': [median(scene.values('observed', bounds, day)) for day in dates],
        'valid_pixels': [int(np.isfinite(scene.values('observed', bounds, day)).sum()) for day in dates]})
    result = table, float(median(baseline_nightly))
    scene._series[key] = result
    return result


def format_time(ax, dates, ylabel):
    ax.set_xlim(dates[0], dates[-1])
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.22)


def radiance_figure(scene, dates, bounds=None, area='full Dyer scene'):
    table, reference = radiance_series(scene, dates, bounds)
    fig, ax = plt.subplots(figsize=(11, 3.7), layout='constrained')
    ax.plot(table.date, table.median_radiance, 'o-', color='#0072B2', label='Nightly spatial median')
    ax.axhline(reference, color='#cc6600', ls='--', label='Median of baseline nightly spatial medians')
    format_time(ax, dates, f'Radiance ({UNITS})')
    ax.set_title(f'{scene.title} · {area}')
    ax.legend(fontsize=9)
    return fig, table.assign(baseline_reference=reference)


def residual_figure(scene, method, dates, zoom):
    table = pd.DataFrame({'date': dates,
        'median_residual': [median(scene.values('residual', zoom, day, method)) for day in dates],
        'valid_pixels': [int(np.isfinite(scene.values('residual', zoom, day, method)).sum()) for day in dates]})
    fig, ax = plt.subplots(figsize=(11, 3.7), layout='constrained')
    ax.plot(table.date, table.median_residual, 'o-', color='#6a3d9a')
    ax.axhline(0, color='black', lw=1, ls='--')
    format_time(ax, dates, f'Residual ({UNITS})')
    ax.set_title(f'{scene.title} · {method} · median residual in the built-up-area zoom')
    return fig, table


def missing_kmeans(values, k=4, min_dates=4, max_iter=40, seed=7, chunk=100000):
    """Lloyd clustering of observed entries, normalizing each pixel by its count.

    Minimize sum_i [sum_{t observed} (y_it - c_label,t)^2 / n_i]. Missing entries
    have zero weight, never act as zero-valued radiance, and are not interpolated.
    All eligible pixels enter every fitting iteration; only initialization uses
    a bounded reproducible sample. Centroid coordinates without observations in
    an iteration retain their previous value (they add no term to that update).
    """
    y = np.asarray(values, dtype='float32')
    if y.ndim != 2 or np.isinf(y).any():
        raise ValueError('Expected pixel × date residuals with NaN missingness')
    observed_dates = np.isfinite(y).any(axis=0)
    eligible = np.isfinite(y).sum(axis=1) >= min_dates
    labels = np.full(len(y), -1, dtype='int16')
    if not eligible.any() or observed_dates.sum() < min_dates:
        return labels, {'eligible_pixels': 0, 'clusters': 0, 'iterations': 0, 'converged': True}
    indices = np.flatnonzero(eligible)
    x = np.ascontiguousarray(y[eligible][:, observed_dates])
    k = min(k, len(x))
    rng = np.random.default_rng(seed)
    sample = x[rng.choice(len(x), min(20000, len(x)), replace=False)]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered')
        centers = np.nanquantile(sample, np.linspace(.1, .9, k), axis=0)
    for j in range(x.shape[1]):
        if not np.isfinite(centers[:, j]).all():
            centers[:, j] = median(x[:, j])
    centers = centers.astype(float)
    assigned = np.zeros(len(x), dtype='int16')
    converged = False

    def assign(block):
        valid = np.isfinite(block)
        z = np.where(valid, block, 0).astype(float)
        count = valid.sum(axis=1)
        distances = (np.sum(z*z, axis=1)[:, None] - 2*z@centers.T + valid.astype(float)@(centers*centers).T)
        distances /= count[:, None]
        return np.argmin(distances, axis=1), z, valid, count

    for iteration in range(max_iter):
        numerator = np.zeros_like(centers)
        denominator = np.zeros_like(centers)
        changes = 0
        for start in range(0, len(x), chunk):
            end = min(start + chunk, len(x))
            lab, z, valid, count = assign(x[start:end])
            changes += int(np.count_nonzero(lab != assigned[start:end]))
            assigned[start:end] = lab
            weights = 1. / count
            for j in range(x.shape[1]):
                numerator[:, j] += np.bincount(lab, weights=z[:, j]*weights, minlength=k)
                denominator[:, j] += np.bincount(lab, weights=valid[:, j]*weights, minlength=k)
        updated = np.divide(numerator, denominator, out=centers.copy(), where=denominator > 0)
        shift = np.linalg.norm(updated-centers) / max(1., np.linalg.norm(centers))
        centers = updated
        if changes == 0 and iteration > 0 or shift < 1e-5:
            converged = True
            break
    # Assign once more against the final centroids, then order groups by mean residual.
    for start in range(0, len(x), chunk):
        assigned[start:start + chunk] = assign(x[start:start + chunk])[0]
    present = np.unique(assigned)
    order = present[np.argsort(np.mean(centers[present], axis=1))]
    remap = np.full(k, -1, dtype='int16')
    remap[order] = np.arange(len(order), dtype='int16')
    labels[indices] = remap[assigned]
    return labels, {'eligible_pixels': int(eligible.sum()), 'clusters': len(order),
                    'iterations': iteration + 1, 'converged': converged,
                    'minimum_observed_dates': min_dates, 'seed': seed,
                    'used_date_columns': np.flatnonzero(observed_dates).tolist()}


def cluster_scene(scene, method, zoom, dates, k=4, min_dates=4, max_iter=40, seed=7):
    _, _, footprint, _ = scene.region(zoom)
    rows = [scene.values('residual', zoom, day, method)[footprint] for day in dates]
    values = np.stack(rows, axis=1)
    del rows
    labels, info = missing_kmeans(values, k=k, min_dates=min_dates, max_iter=max_iter, seed=seed)
    image = np.full(footprint.shape, np.nan, dtype='float32')
    image[footprint] = np.where(labels >= 0, labels + 1, np.nan)
    records = []
    for group in range(info['clusters']):
        subset = values[labels == group]
        medians = median(subset, axis=0)
        counts = np.isfinite(subset).sum(axis=0)
        for day, value, n in zip(dates, medians, counts):
            records.append({'cluster': group + 1, 'date': day, 'median_residual': value,
                            'observed_pixels': int(n), 'cluster_pixels': len(subset), 'color': COLORS[group]})
    return image, pd.DataFrame(records), info


def cluster_figures(scene, method, zoom, display_dates, image, table, info):
    fig, ax = plt.subplots(figsize=(11, 4), layout='constrained')
    for group in range(1, info['clusters'] + 1):
        part = table[table.cluster == group]
        size = int(part.cluster_pixels.iloc[0])
        ax.plot(part.date, part.median_residual, 'o-', color=COLORS[group-1], label=f'Group {group} (n={size:,})')
    ax.axhline(0, ls='--', color='black', lw=1)
    format_time(ax, display_dates, f'Median residual ({UNITS})')
    ax.set_title(f'{scene.title} · {method} · residual profiles by cluster')
    if info['clusters']:
        ax.legend(ncol=info['clusters'], fontsize=9)
    spatial, sax = plt.subplots(figsize=(8, 6), layout='constrained')
    k = max(info['clusters'], 1)
    colors = ListedColormap(COLORS[:k])
    im = paint(sax, scene, image, zoom, title=f'{scene.title} · {method}\nClusters of residual time series',
               cmap=colors, norm=BoundaryNorm(np.arange(.5, k + 1.5), k), missing_label=False)
    if info['clusters']:
        cb = spatial.colorbar(im, ax=sax, ticks=np.arange(1, k+1), shrink=.65)
        cb.ax.set_yticklabels([f'Group {i}' for i in range(1, k+1)])
    threshold = info.get('minimum_observed_dates', 4)
    sax.set_xlabel(f'Gray: outside footprint or fewer than {threshold} observed residual dates.\nGroups are not outage labels.',
                   fontsize=8)
    return fig, spatial


def close_comparison(a2, fine, method, bounds, dates, residual_limit, *, norm=None):
    fig, axes = plt.subplots(len(dates), 2, figsize=(11, 4*len(dates)), layout='constrained')
    for row, day in enumerate(dates):
        for col, scene in enumerate((a2, fine)):
            image = paint(axes[row, col], scene, scene.values('residual', bounds, day, method), bounds,
                          title=f'{day:%b %d} · {"A2 ~500 m" if col == 0 else "Reallocation 10 m"} · {method}',
                          cmap='RdBu_r', vmin=-residual_limit, vmax=residual_limit, max_side=1200, norm=norm)
    map_colorbar(fig, image, ax=axes.ravel().tolist(), shrink=.7, label=f'Residual ({UNITS})')
    fig.suptitle(f'Dyersburg close-up · {method}\nIdentical map bounds and color limits; each method keeps its own valid support', fontsize=13)
    return fig
