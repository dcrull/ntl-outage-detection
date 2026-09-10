# M0 and regularized M1: implementation and first results

Implemented 2026-09-09 from the [counterfactual plan](COUNTERFACTUAL_DETECTOR_PLAN.md). This first version predicts native-cell normal radiance and reports deviations. It does not classify outages, infer outage duration, or implement reallocation. No TING data are read.

## Run and inspect

From the repository root, use the project environment:

```bash
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python scripts/run_counterfactual.py
```

The run is offline and reads the existing `data/viirs/VNP46A2.002/` cache. Open [03_inspect_counterfactual.ipynb](../notebooks/03_inspect_counterfactual.ipynb) with the project's `.venv` kernel. It reads the completed outputs, shows validation, and allows individual-cell inspection by longitude/latitude.

Outputs are in `data/published/counterfactual_dyer/`, visible on the Mac at **`~/mnt/ntl/published/counterfactual_dyer/`**. Rerunning the same configuration replaces generated files; close and reload QGIS layers afterward. A changed configuration requires a separate `--output` directory so experiments cannot silently mix.

For example, a later sensitivity run can change the permitted level movement:

```bash
.venv/bin/python scripts/run_counterfactual.py \
  --process-ratio 0 \
  --output data/published/counterfactual_dyer_static_m1
```

That example has not been selected as the primary model. The current defaults were fixed before examining held-out scores or event maps.

## Inputs and spatial support

- Only `radiance.tif`, tagged as **DNB_BRDF-Corrected_NTL**, is used for radiance. The reader does not open the gap-filled layer.
- Reuses the ingestion's strict usability rule and independently recomputes it from cloud, snow, mandatory QA, and radiance validity. A disagreement fails explicitly. Valid zero radiance remains an observation.
- Checks consumed raster hashes against their extraction metadata, grid consistency across dates, source-date correspondence, and native tile alignment.
- Retains every calendar date, including missing source dates and cloudy nights. No interpolation is performed.
- Uses 7,868 cell centers inside Dyer. IDs are `h09v05_rRRRR_cCCCC`, using original tile rows and columns. No reprojection or resampling is applied to outputs.
- Fits February 5–April 1, then forecasts April 2–17. Product dates are acquisition-day labels; matching to local event times remains unresolved.

## Models and fixed regularization

Both models require at least **10 usable training observations**. This provisional support rule gives 7,688 supported cells in the final fit. Dark cells remain present and identifiable; support is not a claim that outages can be detected in those cells.

Scale pooling uses quartiles of training median brightness among supported cells. Each group's reference scale is the median of its cells' `1.4826 × MAD`, with a **0.1 nW cm⁻² sr⁻¹ minimum scale**. The floor is a provisional regularization choice, not a NASA instrument specification or a detection threshold. It needs sensitivity assessment before classification. All grouping and scale estimation are repeated inside each training fold; no holdout/event values enter them.

### M0: constant median reference

The central prediction is the cell's training median. Its variance estimate blends squared robust cell scale and squared group scale with weights `n/(n+20)` and `20/(n+20)`, then applies the scale floor. This avoids zero-width or unstable intervals for sparse, nearly constant histories.

The predictive reference is Student-t with `n−1` degrees of freedom and scale:

`sqrt(shrunk_variance × (1 + pi/(2n)))`

The second factor approximates uncertainty in a normal-sample median. This is an **approximate predictive reference**, not an exact Bayesian posterior for the median or an empirically calibrated interval.

### M1: local level with strongly limited movement

The observation model is `y_t = level_t + error_t`, with `error_t ~ Normal(0, theta)`. The level evolves by `level_t = level_(t−1) + eta_t`, with `eta_t ~ Normal(0, q × theta)`.

- **q = 0.001 per calendar day**, fixed across cells. The daily level-innovation standard deviation is approximately 3.2% of the observation-noise standard deviation.
- No trend, seasonal terms, external controls, or cell-specific tuning of q.
- Observation variance `theta` has an inverse-gamma prior with shape `20/2` and scale `20 × group_scale²/2`. These 20 prior degrees of freedom shrink variance estimation toward the training brightness group.
- A diffuse location prior is conditioned on the first valid observation, initializing the level to that observation and its conditional variance to `theta`.
- Subsequent valid observations use the normal/inverse-gamma filter updates. Missing dates advance state uncertainty without a measurement update. This version does not clip large innovations or implement a robust outlier-mixture observation model.
- After April 1, the target-cell filter receives **no further observations**. The predicted level stays constant because there is no drift term, while state uncertainty increases with elapsed forecast days.

The Student-t posterior predictive distribution includes conditional level uncertainty, observation noise, and uncertainty in `theta`. It is conditional on the fixed q and empirical-Bayes pooled prior; uncertainty in those hyperparameters is not integrated. The model is therefore not presented as a fully hierarchical Bayesian analysis or a calibrated causal estimate.

M0 uses a median, whereas M1 uses a filtered mean. Differences between them do not isolate the effect of temporal dynamics alone; the zero-process-variance M1 option provides a possible later ablation.

## Validation and initial findings

Four seven-day folds train on the first 28, 35, 42, and 49 baseline calendar days. Two fourteen-day folds train on the first 28 and 42 days. Each forecast is frozen for its entire block, which tests prediction across contiguous withholding intervals. Naturally cloudy observations remain missing and cannot be scored.

Within each design the test dates do not overlap. The two designs reuse dates, so report them separately. Both models use identical eligibility rules and scored observations. Output tables include support and sample counts because early training folds can exclude many cells.

Initial results below pool scored cell-nights; radiance errors are in nW cm⁻² sr⁻¹. They are descriptive, without standard errors that would assume neighboring cells are independent.

| Forecast design | Model | MAE | RMSE | Nominal 90% interval coverage |
| --- | --- | ---: | ---: | ---: |
| 7 days | M0 | 0.7701 | 2.9377 | 87.0% |
| 7 days | M1 | 0.7805 | 2.8674 | 88.9% |
| 14 days | M0 | 0.7808 | 2.9449 | 86.5% |
| 14 days | M1 | 0.7936 | 2.8889 | 88.0% |

M1 modestly improves RMSE and overall interval coverage; M0 has slightly lower absolute error. This is not a decisive model-selection result. Importantly, for cells with training median radiance ≥5, seven-day interval coverage is only **72.1% for M0 and 77.2% for M1**. Their lower-bound exceedance rates are approximately 10–11%, versus a nominal 5%. County-wide averages conceal this problem because dim cells dominate.

The intervals are not recalibrated using these same held-out outcomes. Additional calibration would need a separate or nested temporal evaluation. Scale-floor/prior sensitivity, bright-cell residual inspection, and potentially more pre-event history should precede converting tail scores into detection rules. Viewing geometry remains a possible later explanation to investigate, not an established cause of these errors.

April 2–6 have no evaluable satellite observations. Forecast rasters exist, but residuals and tail scores are blank. These nights cannot support an observed outage or normal-power conclusion.

## Output reference

| Output | Meaning |
| --- | --- |
| `manifest.json` | Run status, settings, versions, code hashes, consumed-input hashes, and output hashes |
| `model_details.json` | Windows, support, scale groups, regularization, and interval caveats |
| `cells.csv` | Stable native IDs, coordinates, support, robust scales, and fitted references |
| `source_dates.csv` | Acquisition dates, source availability, and source filenames |
| `analysis_arrays.npz` | Date × county-cell inputs/QA plus event predictions, grid mask, indices, and transform; missing observations are NaN |
| `validation_summary.csv`, `validation_daily.csv` | Fold, horizon, and training-brightness-stratum scores, including prediction error, interval coverage, tail rates, interval score, and log score |
| `validation_predictions.npz` | Saved fold-specific predictive distributions for reproducible inspection |
| `event_cell_nights.csv.gz` | Both methods, every event cell-night, QA, support/status, prediction, interval, residual, and lower-tail score |
| `event_daily_summary.csv` | Evaluable counts and descriptive residual/tail summaries; proportions refer to evaluable cells, not county area |
| `baseline_*.tif` | Native-grid median, training counts, and support mask |
| `observed/YYYY-MM-DD.tif` | QA-screened, non-gap-filled observed radiance |
| `M0/daily/YYYY-MM-DD/`, `M1/daily/YYYY-MM-DD/` | `center`, `q05`, `q95`, `residual`, `standardized_residual`, and `lower_tail` COGs |

`center` is also the predictive distribution's median. `q05/q95` form a nominal central 90% predictive interval. The unbounded statistical approximation can yield negative lower bounds; these are retained rather than silently clipped. `residual = observed − center`. Standardized residual divides by predictive standard deviation. `lower_tail` is the reference-model CDF at observed radiance, **not outage probability**.

The table separately retains support, usability, and source availability. Status distinguishes an observed supported cell, insufficient baseline, and no usable observation. Insufficient baseline takes precedence in the text status; the other columns preserve simultaneous missingness. Blank scores never mean normal power. There are no detected-outage polygons, onset/recovery estimates, or TING-derived labels in these outputs.

## Verification

Run `.venv/bin/python -m unittest discover -s tests -v`. Tests cover the analytic zero-process-variance posterior, robust M0 behavior, valid zeros, unsupported cells, calendar-gap uncertainty, absence of scores on missing observations, held-out/event leakage, validation fold boundaries, and simulated interval coverage. A synthetic cache test checks missing-source dates, independent QA screening, native IDs, and rejection of a gap-filled source masquerading as the radiance layer.

The actual pipeline has been run on the full cached window, and the inspection notebook has been executed using the project environment. The first results above are intended for review before extending the model or building the reallocation version.
