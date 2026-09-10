# M0/M1 on the 10 m reallocated proxy

Run the existing model specifications on the retained nightly structural allocation:

```bash
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python scripts/run_counterfactual_fine.py
```

Inputs are the verified `data/published/reallocation_check/fine/` cube and mask.
The allocation is not recomputed, and neither TING nor gap-filled A2 enters this
run. Output is `data/published/counterfactual_reallocation_10m/`. The manifest is
marked complete only after validation, model fitting, raster export, and output
hashing finish. An existing complete, unchanged run is verified and reused;
an incomplete or changed output requires a separate `--output` directory.

The full run completed on September 10, 2026. Of 13,671,278 fine pixels in the
analysis footprint, 12,575,023 meet the baseline observation requirement. Outputs
include 192 daily model COGs, three baseline COGs, 12 memory-mapped model arrays,
validation and daily support tables, and provenance (about 25.3 GB total).
April 2–6 have no usable observations and therefore no residual/tail evidence.
April 15 has only 4,311 evaluable fine pixels; the fine grid does not remove the
storm's observation gaps.

Verification passed: 17 relevant model/reallocation tests; all 210 output hashes;
all 192 daily COG grids and distributed data windows; exact saved predictions
and scores at 2,048 sampled AOI pixels; and full-array missingness checks for
every field. The audit is saved as
`data/reallocation/counterfactual_10m_validation.json`. It can be rerun with
`.venv/bin/python scripts/verify_counterfactual_fine.py`.

## Model and observation choices

- Baseline: February 5–April 1, 2025 (56 calendar nights).
- Open-loop event forecasts: April 2–17 (16 nights). No event observations update
  the normal-state forecast.
- M0: baseline median with MAD variance shrunk toward the training brightness
  group's pooled scale, plus approximate median uncertainty.
- M1: the existing local-level model with daily process/observation variance ratio
  0.001 and inverse-gamma scale prior weight 20; no trend.
- Both require at least 10 valid training observations. The original scale floor
  of 0.1 nW cm^-2 sr^-1 and nominal central 90% reference intervals are retained.
- The AOI contains fine pixels with positive area overlap with the selected native
  Dyer analysis cells. The surrounding processing margin is excluded from pooling,
  fitting, and published model fields.
- Both methods use the same fine-grid observations. Fine pixels retain their own
  complete-kernel support; the additional native complete-footprint aggregation
  mask is not imposed. Fine and native fits therefore need not use identical
  nights or be numerically equivalent.

The cube is read in spatial chunks (`--chunk-rows`, default 32). Each training
origin first computes brightness quartiles and pooled MAD scales over the entire
eligible AOI, using only observations before that origin. Those common priors
are passed to the existing model functions for every chunk. Changing chunk sizes
does not change the model. Repeated Student-t quantiles are evaluated once per
unique degree of freedom; tests verify equality with the original dense formula.

## Validation

The same pre-event validation designs are retained: seven-day forecasts after
28, 35, 42, and 49 training days, and fourteen-day forecasts after 28 and 42 days.
The two designs reuse some target dates, so examine them separately.

`validation_summary.csv` contains per-fold M0/M1 error, interval coverage, and
log-density scores for all eligible pixels and the original brightness strata
(less than 1, 1–5, and at least 5). Additive scores are accumulated across chunks;
RMSE is calculated from pooled squared errors, not averaged chunk RMSEs.

Fine pixels share both native radiance and structural inputs. Counts in these
tables are descriptive pixel-night counts, not independent validation samples.
Uncertainty coverage must be evaluated anew at this scale; native calibration
does not transfer automatically.

The completed baseline validation gives the following scores over all eligible
pixel-nights, weighted by observation count across folds within each design:

| Design | Model | MAE | RMSE | Nominal 90% interval coverage |
| --- | --- | ---: | ---: | ---: |
| 7 days | M0 | 0.5230 | 2.2052 | 82.71% |
| 7 days | M1 | 0.5157 | 2.1262 | 85.02% |
| 14 days | M0 | 0.5150 | 2.0706 | 81.04% |
| 14 days | M1 | 0.5023 | 1.9575 | 83.55% |

MAE and RMSE use radiance units (nW cm^-2 sr^-1). These are errors predicting the
fine structural proxy, not errors against an independent outage benchmark.
Coverage remains below the nominal 90%, so the uncertainty references still need
calibration. M1's lower errors here do not establish better outage detection.

## Outputs and inspection

Each method has the following files:

```text
M0/                         # M1 has the same structure
  center.npy
  q05.npy
  q95.npy
  residual.npy
  standardized_residual.npy
  lower_tail.npy
  daily/YYYY-MM-DD/
    center.tif
    q05.tif
    q95.tif
    residual.tif
    standardized_residual.tif
    lower_tail.tif
```

The NumPy arrays have shape `(16, 4097, 5533)` with axes `(event night, row, column)`.
Their event dates are in `manifest.json`; index zero is April 2, not February 5.
They preserve the full input array layout but contain NaN outside the analysis
footprint and where baseline support is insufficient. Forecasts can exist on
cloudy dates; residuals, standardized residuals, and tails require an observation.

To inspect one model without reading the full array into RAM:

```python
import numpy as np
from pathlib import Path

model_dir = ROOT / "data/published/counterfactual_reallocation_10m"
m0_residual = np.load(model_dir / "M0/residual.npy", mmap_mode="r")
m1_residual = np.load(model_dir / "M1/residual.npy", mmap_mode="r")
# m0_residual[6] is April 8: one full height × width scene.
```

The daily rasters are 10 m EPSG:32616 COGs for QGIS. On the Mac, an example is:

```text
~/mnt/ntl/published/counterfactual_reallocation_10m/M0/daily/2025-04-08/residual.tif
```

Use a diverging, zero-centered color scale for residuals. Negative values mean
observed radiance is below its model prediction. The baseline support/count/median
rasters and `event_daily_summary.csv` help distinguish unavailable evidence from
low radiance.

These outputs are diagnostics of a structural proxy, not independently resolved
10 m outages or causal attribution. The same contemporary 2026 structural data
underlie the 2025 allocation. `lower_tail` is an uncalibrated predictive CDF, not
an outage probability. No outage-state inference or TING validation is added here.
