# Reallocation at native VIIRS resolution

This experiment reproduces the original structural allocation, aggregates its output to the actual VNP46A2 grid, and compares M0/M1 on matched observations. Fine-scale outage analysis is deferred. Internal 10 m arrays are numerical working inputs, not validated 10 m radiance or outage products.

## Original implementation and the conservation question

The source checkout is `~/research/ntl-psf-disaggregation`, corresponding to [dcrull/ntl-psf-disaggregation](https://github.com/dcrull/ntl-psf-disaggregation). The inspected commit is `73ce981308ff40de8b66db3a97498f899cec7714`. The reproduction records hashes of the actual source files and final-results document as well as the commit, because some explanatory documents are untracked in that checkout. The source repository is read without modification.

The original primary operator is:

`allocated(x) = h(x) × mean_kernel(L)(x) / max(mean_kernel(h)(x), epsilon)`

Here `h` is a normalized structural proxy and the kernel is a circular mean with **500 m radius**. This neighborhood spans multiple native cells. It is distinct from normalizing structural weights separately inside each native VIIRS footprint. The upstream operator explicitly records `exact_conservation_assumed: false`.

Therefore equality after aggregation is a hypothesis to test, not an algebraic property of this implementation. Even the uniform-proxy case averages neighboring radiance, changing a nonconstant coarse field. The check does not rescale the allocated field to force agreement with its source.

## Reproduced choices

The original final-results document selects built-form allocation, no water prior, and the circular reference kernel. It supersedes the older YAML's combined-water reporting choice. We follow that final structural configuration, while retaining this project's strict nightly A2 screening rather than adopting the original demonstration's broad-QA temporal composite.

- Overture buildings and transportation segments: pinned release **2026-07-22.0**.
- Structural formula: `0.7 × sqrt(building_fraction) + 0.3 × normalized_road_density`.
- Building fraction: union occupancy sampled at 2 m, averaged into the internal 10 m grid.
- Roads: original class weights, segmentized midpoint length allocation with maximum 1 m segments; density saturates at 20,000 weighted metres per km².
- Positive proxy floor: `h_raw = 0.05 + 0.95 × structural_formula`, followed by mean-one normalization.
- No water prior; denominator epsilon is `1e-6 × normalized_proxy_mean`.
- Internal grid: 10 m UTM zone 16N, covering the native county crop plus an approximately 1,100 m projected halo. The raw VIIRS source crop is derived from this entire envelope, including its rotated corners.
- Radiance: nearest-cell assignment of each night's **non-gap-filled A2 corrected radiance**, using the existing strict quality rule. Raw HDF5 halo values are checked against cached county values.
- Allocation requires complete valid kernel support, within a `1e-6` fractional tolerance. Cloudy neighboring cells can therefore make an otherwise clear target unavailable.

The Overture release is contemporary relative to the April 2025 event. This is a documented reproduction of the original structural assumptions, **not an event-date inventory of buildings and roads**. It must not be used to claim historical fine-scale validation.

The wrapper imports the original building/road rasterizers, kernel builder, and numerical convolution. Its nightly calculation caches the static structural gain and computes only the changing radiance component. A test compares this optimized path with the original complete allocation function, including cloudy neighborhoods and valid zero observations.

## Native aggregation and controls

All published rasters use the original VIIRS transform and county-center mask. Native identifiers retain source tile rows and columns.

Each native cell polygon is densified in geographic coordinates and projected to the working CRS. Its intersections with 10 m working pixels define an area-weighted mean, preserving radiance units. Overlap weights must sum to one; an incompletely covered footprint fails instead of being silently normalized. A daily aggregate is retained only when effectively its whole footprint has valid fine-grid values.

Four inputs are retained:

| Input | Purpose |
| --- | --- |
| `direct` | Original strict-screened native A2 radiance |
| `upsample` | Direct nearest upsampling followed by area aggregation; measures projected-grid round-trip error |
| `uniform` | Kernel-smoothed radiance with uniform structural gain; isolates smoothing |
| `allocated` | Original built-form allocation aggregated by native-cell area |

The UTM 10 m grid does not nest inside the geographic VIIRS grid. Fine pixels crossing a native boundary can introduce round-trip differences even in `upsample`. This control is essential when interpreting allocation differences. Aggregation is a mean, not a sum of radiance values interpreted as physical flux.

Daily comparison tables report common support, bias, MAE, RMSE, maximum absolute difference, area-weighted bias, and the fraction satisfying `abs(proxy − direct) <= 1e-5 + 1e-5 × abs(direct)`. The tolerance is a numerical equality diagnostic, not a scientifically established acceptable-error threshold.

## M0/M1 comparison

For both fitting and evaluation, all four variants use the **intersection of their valid cell-night masks**. M0/M1 settings and pre-event/event windows are unchanged. Direct A2 is refit on this common mask, so model differences can be attributed to input values rather than differing training dates.

This matched direct fit is not interchangeable with the earlier direct-only run, which had greater observation coverage. The new support table explicitly reports that loss. Pre-event seven-day and fourteen-day validation designs remain separate; TING never enters fitting, support selection, or evaluation here.

The pipeline reports differences in model references, predictive intervals, residuals, and lower-tail scores. Those scores retain the calibration limitations of the initial M0/M1 experiment and are not outage probabilities.

## Reproduction

Use the project environment:

```bash
.venv/bin/python -m pip install -r requirements-reallocation.txt
```

Retrieve the pinned structural data if absent. The envelope covers the complete processing halo. The [official Overture Python client](https://docs.overturemaps.org/getting-data/overturemaps-py/) supports bounded GeoParquet retrieval and an explicit release.

```bash
mkdir -p data/reallocation/overture/2026-07-22.0
.venv/bin/overturemaps download \
  --bbox=-89.7576695691,35.8571162283,-89.1325306447,36.2386742360 \
  -f geoparquet --type building --release 2026-07-22.0 \
  --output data/reallocation/overture/2026-07-22.0/buildings.parquet
.venv/bin/overturemaps download \
  --bbox=-89.7576695691,35.8571162283,-89.1325306447,36.2386742360 \
  -f geoparquet --type segment --release 2026-07-22.0 \
  --output data/reallocation/overture/2026-07-22.0/segments.parquet
.venv/bin/python scripts/run_reallocation_check.py
```

The download `.state` sidecars are required to verify release and spatial coverage. The current extract has 31,896 building features and 13,783 transportation segments; non-road segments are excluded during rasterization.

`--upstream` selects the original checkout. `--work` selects the intermediate cache; `--output` selects published results. `--prepare-only` constructs just the structural workspace and overlap weights. Changed structural contracts require a separate work directory. Complete nightly intermediates are reused only when their source, code, and static-workspace hashes match.

## Outputs

Published files are under `data/published/reallocation_check/`, or **`~/mnt/ntl/published/reallocation_check/`** on the Mac:

- `native_equivalence_by_date.csv`: reaggregation comparisons on common valid cells.
- `support_by_date.csv`: original versus common observation counts.
- `validation_summary.csv`: M0/M1 pre-event scores for each input variant.
- `model_differences.csv`: paired event prediction and diagnostic differences relative to matched direct A2.
- `cells.csv`: native IDs, centers, and projected footprint areas.
- `comparison_arrays.npz`: all coarse time series, shared masks, model outputs, and grid metadata.
- `native/{direct,upsample,uniform,allocated}/YYYY-MM-DD.tif`: native-resolution radiance inputs.
- `differences/YYYY-MM-DD.tif`: allocated-minus-direct A2 radiance on common valid support. Positive means allocation is brighter; negative means dimmer. These companion layers have their own provenance manifest and can be regenerated with `.venv/bin/python scripts/export_reallocation_differences.py` after a completed comparison run.
- `{variant}/{M0,M1}/YYYY-MM-DD/{center,q05,q95,residual,lower_tail}.tif`: native model diagnostics.
- `manifest.json`: run status, original-source provenance, parameters, and output checksums.

The working cache retains the dimensionless structural gain, fine-to-source mapping, overlap weights, and nightly coarse aggregates. No fine-resolution outage classifications are published.

## Completed check: results

The full February 5–April 17 run is complete. **The aggregated structural proxy is not numerically equivalent to direct A2.** Across 176,412 common valid cell-nights:

| Native-grid input versus direct A2 | MAE | RMSE | Maximum absolute difference |
| --- | ---: | ---: | ---: |
| Upsample-and-aggregate control | 0.00581 | 0.02689 | 0.89036 |
| Uniform smoothing | 0.40608 | 1.86445 | 57.57939 |
| Structural allocation | 0.42100 | 1.96563 | 70.78198 |

All errors are in nW cm⁻² sr⁻¹. Allocation differences substantially exceed grid round-trip error. Uniform smoothing has a similar aggregate error magnitude, consistent with neighborhood averaging being a major contributor; this is not a formal decomposition of every cell's error. Small county-wide signed bias would not establish cellwise conservation.

Common observation support retains 136,391 of 198,364 direct baseline cell-nights and 40,021 of 56,595 event cell-nights. Requiring ten training observations leaves **6,938 supported cells** in every matched model, versus 7,688 in the original direct-only fit. April 2–6 remain without usable event evidence. April 15 also has no complete common-support footprint.

Both M0 and M1 were run on all four inputs. Relative to the matched direct reference, the allocated normal-prediction MAE is 0.3760 for M0 and 0.3278 for M1. Event residual differences have MAE 0.2307 and 0.2188 respectively. These are differences between methods, not errors against an independent outage benchmark.

For the seven-day validation design, allocated M0/M1 have MAE 0.4971/0.4864 versus matched-direct 0.6446/0.6407. Uniform smoothing is nearly the same as allocation at 0.4977/0.4875. Because each model predicts a different input field, these smaller errors do **not** demonstrate better outage detection. Allocated nominal 90% interval coverage is only 80.8%/82.8%; confidence calibration remains unresolved.

Open [04_inspect_reallocation_check.ipynb](../notebooks/04_inspect_reallocation_check.ipynb) for period-specific equivalence metrics, support loss, matched validation, individual-cell histories, and native-grid maps. The next scientific choice is whether to retain this cross-cell allocation as a distinct method or investigate a separately named within-footprint conserving alternative. This run preserves the original method and does not make that substitution.

### Retained 10 m scenes for exploration

`scripts/export_reallocation_fine.py` reconstructs the original nightly allocation
from the cached raw A2 inputs and structural gain without rerunning M0/M1. It
verifies every night's native aggregate against the published `allocated` array,
including missing values, and writes `fine/reallocation_10m.npy` and a separate
provenance manifest under `data/published/reallocation_check/`.

The float32 cube has shape `(72, 4097, 5533)` with axes `(time, row, column)` and
occupies about 6.5 GB. It retains the full original 10 m projected workspace,
including the processing margin outside Dyer. `fine/native_aoi_mask.npy` identifies
fine pixels with positive area overlap with the selected native county cells;
it is available separately and does not mask the saved cube. This is a footprint
mask for the native analysis cells, not a newly rasterized county boundary.

[05_explore_nightly_radiance.ipynb](../notebooks/05_explore_nightly_radiance.ipynb)
loads native A2, native aggregated reallocation, and the fine cube with integer
image indices. The fine cube is a read-only memory map; individual fine dataframes
can be opened on demand. The projected fine grid is not an exact 50-by-50
subdivision of each geographic A2 cell. Retaining these numerical outputs does
not establish independent 10 m radiance observations or perform fine-scale outage
analysis.

Verification passed: 24 tests across the ingestion, counterfactual, and reallocation components, including agreement with the original operator using the production 101 × 101 FFT kernel. The executed notebook contains no cell errors. All 934 listed output hashes and 928 native-grid COGs were checked; raster values match the saved arrays. Native overlap weights sum to one with maximum error below `1.2e-13`. Shared training counts and absence of residual/tail evidence on missing nights were also verified. The audit is saved locally as `data/reallocation/validation_report.json`.
