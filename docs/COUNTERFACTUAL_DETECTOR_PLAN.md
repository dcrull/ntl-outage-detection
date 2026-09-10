# Native-grid counterfactual radiance detector: literature review and draft design

Draft, 2026-09-09. This is a focused literature review and proposed experiment, not a fitted detector or evidence of outage-detection accuracy. It develops the counterfactual component of the [analysis outline](../ANALYSIS_PLAN_DRAFT.md). TING remains reserved for independent evaluation.

## 1. Analysis unit and prediction target

Use the **actual VNP46A2 Collection 002 grid**, preserving its transform, cell boundaries, and values without resampling. It is a 15-arc-second geographic grid, conventionally called 500 m; it is not a grid of exact 500 × 500 m squares. Store stable cell identifiers using tile plus source row/column, rather than crop indices alone. Native grid spacing also does not establish independent resolving power at that spacing.

The current county crop has 7,868 cell centers in Dyer County. Retain the existing center-based inclusion rule for the primary analysis, with county-intersection fractions for boundary sensitivity and area calculations. Calculate area geodesically or in an appropriate equal-area projection, not square degrees. See [ingestion details](VIIRS_INGESTION.md).

The eventual sensor benchmark would assign sensors to these same cells and retain reporting denominators. **The TING extract currently contains outage polygons, not sensor locations or complete power-state telemetry.** For this extract, build a separate polygon-to-cell intersection table, retaining event IDs, source timestamps, and overlap fractions. Union overlapping geometries at the chosen time before calculating area fractions. Polygon overlap measures recorded footprint coverage, not the fraction of customers without power. Absence of a polygon is not a confirmed negative. Timestamp and restoration semantics need resolution before constructing intervals; see [TING caveats](TING_INGESTION.md).

For cell i and night t, estimate the distribution of corrected radiance that would be expected **if its preceding lighting regime continued without a new disruption**, conditional on usable observation conditions. Denote its central prediction by m(i,t), and the observed radiance by y(i,t). The residual is:

`r(i,t) = y(i,t) - m(i,t)`

Negative residuals are evidence of dimming relative to that reference. Store predictive quantiles as well as a central estimate. A lower-tail probability describes unusual radiance under the reference model; it is not automatically a probability of a power outage.

The desired causal target is radiance with power maintained under otherwise comparable conditions. Our historical prediction only approximates it: evacuation, destruction, changes in operating hours, and residual observation artifacts can also create deficits. No candidate below uniquely identifies power loss from radiance alone.

## 2. What the literature contributes

The papers below separate direct nighttime-light applications from statistical methods we would adapt. The proposed choices in later sections are our design judgments, not reported findings for Dyer County.

| Study | Relevant contribution | Implication for this experiment |
| --- | --- | --- |
| [Román et al. (2019), *Satellite-based assessment of electricity restoration efforts in Puerto Rico after Hurricane Maria*, PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218883) | Uses pre-hurricane radiance as a recovery reference; native-grid metrics and multi-date aggregation address large, prolonged disruption. | A simple pre-event reference is an essential comparator. Its recovery metrics do not establish a calibrated nightly counterfactual or validate short outages obscured by cloud. |
| [Li et al. (2022), *Continuous monitoring of nighttime light changes based on daily NASA’s Black Marble product suite*, Remote Sensing of Environment](https://doi.org/10.1016/j.rse.2022.113269) | Introduces VZA-COLD: viewing-zenith-angle stratification combined with continuous change monitoring of Black Marble time series. | Viewing geometry should be considered when predicting normal radiance, even when A2 provides the radiance input. |
| [Li et al. (2026), *Satellite imagery reveals increasing volatility in human night-time activity*, Nature](https://www.nature.com/articles/s41586-026-10260-w) | Uses robust harmonic models within viewing-angle strata, A2 radiance, and A1 geometry. Its long-record change analysis requires 14 consecutive anomalous observations and excludes changes returning within one month. | Borrow the observation model, not these persistence rules: our target includes short outages and our record cannot support this full long-term configuration. |
| [Xu, Cao and Chen (2026), *Assessing Power System Reliability Using Anomaly Detection in Daily Nighttime Light Data*, Remote Sensing](https://doi.org/10.3390/rs18091417) | Groups observations by viewing angle and derives adaptive per-pixel references from high-brightness samples. The authors acknowledge validation favoring large, long outages and limitations under persistent cloud. | Directly relevant outage literature. Test an angle-aware reference, but do not import its thresholds or presume its validation transfers to short Dyer events. |
| [Brodersen et al. (2015), *Inferring causal impact using Bayesian structural time-series models*, Annals of Applied Statistics](https://kaybrodersen.github.io/publications/Brodersen_2015_AOAS.pdf) | Combines state-space components and regression controls to generate posterior counterfactual trajectories. Controls must be unaffected by the intervention. | Basis for a compact state-space forecast. A version using only cell history remains a model-based reference, with weaker causal identification than a credible controlled design. |
| [Ben-Michael, Feller and Rothstein (2021), *The Augmented Synthetic Control Method*, JASA](https://arxiv.org/abs/1811.04170) | Builds a weighted comparison trajectory and uses an outcome model, such as ridge regression, to correct imperfect pre-treatment balance. | Candidate for removing shared variation, conditional on credible donor cells and adequate overlapping observations. |
| [Athey et al. (2021), *Matrix Completion Methods for Causal Panel Data Models*, JASA](https://arxiv.org/abs/1710.10251) | Estimates untreated outcomes using a regularized low-rank panel structure, linking panel regression and synthetic-control approaches. | A later extension if the panel and donor assumptions support it. Missing observations and unknown outage exposure are different problems; completing a matrix does not resolve both. |

This is a targeted shortlist, not a systematic or exhaustive review. Some publisher pages were access-limited; the review uses accessible primary article text, author-hosted papers, and publisher-indexed abstracts/method excerpts. No external study’s reported accuracy is treated as expected performance here.

## 3. Proposed methods to compare

### M0. Robust constant reference — required comparator

Predict each cell’s pre-event median on every subsequent night. Estimate its normal variation robustly, using the median absolute deviation or empirical residual distribution, with pooled calibration where individual records are too sparse.

This uses the existing cache and provides an interpretable test of whether additional modeling earns its complexity. It assumes stable normal lighting and observation effects over the window. A robust spread alone is not a calibrated predictive interval, and a zero spread needs a defensible noise floor. Do not treat near-zero baseline cells as reliably detectable through percentage drops.

### M1. Parsimonious state-space forecast — first main candidate

Model normal brightness as a slowly changing latent level with observation noise:

`y(i,t) = level(i,t) + beta(i) × observation_covariates(i,t) + error(i,t)`

`level(i,t) = level(i,t-1) + process_noise(i,t)`

Start with the local-level model and no extra covariates. Strongly regularize how quickly the level can change; consider sharing variance information among cells with similar pre-event brightness and variability. Robust errors are a proposed extension to reduce sensitivity to occasional poor retrievals. Add a damped trend only if held-out forecasting improves.

Keep the daily calendar intact: a missing observation skips the measurement update, while the latent state still advances through elapsed days. For this retrospective pilot, fit using the pre-event period and propagate the normal trajectory forward without updating it from target-cell event radiance. Otherwise an outage can be absorbed into the estimated normal level. Prediction intervals must include observation noise and parameter/state uncertainty.

This is the strongest initial candidate for explicit uncertainty with the existing cache, **not a presumption that it will outperform M0**. With only a few dozen usable observations, variance estimates and forecast drift may be unstable. No annual harmonic model is justified by eight weeks of history.

### M2. Viewing-angle-conditioned robust regression — main observation-model alternative

Predict each cell’s radiance from a simple viewing-angle relationship, with a tightly constrained time term only if supported. Use a low-complexity curve or partially pooled coefficients; independently fitting several angle bins per cell would divide an already short record into very small samples.

This adapts the angle-aware nighttime-light literature above. It requires adding matched A1 viewing geometry while retaining **A2 corrected radiance as the outcome**. Match dates, platform, source grid, and valid acquisition support; verify that ancillary geometry corresponds to the A2 retrieval. Avoid extrapolating confidently to angles absent from training.

Test whether angle explains held-out variation before retaining it. M2 can initially compete with M1 as a simpler regression; an eventual state-space model can include the useful angle terms. An upper-brightness reference is a sensitivity analysis, not the default: preferentially selecting bright samples can raise the expected reference and make ordinary nights look anomalously dim.

### M3. Synthetic-control prediction — conditional extension

For each target cell, learn a small weighted combination of other cells that reproduces its pre-event trajectory. Predict the target during the event from those donors’ contemporaneous, quality-screened radiance. Compare constrained weights with a regularized augmented version only if needed.

This can capture shared nightly variation that target history misses. Its usefulness depends on donors remaining unaffected and their relationship to the target remaining stable. Nearby cells may share power infrastructure and storm impacts. Good pre-event fit does not prove they are valid controls.

Select donors using pre-event satellite behavior and a prespecified spatial rule, without TING or post-event target dimming. Under the satellite-only constraint, donor non-exposure remains an assumption to stress-test. The cached full h09v05 tiles allow a limited search beyond Dyer without downloading the whole state; the output AOI can remain Dyer.

Require sufficient common valid nights and define missing-donor handling in advance. Do not silently reweight to whichever donors remain clear. A predetermined fallback to M0/M1 must be identified in the outputs. Test sensitivity to donor sets and simulated donor contamination. If credible donors cannot be supported, omit causal claims and defer this method.

## 4. Model comparison before inspecting TING agreement

The [completed cache audit](../data/viirs/VNP46A2.002/validation_summary.json) covers February 5–April 1 (56 baseline dates) and April 2–17 (16 subsequent dates). The median county cell has **26 usable baseline nights and 8 subsequent nights**. Initial notebook support of ten baseline nights is descriptive, not an established modeling threshold.

1. Build one cell-night table containing raw corrected radiance, QA, usability, source provenance, and pre-event support. Preserve valid zeros. Keep TING in a separate evaluation table.
2. Compare M0 and M1 using forward, calendar-block holdouts entirely before April 2. Start with roughly four weeks for training; evaluate subsequent blocks, including multi-day forecast horizons. Fit preprocessing, pooling, and model choices within each training fold. Avoid random date splits and future-data interpolation.
3. Score prediction error and predictive-distribution calibration, including lower-tail behavior, interval coverage, width, and forecast horizon. Stratify by pre-event brightness and support. Evaluate all models on common eligible observations and also report each model’s coverage; dropping hard cells must not look like improved accuracy.
4. Hold out whole nights/contiguous gaps to reproduce cloud-related loss of information. Preserve spatial dependence when assessing uncertainty; thousands of correlated cells are not thousands of independent experiments. With this short record, extreme-tail calibration will remain weak. Additional pre-event history is preferable if results are unstable.
5. Add M2 after geometry ingestion and retain it only if it improves held-out performance or calibration. Investigate M3 separately with explicit donor assumptions. Prefer the simpler model when gains are unclear.
6. Freeze the forecast specification before event evaluation. Tune subsequent anomaly/state thresholds with pre-event placebo windows and explicit simulated dimming experiments, never TING. Placebo alerts are not proven false positives, because historical outages may be unknown; simulated deficits test algorithm behavior rather than real-world attribution.

The fixed April 2 cutoff defines this pilot’s forecast origin, not known cell-level treatment. Confirm acquisition times before matching the local-evening event. A deployable detector needs a separate rolling-training policy that avoids absorbing newly emerging disruptions; report live filtering separately from retrospective smoothing.

## 5. Interface to change-point/state inference

Each model should export cell ID, acquisition date/time support, observed radiance and QA, predicted central radiance, predictive quantiles, residual, normal-model lower-tail score, forecast horizon, training support, method, and any fallback/insufficient-support flag.

Only usable observed radiance contributes new dimming evidence. Cloud-masked observations remain unknown, and the A2 gap-filled layer is not substituted as an outage measurement. The later state model can propagate uncertainty across missing nights, but cannot turn an unobserved interval into confirmed outage or recovery. Keep relative deficits undefined where normal brightness is too low for a stable ratio.

Extent maps must distinguish detected dimming, no detected dimming, and insufficient evidence. Onset/recovery are bounded by observed nights; a proposed state model must explain how its assumptions affect duration estimates. Spatial persistence may support a candidate disruption, but spatially coherent cloud errors can also produce apparent agreement.

**Critical event limitation:** April 2–6 have no usable Dyer observations under current screening, while the intersecting TING records are concentrated on April 2–3. An outage wholly within that blind interval may leave no observable satellite deficit. Counterfactual forecasting can be assessed on baseline holdouts, but this event may not support meaningful satellite sensitivity or duration validation against those short TING records. Classify this as an observation limitation, not a normal-power result.

## 6. Decisions for the next implementation step

- Implement M0 and a strongly regularized M1 first; select on held-out prediction and calibration, not a preferred causal label.
- Keep the output on the native grid and retain uncertain/dark cells explicitly.
- Decide whether to add A1 geometry for M2 and whether sparse folds warrant extending the pre-event period.
- Treat M3 as conditional on a documented donor study; defer matrix completion and large machine-learning models.
- Define minimum support, predictive interval calibration, meaningful dimming magnitude, and the later state model before evaluating TING agreement.
- Resolve the TING timestamp semantics and benchmark coverage limits independently of detector fitting. The current polygon extract does not support full cell-level false-positive/false-negative accounting.
