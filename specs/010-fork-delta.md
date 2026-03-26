# Fork Delta: OUM And Experimental BMM Work

## Purpose

This spec is the main onboarding document for future threads. It describes, in one place, the major work added on top of the upstream `mvMORPH` fork baseline, how that work evolved, what the current public surface is, and which files matter most.

## Status

- Status: active
- Last updated: 2026-03-26
- Applies to branch: `master`

## Executive Summary

Since the fork diverged from upstream/master, two major workstreams were built and then merged back together:

1. `mvgls` OUM support was extended so user-supplied covariates are retained alongside regime optima, and that work was backed by new validation and simulation harnesses.
2. A substantial experimental BMM correlation-aware modeling line was developed:
   - started as `corrshrink`
   - became `corrstrength`
   - was replaced by `corrpower`
   - then stabilized
   - then extended with a simpler correlation-only nested mode using `bmm.scale=FALSE`

The current integrated recommendation is:

- for OUM work, use the merged covariate-capable `mvgls(..., model="OUM")` path
- for experimental BMM correlation-aware work, use `bmm.structure="corrpower"`
- for the nested correlation-only variant, use `bmm.structure="corrpower", bmm.scale=FALSE`

## Baseline

Relative to current `origin/master` / `upstream/master`:

- baseline commit: `65228fc` `Update NAMESPACE`

The local `master` now includes all merged work described below.

## Workstream 1: OUM Covariates And Validation

### What Changed

The fork added support for OUM fits in `mvgls` that preserve user-supplied covariates alongside regime-specific optima, instead of collapsing the fit to regime terms only.

This change matters because downstream helpers should continue to work for regression-style OUM fits, not just intercept-only regime models.

### Main User-Facing Effect

Users can now fit:

- `mvgls(Y ~ x, ..., model="OUM")`

and still expect the fit object to preserve both:

- regime optimum information
- covariate/regression information

for downstream summaries and model-comparison helpers.

### Important Files

- `R/mvgls.r`
- `R/classes_methods.r`
- `tests/testthat/test-mvgls-oum-covariates.R`

### Validation And Simulation Additions

Beyond the regression tests, a larger OUM simulation layer was added:

- `tests/experimental_oum_detectability_grid.R`
- `tests/experimental_oum_theta_recovery_grid.R`
- `tests/experimental_oum_null_pathology_grid.R`
- `tests/experimental_bmm_oum_model_selection_grid.R`
- `tests/reports/oum_simulation_summary_report.md`
- `tests/reports/oum_simulation_summary_report.html`

These scripts support questions like:

- when OU vs OUM is distinguishable
- when OUM theta recovery is reliable
- when null-pathology or false-selection problems appear
- how BMM and OUM compare in model-selection grids

### Important Commit Landmarks

- `55ba8af` `Support OUM covariates in mvgls`
- `0e58e6e` `Add OUM regression validation and reporting`
- `3cc3539` `Add OUM model-selection simulation harnesses`
- `e4d84ae` `Add OUM theta recovery simulation harness`
- `e37d795` `Add OUM simulation summary report`
- `89ce1d3` `Add OUM high-dimensional guardrails`
- `cf8b072` `Fix OUM merge cleanup on master`

## Workstream 2: Experimental BMM Correlation-Aware Models

### Phase A: `corrshrink`

The first experimental correlation-aware BMM work was introduced under the `corrshrink` name. This added:

- a regime reference control
- simulation and `EIC` support
- identifiability, anchor-sensitivity, and hardening harnesses
- derived regime summaries and diagnostics

Important early commits:

- `31b2b9b` `Add experimental corr-shrink mvgls simulate and EIC support`
- `b893392` `Add corr-shrink identifiability study script`
- `308033e` `Add corr-shrink BMM reference regime control`
- `4d54aac` `Add corr-shrink anchor sensitivity workflow`
- `4a4e13e` `Add corr-shrink regime summary reporting`
- `2ce5910` `Add corr-shrink diagnostics and corr-strength utilities`

### Phase B: `corrstrength`

The original `corrshrink` surface was then replaced by `corrstrength`, using a normalized shared-template parameterization. The intent was to separate:

- regime-specific overall rate scale
- regime-specific correlation strength

This version worked technically but turned out to be practically unstable in several scenarios, especially weak-signal and higher-dimensional settings.

Key transition commit:

- `11dd36f` `Replace experimental corrshrink with corrstrength BMM`

### Phase C: `corrpower`

The next redesign replaced linear correlation-strength scaling with a matrix-power formulation:

- shared baseline correlation matrix
- regime-specific `corr_power`
- optionally regime-specific `scale`

This became the experimental `corrpower` family.

Important files:

- `R/mvgls_bmm_corrpower.r`
- `R/mvgls_bmmcorr_shared.r`
- `R/classes_methods.r`
- `man/corrpower_diagnostics.Rd`

Key initial commit:

- `2d50cec` `Add experimental corrpower BMM option`

### Phase D: Stabilization

Simulation work showed that the first `corrpower` version was better than `corrstrength` at avoiding catastrophic failures, but boundary behavior was still too high. A stabilization pass then added weak regularization and smoother scale handling, which materially improved performance.

Key milestones:

- `2087124` `Add corrpower campaign harnesses`
- `aea3ff2` `Stabilize corrpower and add fossil recovery grid`
- `86ab22e` `Add corrpower simulation report`

### Phase E: Public API Cleanup

After the simulation work supported `corrpower` over `corrstrength`, the experimental public surface was cleaned up:

- `corrstrength` was retired from the public API
- internal helper names were moved toward neutral `bmmcorr` naming
- the nested correlation-only mode was added
- the original `corrpower_coronly` flag was then simplified to:
  - `bmm.structure="corrpower", bmm.scale=FALSE`

Key milestones:

- `4db9e5e` `Retire experimental corrstrength public API`
- `dee58dc` `Rename internal BMM correlation helpers`
- `4e2ebc0` `Add experimental corrpower correlation-only BMM mode`
- `dc3ddaa` `Simplify corrpower correlation-only API`

## Current Experimental BMM Public Surface

### Supported

Primary experimental entry point:

- `mvgls(..., model="BMM", bmm.structure="corrpower")`

Nested correlation-only form:

- `mvgls(..., model="BMM", bmm.structure="corrpower", bmm.scale=FALSE)`

Supported downstream helpers in the current experimental path:

- `simulate()`
- `AIC()`
- `BIC()`
- `GIC()`
- `EIC()`
- `predict()`
- `ancestral()`
- `residuals(type="normalized")`
- `confint()`
- `corrpower_diagnostics()`

### Retired Or Explicitly Unsupported

Retired public selectors:

- `bmm.structure="corrshrink"`
- `bmm.structure="corrstrength"`

Still unsupported in the current experimental BMM correlation-aware path:

- penalized estimation / PL workflows
- `REML`
- measurement error
- `FCI`
- `vcov(type="coef")`

### Main Implementation Files

- `R/mvgls.r`
  - dispatch and argument parsing, including `bmm.scale`
- `R/mvgls_bmm_corrpower.r`
  - core `corrpower` fitter and nested no-scale mode
- `R/mvgls_bmmcorr_shared.r`
  - archived corrstrength/shared helper path
- `R/classes_methods.r`
  - diagnostics, `confint`, summaries, simulation/prediction plumbing
- `R/fun.r`
  - simulation wrapper support for the experimental BMM family

## Simulation And Reporting Layer Added In This Fork

### Corrpower / BMM Reports

- `reports/corrpower_simulation_report.md`
- `reports/corrpower_simulation_report.html`

These summarize the main arc of experimental BMM development:

- corrstrength instability
- first fossil comparisons
- unstabilized corrpower vs corrstrength
- stabilized corrpower vs corrstrength
- stabilized fossil recovery results

### Corrpower Harnesses

- `tests/experimental_corrpower_mvgls.R`
- `tests/experimental_corrpower_hardening.R`
- `tests/experimental_corrpower_identifiability.R`
- `tests/experimental_corrpower_fossil_grid.R`
- `tests/experimental_corrpower_fossil_recovery_grid.R`
- `tests/experimental_corrpower_coronly_mvgls.R`
- `tests/experimental_corrpower_coronly_hardening.R`
- `tests/experimental_corrpower_coronly_identifiability.R`
- `tests/experimental_bmmcorr_model_comparison.R`
- `tests/experimental_bmmcorr_family_comparison.R`

### Remote Campaign Launchers Added

- `tests/run_corrpower_campaign_parallel.sh`
- `tests/run_bmmcorr_targeted_parallel.sh`
- `tests/run_corrpower_fossil_recovery_parallel.sh`

These were added to support large remote runs via GNU `parallel` and one-time installed-package workflows.

## Current Interpretation For New Threads

If a new thread is about the experimental BMM work, it should assume:

- `corrpower` is the preferred experimental model family
- `bmm.scale=FALSE` is the nested correlation-only comparator
- `corrstrength` and `corrshrink` are historical only
- simulation evidence exists locally in `reports/corrpower_simulation_report.md`

If a new thread is about OUM:

- the merged `master` already includes the OUM covariate support work
- the OUM simulation and summary-report infrastructure is already present under `tests/experimental_oum_*.R` and `tests/reports/`

## New Thread Fast Start

For the experimental BMM line, read these in order:

1. `specs/001-repo-context.md`
2. `reports/corrpower_simulation_report.md`
3. `R/mvgls_bmm_corrpower.r`
4. `man/corrpower_diagnostics.Rd`
5. `tests/experimental_corrpower_mvgls.R`
6. `tests/experimental_bmmcorr_family_comparison.R`

For the OUM line, read these in order:

1. `tests/testthat/test-mvgls-oum-covariates.R`
2. `tests/reports/oum_simulation_summary_report.md`
3. `tests/experimental_oum_detectability_grid.R`
4. `tests/experimental_oum_theta_recovery_grid.R`

## Technical Debt To Keep In Mind

- `R/mvgls_bmmcorr_shared.r` still contains legacy corrstrength-oriented internal names and aliases for archived benchmark material.
- `R/classes_methods.r` still contains some compatibility wrappers and internal corrstrength naming.
- The experimental BMM line is merged locally but remains scientifically experimental.
- The current simulation layer is rich, but not all scenarios are reduced to lightweight unit tests.
