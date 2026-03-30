# Validation And Open Questions

## Purpose

This spec summarizes what has already been validated in the merged fork and what the most likely next questions are. It is intentionally practical: it tells a future thread what is already known well enough not to rediscover from scratch.

## Status

- Status: active
- Last updated: 2026-03-30
- Applies to branch: `master`

## Current Validation Picture

### OUM

Validated locally through:

- the filtered `testthat` check for `tests/testthat/test-mvgls-oum-covariates.R`

Supported by larger simulation harnesses and reports:

- `tests/experimental_oum_detectability_grid.R`
- `tests/experimental_oum_theta_recovery_grid.R`
- `tests/experimental_oum_null_pathology_grid.R`
- `tests/experimental_bmm_oum_model_selection_grid.R`
- `tests/reports/oum_simulation_summary_report.md`

In the merged repo cleanup pass, a tiny smoke run of:

- `tests/experimental_oum_theta_recovery_grid.R`

was also run successfully on a single focused-null scenario.

### Experimental BMM Corrpower

The corrpower line has been validated locally through dedicated harnesses, including:

- `tests/experimental_corrpower_mvgls.R`
- `tests/experimental_corrpower_hardening.R`
- `tests/experimental_corrpower_identifiability.R`
- `tests/experimental_corrpower_coronly_mvgls.R`
- `tests/experimental_corrpower_coronly_hardening.R`
- `tests/experimental_corrpower_coronly_identifiability.R`
- `tests/experimental_bmmcorr_model_comparison.R`
- `tests/experimental_bmmcorr_family_comparison.R`

The merged repo cleanup pass also reran:

- `tests/experimental_corrpower_mvgls.R`
- `tests/experimental_corrpower_coronly_mvgls.R`

and both passed.

### Experimental `mvSIMMAP`

Validated locally through:

- the filtered `testthat` check for `tests/testthat/test-mvSIMMAP.R`

That regression file currently covers:

- BM/BMM, OU/OUM, and EB special-case likelihood equivalence where the standalone fitter overlaps existing model families
- `process.groups` semantics for shared OU, OUM, and EB parameter blocks
- compact summary/print behavior, Hessian diagnostics, retry handling, and standard `logLik` / `BIC()` compatibility
- scaffold-based mixed-process simulation through `simulate(mvSIMMAP_object, ...)`, including manual parameter overrides and `data=NULL` fixed-scaffold construction

There is not yet a larger simulation campaign or summary report for `mvSIMMAP()` comparable to the OUM and corrpower workstreams.

What changed recently:

- mixed-process SIMMAP simulation is now available through `simulate()` on `mvSIMMAP` objects
- `mvSIMMAP(..., data = NULL, optimization = "fixed", ...)` can now be used as a lightweight simulation scaffold constructor
- `mvSIM()` itself is still limited to the legacy single-family model selectors and is not yet a direct mixed-process SIMMAP simulator

### Documentation

This docs refresh confirmed Rd parsing for:

- `man/mvgls.Rd`
- `man/corrpower_diagnostics.Rd`
- `man/EIC.Rd`
- `man/mvSIMMAP.Rd`

## What The BMM Simulation Work Already Established

The local report:

- `reports/corrpower_simulation_report.md`

should be treated as the canonical summary of the experimental BMM campaigns run so far.

The main conclusions already supported by simulation are:

1. `corrstrength` was scientifically interesting but too unstable to keep as the preferred experimental path.
2. `corrpower` improved robustness relative to `corrstrength`.
3. A stabilization pass materially improved `corrpower`.
4. The nested no-scale version remains a useful comparator:
   - `bmm.structure="corrpower", bmm.scale=FALSE`
5. Under stabilized `corrpower`, fossil terminals improved recovery of the true regime covariance structure in a large fraction of simulated fossil scenarios.

## Current Working Assumptions

If a new thread starts today, it should assume:

- `corrpower` is the default experimental BMM correlation-aware family
- `bmm.scale=FALSE` is the simplest nested comparator
- `corrstrength` and `corrshrink` are historical only
- `mvSIMMAP()` is the current mixed-process SIMMAP entry point, while `mvgls` integration remains a draft design
- the most valuable next work should probably build on `corrpower`, not revisit `corrstrength`

Operational runbooks for resuming simulation work now live in:

- `specs/030-bmm-simulation-runbook.md`
- `specs/040-oum-simulation-runbook.md`

## Main Open Questions

### OUM

- How much of the OUM simulation surface should be converted into lighter-weight regression tests?
- Are there additional high-dimensional edge cases that should become explicit guardrails rather than being handled only in simulation scripts?

### Experimental BMM Corrpower

- Should the next methodological step be penalized/PL support for `corrpower`?
- Is the nested no-scale mode enough as the simpler comparator, or is another model-family comparison still needed?
- How much additional cleanup is worth doing in the archived `corrstrength` internal layer versus leaving it as historical scaffolding?

### Experimental `mvSIMMAP`

- How much of the current `mvSIMMAP()` surface needs simulation-based validation beyond the existing regression tests?
- Which moderate tree sizes and trait dimensions remain practical under the dense likelihood implementation?
- How much shared helper logic should be extracted before any `mvgls` integration work begins?

## Recommended Next Checks Before Big New Work

If starting a new feature branch, a good lightweight confidence pass is:

1. `Rscript tests/experimental_corrpower_mvgls.R`
2. `Rscript tests/experimental_corrpower_coronly_mvgls.R`
3. `Rscript -e 'testthat::test_local(".", filter = "mvgls-oum-covariates", reporter = "summary")'`

If working on the corrpower fitter or diagnostics, also run:

4. `Rscript tests/experimental_corrpower_hardening.R`
5. `Rscript tests/experimental_bmmcorr_family_comparison.R`

If working on `mvSIMMAP()` or mixed-process SIMMAP docs/logic, also run:

6. `Rscript -e 'testthat::test_local(".", filter = "mvSIMMAP", reporter = "summary")'`

If working on OUM simulation logic, also run:

7. a tiny env-restricted smoke run of `tests/experimental_oum_theta_recovery_grid.R`

## Notes For New Threads

- Do not assume there is still a separate corrpower worktree; that cleanup is already done.
- Do not assume the repo uses `main`; the integrated local branch is `master`.
- Read `specs/010-fork-delta.md` before changing public experimental BMM behavior.
- Mixed-process SIMMAP support currently lives in standalone `mvSIMMAP()`, not in `mvgls()`.
