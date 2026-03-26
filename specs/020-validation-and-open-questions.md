# Validation And Open Questions

## Purpose

This spec summarizes what has already been validated in the merged fork and what the most likely next questions are. It is intentionally practical: it tells a future thread what is already known well enough not to rediscover from scratch.

## Status

- Status: active
- Last updated: 2026-03-26
- Applies to branch: `master`

## Current Validation Picture

### OUM

Validated locally through:

- `tests/testthat/test-mvgls-oum-covariates.R`

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

### Documentation

The merged repo cleanup pass confirmed Rd parsing for:

- `man/mvgls.Rd`
- `man/corrpower_diagnostics.Rd`
- `man/EIC.Rd`

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
- the most valuable next work should probably build on `corrpower`, not revisit `corrstrength`

## Main Open Questions

### OUM

- How much of the OUM simulation surface should be converted into lighter-weight regression tests?
- Are there additional high-dimensional edge cases that should become explicit guardrails rather than being handled only in simulation scripts?

### Experimental BMM Corrpower

- Should the next methodological step be penalized/PL support for `corrpower`?
- Is the nested no-scale mode enough as the simpler comparator, or is another model-family comparison still needed?
- How much additional cleanup is worth doing in the archived `corrstrength` internal layer versus leaving it as historical scaffolding?

## Recommended Next Checks Before Big New Work

If starting a new feature branch, a good lightweight confidence pass is:

1. `Rscript tests/experimental_corrpower_mvgls.R`
2. `Rscript tests/experimental_corrpower_coronly_mvgls.R`
3. `Rscript tests/testthat/test-mvgls-oum-covariates.R`

If working on the corrpower fitter or diagnostics, also run:

4. `Rscript tests/experimental_corrpower_hardening.R`
5. `Rscript tests/experimental_bmmcorr_family_comparison.R`

If working on OUM simulation logic, also run:

6. a tiny env-restricted smoke run of `tests/experimental_oum_theta_recovery_grid.R`

## Notes For New Threads

- Do not assume there is still a separate corrpower worktree; that cleanup is already done.
- Do not assume the repo uses `main`; the integrated local branch is `master`.
- Read `specs/010-fork-delta.md` before changing public experimental BMM behavior.
