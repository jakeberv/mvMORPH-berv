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

There is still not yet a large automated simulation campaign for `mvSIMMAP()` at the scale of the OUM and corrpower workstreams.

There is now, however, a durable local benchmark summary for the standalone mixed fitter in:

- `specs/060-mvSIMMAP-simulation-validation.md`

That spec records the current 1D, 2D, and 3D matched-tree recovery checks, including the main takeaway that `theta` and background `BM` covariance recover more cleanly than derived-regime `OU` pull strengths, plus the local `rpf` versus `inverse` timing comparison for the dense implementation.

What changed recently:

- mixed-process SIMMAP simulation is now available through `simulate()` on `mvSIMMAP` objects
- `mvSIMMAP(..., data = NULL, optimization = "fixed", ...)` can now be used as a lightweight simulation scaffold constructor
- the default OU/OUM alpha decomposition is now `scalarPositive`, with `diagonalPositive` and `cholesky` available as opt-in richer alternatives
- the current mixed fitter uses one global GLS-estimated root vector, corresponding to the fixed-root side of the older `mvOU` root semantics
- reusable exact-tree benchmark drivers now exist for the canonical `BM -> OU` case, the shared-OUM case, and a mixed `BM + EB + OUM` case
- in the first exact-tree `BM + EB + OUM` checks, the `OUM` signal was much easier to distinguish than the `EB` signal; `EB` can still collapse toward its `beta = 0` BM boundary in some replicates even when the generating model includes a stronger `EB` regime
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
- How strong and how deep does an `EB` regime need to be before the mixed fitter reliably prefers `EB` over a regime-specific BM block instead of letting `beta` collapse to `0`?
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
7. `MVSIMMAP_BENCH_REPS=1 Rscript tests/experimental_mvsimmap_recovery_benchmark.R`
8. `MVSIMMAP_OUM_EB_BENCH_REPS=1 Rscript tests/experimental_mvsimmap_oum_eb_benchmark.R`

If working on OUM simulation logic, also run:

9. a tiny env-restricted smoke run of `tests/experimental_oum_theta_recovery_grid.R`

## Notes For New Threads

- Do not assume there is still a separate corrpower worktree; that cleanup is already done.
- Do not assume the repo uses `main`; the integrated local branch is `master`.
- Read `specs/010-fork-delta.md` before changing public experimental BMM behavior.
- Mixed-process SIMMAP support currently lives in standalone `mvSIMMAP()`, not in `mvgls()`.
