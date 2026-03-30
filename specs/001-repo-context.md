# Repo Context

## Purpose

This spec gives a new thread the minimum local context needed to work safely in this forked `mvMORPH` repo after the OUM, experimental corrpower BMM, and standalone `mvSIMMAP()` workstreams were integrated on `master`.

## Status

- Status: active
- Last updated: 2026-03-30
- Applies to branch: `master`

## Repository Identity

- Local checkout: `<local checkout path>`
- Remote fork: `origin = https://github.com/jakeberv/mvMORPH-berv.git`
- Upstream source: `upstream = https://github.com/JClavel/mvMORPH.git`

At the time this spec was written:

- the local integration branch is `master`
- the extra corrpower worktree has been removed
- the current repo uses a single local checkout again

## Baseline And Fork Delta

The meaningful comparison baseline for this fork is current `origin/master` and `upstream/master`, which match at:

- `65228fc` `Update NAMESPACE`

The local `master` branch contains merged work on top of that baseline, including:

- OUM covariate support and OUM simulation/validation harnesses
- experimental BMM correlation-aware work ending in the `corrpower` family
- a standalone mixed-process `mvSIMMAP()` fitter for painted trees
- simulation summary reports for the OUM and corrpower workstreams

See `specs/010-fork-delta.md` for the detailed change inventory.

## Important Directories

- `R/`
  - package source code
- `man/`
  - Rd help files
- `tests/`
  - experimental harnesses and launch scripts
- `tests/testthat/`
  - regression tests used as lightweight package checks
- `tests/reports/`
  - OUM simulation summary outputs and figure assets
- `reports/`
  - corrpower simulation summary outputs
- `specs/`
  - durable branch-local context for future threads

## Important Current Surfaces

### OUM In `mvgls`

Main user-facing addition:

- `mvgls(..., model="OUM")` now retains user-supplied covariates alongside regime optima

Important paths:

- `R/mvgls.r`
- `R/classes_methods.r`
- `tests/testthat/test-mvgls-oum-covariates.R`
- OUM simulation harnesses under `tests/experimental_oum_*.R`

### Experimental BMM Correlation-Aware Work

Current public experimental path:

- `mvgls(..., model="BMM", bmm.structure="corrpower")`
- correlation-only nested form:
  - `mvgls(..., model="BMM", bmm.structure="corrpower", bmm.scale=FALSE)`

Retired public selectors:

- `bmm.structure="corrshrink"`
- `bmm.structure="corrstrength"`

Important paths:

- `R/mvgls_bmm_corrpower.r`
- `R/mvgls_bmmcorr_shared.r`
- `R/classes_methods.r`
- `man/corrpower_diagnostics.Rd`
- `tests/experimental_corrpower_*.R`
- `tests/experimental_bmmcorr_*.R`

### Experimental `mvSIMMAP()`

Current public mixed-process SIMMAP path:

- `mvSIMMAP(tree, data, process, process.groups = ..., ...)`
- `simulate(mvSIMMAP_object, nsim = ..., param = list(...))` for mixed-process SIMMAP draws from a fitted or fixed scaffold

Current simulation-first scaffold path:

- `mvSIMMAP(tree, data = NULL, process = ..., optimization = "fixed", param = list(ntraits = ..., ...))`

Important paths:

- `R/mvSIMMAP.r`
- `man/mvSIMMAP.Rd`
- `tests/testthat/test-mvSIMMAP.R`
- `specs/050-mvgls-simmap-mixed-design.md`

## Working Conventions

- Use `master` as the integrated local branch unless there is a reason to branch off again.
- Prefer `tests/testthat` for lightweight regression checks.
- Prefer the dedicated experimental harnesses in `tests/` for simulation or stress-test work.
- Treat the files in `reports/` and `tests/reports/` as generated summaries that capture prior campaign outcomes.
- Treat the BMM correlation-aware code as experimental even though it is merged locally.
- Treat `mvSIMMAP()` as the current mixed-process SIMMAP entry point; `mvgls` integration is still only a draft design.
- Treat `simulate()` on `mvmorph.mixed` objects as the current mixed-process SIMMAP simulation entry point; `mvSIM()` itself is still not mixed-aware.

## First Files To Read For New Work

If the new thread is about OUM:

1. `tests/testthat/test-mvgls-oum-covariates.R`
2. `tests/experimental_oum_detectability_grid.R`
3. `tests/reports/oum_simulation_summary_report.md`

If the new thread is about the new experimental BMM model:

1. `R/mvgls_bmm_corrpower.r`
2. `man/corrpower_diagnostics.Rd`
3. `reports/corrpower_simulation_report.md`
4. `specs/010-fork-delta.md`

If the new thread is about mixed-process SIMMAP models:

1. `man/mvSIMMAP.Rd`
2. `tests/testthat/test-mvSIMMAP.R`
3. `specs/050-mvgls-simmap-mixed-design.md`
4. `specs/010-fork-delta.md`
