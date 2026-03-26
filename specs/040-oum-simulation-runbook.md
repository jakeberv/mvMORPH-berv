# OUM Simulation Runbook

## Purpose

This runbook is the operational handoff for the OUM simulation and validation work in this fork. It is meant to let a new thread resume OUM model-selection, theta-recovery, and null-pathology studies without having to rediscover which scripts to run or how to shrink them down for smoke tests.

## Status

- Status: active
- Last updated: 2026-03-26
- Applies to branch: `master`

## Scope

This runbook covers the OUM extensions added in this fork:

- `mvgls(..., model="OUM")` with retained user covariates
- OUM-vs-OU detectability studies
- OUM theta recovery studies
- OUM null-pathology studies
- broader BMM-vs-OUM model-selection calibration grids

It does not cover:

- the experimental corrpower BMM campaigns
- upstream `mvOU()` workflows outside the `mvgls`-based fork extensions

## Core OUM Surfaces

### Main User-Facing Extension

The key fork-level user-facing change is:

- `mvgls(..., model="OUM")` now retains user-supplied covariates alongside regime-specific optima

Key validation file:

- `tests/testthat/test-mvgls-oum-covariates.R`

### Main OUM Simulation Scripts

- `tests/experimental_oum_detectability_grid.R`
  - OU-vs-OUM detectability and misspecification-focused grid
- `tests/experimental_oum_theta_recovery_grid.R`
  - focused theta-recovery and null runs
- `tests/experimental_oum_null_pathology_grid.R`
  - false-selection / pathology-oriented null grid
- `tests/experimental_bmm_oum_model_selection_grid.R`
  - broader calibration grid spanning BM, OU, BMM, and OUM

### Legacy Broad Covariate Grid

- `tests/experimental_oum_covariate_grid.R`
  - early broad covariate-oriented OUM sweep kept mainly for baseline smoke coverage
  - a later unmerged rewrite expanded this script with richer predictor families (`cont_indep`, `cont_confounded`, `disc_binary`, `disc3_sparse_confounded`), explicit near-null effect-size sweeps, and generic `delta50` / `delta80` sharpness summaries
  - those ideas were mostly superseded by the focused scripts above, so future work should extend the focused suites rather than revive the broad-grid rewrite wholesale
  - if a future thread wants a single broad covariate sweep again, the only idea worth selectively reusing is the generic near-null sharpness summary layer

### Summary Report

- `tests/reports/oum_simulation_summary_report.md`
- `tests/reports/oum_simulation_summary_report.html`

This report is the best first summary of what the OUM simulation work already found.

## Recommended Local Validation Ladder

### Minimal Check

```bash
Rscript tests/testthat/test-mvgls-oum-covariates.R
```

### Good Pre-Simulation Check

```bash
Rscript tests/testthat/test-mvgls-oum-covariates.R
OUM_REC_CORES=1 OUM_REC_PHASES=focused_null OUM_REC_METHODS=LL OUM_REC_REGIMES=2 OUM_REC_BALANCE=imbalanced OUM_REC_NULL_PREDICTOR=cont_indep OUM_REC_COVARIATE_EFFECT=0.00 OUM_REC_NULL_NTIPS=60 OUM_REC_NULL_P=8 OUM_REC_NULL_ALPHA=0.50 OUM_REC_NULL_REPS=1 Rscript tests/experimental_oum_theta_recovery_grid.R
```

That tiny `theta_recovery` smoke command was already used successfully during the merged-repo cleanup.

## Script-By-Script Run Notes

### 1. `experimental_oum_detectability_grid.R`

#### Main Question

- when should OUM be distinguishable from OU or other nearby alternatives?

#### Main Env Vars

- `OUM_DET_CORES`
- `OUM_DET_SEED`
- `OUM_DET_PHASES`
- `OUM_DET_METHODS`
- `OUM_DET_BALANCE`
- `OUM_DET_NULL_PREDICTOR`
- `OUM_DET_FRONTIER_PREDICTOR`
- `OUM_DET_MISSPEC_PREDICTOR`
- `OUM_DET_COVARIATE_EFFECT`
- `OUM_DET_REGIMES`
- `OUM_DET_NULL_NTIPS`
- `OUM_DET_NULL_P`
- `OUM_DET_NULL_ALPHA`
- `OUM_DET_FRONTIER_NTIPS`
- `OUM_DET_FRONTIER_P`
- `OUM_DET_FRONTIER_ALPHA`
- `OUM_DET_FRONTIER_THETA`
- `OUM_DET_MISSPEC_NTIPS`
- `OUM_DET_MISSPEC_P`
- `OUM_DET_MISSPEC_ALPHA`
- `OUM_DET_MISSPEC_THETA`
- `OUM_DET_MISSPEC_RATE`
- `OUM_DET_NULL_REPS`
- `OUM_DET_FRONTIER_REPS`
- `OUM_DET_MISSPEC_REPS`

#### Typical Smoke Launch

```bash
OUM_DET_CORES=1 \
OUM_DET_PHASES=oum_null \
OUM_DET_METHODS=LL \
OUM_DET_REGIMES=2 \
OUM_DET_BALANCE=balanced \
OUM_DET_NULL_PREDICTOR=cont_indep \
OUM_DET_COVARIATE_EFFECT=0.00 \
OUM_DET_NULL_NTIPS=60 \
OUM_DET_NULL_P=8 \
OUM_DET_NULL_ALPHA=0.50 \
OUM_DET_NULL_REPS=1 \
Rscript tests/experimental_oum_detectability_grid.R
```

#### Outputs

The script prints summaries to stdout and can also save CSV outputs when the relevant output env vars are set near the bottom of the file.

### 2. `experimental_oum_theta_recovery_grid.R`

#### Main Question

- when is regime-theta recovery reliable for OUM fits?

#### Main Env Vars

- `OUM_REC_CORES`
- `OUM_REC_SEED`
- `OUM_REC_PHASES`
- `OUM_REC_METHODS`
- `OUM_REC_REGIMES`
- `OUM_REC_BALANCE`
- `OUM_REC_NULL_PREDICTOR`
- `OUM_REC_RECOVERY_PREDICTOR`
- `OUM_REC_COVARIATE_EFFECT`
- `OUM_REC_NULL_NTIPS`
- `OUM_REC_NULL_P`
- `OUM_REC_NULL_ALPHA`
- `OUM_REC_NULL_REPS`
- `OUM_REC_RECOVERY_NTIPS`
- `OUM_REC_RECOVERY_P`
- `OUM_REC_RECOVERY_ALPHA`
- `OUM_REC_RECOVERY_THETA`
- `OUM_REC_RECOVERY_REPS`

#### Proven Smoke Command

```bash
OUM_REC_CORES=1 \
OUM_REC_PHASES=focused_null \
OUM_REC_METHODS=LL \
OUM_REC_REGIMES=2 \
OUM_REC_BALANCE=imbalanced \
OUM_REC_NULL_PREDICTOR=cont_indep \
OUM_REC_COVARIATE_EFFECT=0.00 \
OUM_REC_NULL_NTIPS=60 \
OUM_REC_NULL_P=8 \
OUM_REC_NULL_ALPHA=0.50 \
OUM_REC_NULL_REPS=1 \
Rscript tests/experimental_oum_theta_recovery_grid.R
```

#### Success Signal

Operationally:

- script exits `0`
- prints a one-row summary
- issue rows stay at `0` in the smoke case

### 3. `experimental_oum_null_pathology_grid.R`

#### Main Question

- where do null scenarios produce false selection or pathological fits?

#### Main Env Vars

- `OUM_NULL_CORES`
- `OUM_NULL_SEED`
- `OUM_NULL_METHODS`
- `OUM_NULL_NTIPS`
- `OUM_NULL_P`
- `OUM_NULL_REGIMES`
- `OUM_NULL_BALANCE`
- `OUM_NULL_PREDICTOR`
- `OUM_NULL_COVARIATE_EFFECT`
- `OUM_NULL_ALPHA`
- `OUM_NULL_REPS`

#### Typical Smoke Launch

```bash
OUM_NULL_CORES=1 \
OUM_NULL_METHODS=LL \
OUM_NULL_NTIPS=60 \
OUM_NULL_P=24 \
OUM_NULL_REGIMES=2 \
OUM_NULL_BALANCE=balanced \
OUM_NULL_PREDICTOR=cont_indep \
OUM_NULL_COVARIATE_EFFECT=0.00 \
OUM_NULL_ALPHA=0.50 \
OUM_NULL_REPS=1 \
Rscript tests/experimental_oum_null_pathology_grid.R
```

### 4. `experimental_bmm_oum_model_selection_grid.R`

#### Main Question

- across BM, OU, BMM, and OUM, which model wins where?

#### Main Env Vars

- `CAL_GRID_CORES`
- `CAL_GRID_SEED`
- `CAL_GRID_METHODS`
- `CAL_GRID_PHASES`
- `CAL_GRID_BALANCE`
- `CAL_GRID_COVARIATE_EFFECT`
- `CAL_GRID_REGIMES`
- `CAL_GRID_NULL_REPS`
- `CAL_GRID_BMM_REPS`
- `CAL_GRID_OUM_REPS`
- `CAL_GRID_MIXED_REPS`
- `CAL_GRID_NULL_NTIPS`
- `CAL_GRID_NULL_P`
- `CAL_GRID_NULL_OU_ALPHA`
- `CAL_GRID_BMM_NTIPS`
- `CAL_GRID_BMM_P`
- `CAL_GRID_BMM_RATE`
- `CAL_GRID_OUM_NTIPS`
- `CAL_GRID_OUM_P`
- `CAL_GRID_OUM_ALPHA`
- `CAL_GRID_OUM_THETA`
- `CAL_GRID_MIXED_NTIPS`
- `CAL_GRID_MIXED_P`
- `CAL_GRID_MIXED_ALPHA`
- `CAL_GRID_MIXED_THETA`
- `CAL_GRID_MIXED_RATE`

#### Typical Smoke Launch

```bash
CAL_GRID_CORES=1 \
CAL_GRID_PHASES=oum_frontier \
CAL_GRID_METHODS=LL \
CAL_GRID_BALANCE=balanced \
CAL_GRID_COVARIATE_EFFECT=0.00 \
CAL_GRID_REGIMES=2 \
CAL_GRID_OUM_REPS=1 \
CAL_GRID_OUM_NTIPS=60 \
CAL_GRID_OUM_P=24 \
CAL_GRID_OUM_ALPHA=0.25 \
CAL_GRID_OUM_THETA=0.0100 \
Rscript tests/experimental_bmm_oum_model_selection_grid.R
```

## Output Conventions

The OUM simulation scripts are mostly self-contained `Rscript` entry points that:

- print compact summaries to stdout
- optionally write CSV outputs when output env vars are provided

Unlike the corrpower campaigns, these OUM scripts do not currently share one unified launcher pattern in this repo. For a new thread, it is usually simplest to:

1. run a tiny env-restricted smoke test locally
2. then decide whether to build or reuse a dedicated launcher for a bigger campaign

## Practical Remote Notes

Each OUM script already contains its own runtime bootstrap pattern:

- install missing R packages into a user-level library if needed
- clean staged `src` artifacts
- `pkgload::load_all(..., compile = TRUE, recompile = TRUE)`

That makes them fairly self-sufficient, but it also means:

- for large remote campaigns, repeated compile work can be expensive
- if a future thread wants higher-throughput OUM campaigns, it may be worth building one-time installed-package launchers similar to the corrpower shell scripts

## Success Criteria

For smoke tests:

- exit code `0`
- compact summary printed
- no issue rows in the trivial cases

For larger campaigns:

- scenario summary tables are written or printed as expected
- issue rows are interpretable and not dominated by infrastructure problems
- model-selection and theta-recovery summaries match the design question of interest

## What To Read Before Designing New OUM Simulations

1. `tests/reports/oum_simulation_summary_report.md`
2. `tests/testthat/test-mvgls-oum-covariates.R`
3. `tests/experimental_oum_detectability_grid.R`
4. `tests/experimental_oum_theta_recovery_grid.R`
5. `tests/experimental_bmm_oum_model_selection_grid.R`

## Common Pitfalls

- Do not assume the OUM simulation scripts are lightweight; most need env restriction for quick smoke runs.
- Do not jump straight to a large campaign before checking `tests/testthat/test-mvgls-oum-covariates.R`.
- Do not assume the BMM experimental results automatically transfer to OUM questions; keep the workstreams conceptually separate.

## Recommended Next-Simulation Questions

If resuming the OUM line, the most natural next questions are:

1. which OUM simulation surfaces should become lighter-weight automated regression tests
2. whether additional high-dimensional guardrails should be enforced in code rather than only documented in simulation summaries
3. whether the existing OUM scripts should get one-time-installed remote launchers similar to the corrpower campaigns
