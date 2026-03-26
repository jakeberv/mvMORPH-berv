# BMM Simulation Runbook

## Purpose

This runbook is the operational handoff for the experimental BMM correlation-aware work in this fork. It is meant to let a new thread resume simulation work quickly without re-deriving which scripts matter, how they are launched, or what counts as a meaningful result.

## Status

- Status: active
- Last updated: 2026-03-26
- Applies to branch: `master`

## Scope

This runbook covers the current experimental BMM family built around `corrpower`:

- `bmm.structure="corrpower"`
- nested correlation-only mode:
  - `bmm.structure="corrpower", bmm.scale=FALSE`

It covers:

- local smoke checks
- targeted comparison runs
- fossil-recovery campaigns
- the combined historical corrpower campaign scripts

It does not cover:

- the retired public `corrstrength` or `corrshrink` paths except as historical benchmarks
- penalized / PL `corrpower` workflows, which are not implemented

## Current Model Surface

### Public Experimental API

Primary model:

- `mvgls(..., model="BMM", bmm.structure="corrpower")`

Nested no-scale comparator:

- `mvgls(..., model="BMM", bmm.structure="corrpower", bmm.scale=FALSE)`

Retired selectors:

- `bmm.structure="corrshrink"`
- `bmm.structure="corrstrength"`

### Key Files

- `R/mvgls.r`
- `R/mvgls_bmm_corrpower.r`
- `R/mvgls_bmmcorr_shared.r`
- `R/classes_methods.r`
- `man/corrpower_diagnostics.Rd`

### Key Reports

- `reports/corrpower_simulation_report.md`
- `reports/corrpower_simulation_report.html`

That report should be treated as the canonical summary of the major corrpower campaigns run so far.

## Core Harnesses

### Lightweight Functional / Regression Harnesses

Use these first after any code change:

- `tests/experimental_corrpower_mvgls.R`
- `tests/experimental_corrpower_hardening.R`
- `tests/experimental_corrpower_identifiability.R`
- `tests/experimental_corrpower_coronly_mvgls.R`
- `tests/experimental_corrpower_coronly_hardening.R`
- `tests/experimental_corrpower_coronly_identifiability.R`
- `tests/experimental_bmmcorr_family_comparison.R`

### Comparison Harnesses

- `tests/experimental_bmmcorr_model_comparison.R`
  - targeted `corrpower` vs historical `corrstrength` comparison
- `tests/experimental_bmmcorr_family_comparison.R`
  - compares:
    - proportional BMM
    - `corrpower` with `bmm.scale=FALSE`
    - full `corrpower`

### Fossil Harnesses

- `tests/experimental_corrpower_fossil_grid.R`
  - earlier fossil campaign harness
- `tests/experimental_corrpower_fossil_recovery_grid.R`
  - current preferred fossil-truth-recovery harness

## Recommended Local Validation Ladder

### Minimal Fast Check

Run these before anything bigger:

```bash
Rscript tests/experimental_corrpower_mvgls.R
Rscript tests/experimental_corrpower_coronly_mvgls.R
```

### Good Pre-Campaign Check

```bash
Rscript tests/experimental_corrpower_mvgls.R
Rscript tests/experimental_corrpower_hardening.R
Rscript tests/experimental_corrpower_coronly_mvgls.R
Rscript tests/experimental_corrpower_coronly_hardening.R
Rscript tests/experimental_bmmcorr_family_comparison.R
```

### Documentation Check

```bash
Rscript -e 'library(devtools); load_all("<local checkout path>", quiet=TRUE)'
Rscript -e 'tools::parse_Rd("<local checkout path>/man/mvgls.Rd"); tools::parse_Rd("<local checkout path>/man/corrpower_diagnostics.Rd"); tools::parse_Rd("<local checkout path>/man/EIC.Rd")'
```

## Targeted Comparison Campaign

### Purpose

Use the targeted comparison campaign when the question is:

- is the current stabilized `corrpower` fitter better than the old corrstrength benchmark on hard cases?

This is the campaign that historically established stabilized `corrpower` as the preferred experimental path.

### Main Script

- `tests/experimental_bmmcorr_model_comparison.R`

### Launcher

- `tests/run_bmmcorr_targeted_parallel.sh`

### Main Env Vars

- `BMMCORR_TARGETED_CORES`
- `BMMCORR_TARGETED_CHUNKS`
- `BMMCORR_TARGETED_REPS`
- `BMMCORR_TARGETED_OUTPUT_ROOT`
- `BMMCORR_TARGETED_R_LIB`

The launcher internally passes through:

- `CORRPOWER_COMPARISON_FULL=TRUE`
- chunking env vars for the comparison script

### Typical Local / Remote Launch

```bash
bash tests/run_bmmcorr_targeted_parallel.sh
```

### Outputs

Under the output root, expect:

- `comparison/chunk_*.csv`
- `comparison/chunk_*.rds`
- `comparison_combined.csv`
- `campaign_summary.txt`
- `logs/joblog.tsv`

### Success Criteria

At a minimum:

- all chunk jobs finish successfully
- `comparison_combined.csv` is produced
- `campaign_summary.txt` reports `comparison_combined=TRUE`

Scientifically, this campaign is mainly used to compare:

- covariance recovery
- regime-summary recovery
- pathology rate
- boundary rate

## Fossil-Recovery Campaign

### Purpose

Use the fossil-recovery campaign when the question is:

- do fossil terminals improve truth recovery under the stabilized corrpower family?

This is the current preferred fossil simulation run for the BMM workstream.

### Main Script

- `tests/experimental_corrpower_fossil_recovery_grid.R`

### Launcher

- `tests/run_corrpower_fossil_recovery_parallel.sh`

### Main Env Vars

- `CORRPOWER_FOSSIL_RECOVERY_CORES`
- `CORRPOWER_FOSSIL_RECOVERY_CHUNKS`
- `CORRPOWER_FOSSIL_RECOVERY_REPS`
- `CORRPOWER_FOSSIL_RECOVERY_OUTPUT_ROOT`
- `CORRPOWER_FOSSIL_RECOVERY_R_LIB`
- `CORRPOWER_FOSSIL_RECOVERY_INSTALL_NCPUS`
- `CORRPOWER_FOSSIL_RECOVERY_LAMBDA_SCALE`
- `CORRPOWER_FOSSIL_RECOVERY_LAMBDA_CORR_POWER`

### Typical Launch

```bash
bash tests/run_corrpower_fossil_recovery_parallel.sh
```

### Outputs

Under the output root, expect:

- `chunks/chunk_*.csv`
- `chunks/chunk_*.rds`
- `fossil_recovery_combined.csv`
- `campaign_summary.txt`
- `logs/joblog.tsv`
- `logs/chunk_*.log`

### What The Combined CSV Should Contain

This campaign is valuable because it records truth-recovery metrics for full-tree vs extant-only fits. It should be used to compare things like:

- regime covariance recovery
- regime-summary recovery
- boundary behavior
- pathological-scale behavior

### Success Criteria

Operationally:

- all chunk jobs finish with exit code `0`
- `fossil_recovery_combined.csv` exists
- `campaign_summary.txt` reports `fossil_recovery_combined=TRUE`

Scientifically:

- fossil scenarios should be judged against extant-only matched fits
- the most important outputs are the recovery-to-truth metrics, not raw logLik comparisons across different datasets

## Historical Combined Corrpower Campaign

### Purpose

This launcher exists to run a historical two-part campaign:

- hard-case comparison chunks
- earlier fossil-grid chunks

### Launcher

- `tests/run_corrpower_campaign_parallel.sh`

### When To Use It

Usually only when reproducing the earlier branch history or regenerating the campaign summarized in the existing report. For new fossil work, prefer:

- `tests/run_corrpower_fossil_recovery_parallel.sh`

### Outputs

- `comparison_combined.csv`
- `fossil_combined.csv`
- `campaign_summary.txt`
- per-chunk files under `comparison/` and `fossil/`

## Practical Remote Notes

### Package Installation Pattern

The reliable pattern for high-parallel runs in this repo is:

1. clean `src/*.o` and `src/mvMORPH.so`
2. install required R packages into the target user library
3. run `R CMD INSTALL --preclean`
4. execute workers against the installed package

Avoid large-worker campaigns that rely on many concurrent `pkgload::load_all()` compiles against the same source tree.

### User-Level R Library

When running as `tguser` remotely, prefer the user-level library already used in prior campaigns:

- `$HOME/R/x86_64-pc-linux-gnu-library/4.5`

## What To Read Before Designing New BMM Simulations

1. `reports/corrpower_simulation_report.md`
2. `tests/experimental_bmmcorr_family_comparison.R`
3. `tests/experimental_corrpower_fossil_recovery_grid.R`

## Common Pitfalls

- Do not treat `corrstrength` as an active public model; it is historical only.
- Do not interpret raw log-likelihood differences across fossil and extant-only datasets as the main scientific result.
- Do not assume PL / penalized `corrpower` exists yet.
- Do not launch large worker pools against `pkgload::load_all()` on a shared source tree; install once first.

## Recommended Next-Simulation Questions

If resuming simulation work on the experimental BMM line, the most natural next questions are:

1. whether `corrpower` should gain penalized / PL support
2. whether the `bmm.scale=FALSE` nested comparator is sufficient as the simpler alternative
3. whether fossil gains remain strong under more realistic design choices or targeted harder corners
