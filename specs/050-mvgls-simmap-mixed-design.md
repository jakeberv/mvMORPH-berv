# MVGLS Mixed SIMMAP Design

## Purpose

This spec captures the current design for bringing mixed BM/OU/OUM/EB regime support on painted SIMMAP trees into `mvgls` while preserving formula-based regression. It is meant to save a future thread from having to rediscover the architectural constraints, the likely first implementation shape, and the main technical risks.

## Status

- Status: draft
- Last updated: 2026-03-30
- Applies to branch: `master`

This is a forward-looking design note. The currently implemented public mixed-process SIMMAP surface is the standalone `mvSIMMAP()` fitter, not `mvgls(..., model="SIMMAPmixed")`.

## Scope

- In scope:
  - a tentative `mvgls(..., model="SIMMAPmixed", process=..., process.groups=...)` design
  - a first-pass parameterization that fits the existing `mvgls` GLS architecture
  - integration points across `mvgls`, helper code, and penalized methods
  - expected validation work and open design questions
- Out of scope:
  - implementation
  - locking the final public argument names
  - a fully nonseparable `mvgls` model with regime-specific matrix `sigma` and matrix `alpha`
  - a PCMBase-style pruning rewrite

## Background

- An experimental standalone `mvSIMMAP()` fitter now exists in `R/mvSIMMAP.r`.
- `mvSIMMAP()` accepts a painted SIMMAP tree plus a per-regime process assignment among `BM`, `OU`, `OUM`, and `EB`.
- `mvSIMMAP()` also accepts `process.groups`, so painted regimes can either keep separate parameter blocks or share common process blocks.
- `mvSIMMAP(..., data = NULL, optimization = "fixed", param = list(ntraits = ...))` can now act as a simulation scaffold for mixed-process SIMMAP models.
- For `OU`, shared `process.groups` collapse painted regimes into one OU regime for the likelihood, with shared `alpha`, `sigma`, and `theta`.
- For `OUM`, shared `process.groups` share OU dynamics while retaining painted-regime-specific optima.
- The practical standalone default for OU/OUM groups is now `decomp = "scalarPositive"`, so each OU/OUM process group gets one positive pull parameter unless the user opts into a richer `alpha` structure.
- The standalone mixed fitter currently uses a fixed-root-style mean structure with one global root vector estimated by GLS; broader `mvOU`-style `root` / `vcv` options are not yet part of the mixed-process API.
- Internally, `mvSIMMAP()` builds a dense vectorized design matrix `D` and a dense covariance matrix `V`, then evaluates the likelihood through `loglik_mvmorph()`.
- Optimized `mvSIMMAP()` fits now expose `LogLik` as a standard `logLik` object, so generic tools like `BIC()` work without changing the wider `mvMORPH` comparison machinery.
- Fitted objects now also carry richer post-fit diagnostics, including Hessian status labels and optional jittered-restart metadata.
- `simulate()` on `mvSIMMAP` objects now reuses that dense mixed-model structure to generate multivariate Gaussian draws under the fitted or manually overridden mixed-process parameterization.
- `mvgls()` currently supports formula-based regression with a shared trait covariance and a smaller menu of phylogenetic row-covariance models.
- Existing `mvgls` `OUM` support is the closest precedent: it augments the design matrix with regime-weight columns and estimates a scalar OU strength parameter.
- PCMFit and PCMBase are not just analogous here; they are the main conceptual source for the mixed-Gaussian branchwise formulation used in this fork's `mvSIMMAP()` work. What was borrowed is the modeling idea of composing branch-local Gaussian transitions across mapped segments. What was not borrowed is their implementation code or their pruning engine.
- The `mvSIMMAP()` work in this fork should therefore be described as a dense `mvMORPH` implementation inspired by PCMBase/PCMFit, not as an independent rediscovery of the same formulation.

## Current Behavior

### Public API

- There is no `mvgls` API yet for mixed BM/OU/OUM/EB SIMMAP models.
- The current mixed-process fitting entry point is `mvSIMMAP(tree, data, process, ...)`.
- The current mixed-process simulation entry point is `simulate(mvSIMMAP_object, ...)`.
- There is still no direct `mvSIM(..., model="SIMMAPmixed")` selector.

### Internal implementation

- `mvgls()` routes supported models through `.prepModel()`, `.setBounds()`, `.startGuess()`, `.corrStr()`, and `.transformTree()`.
- `.transformTree()` currently assumes that a model can be represented either as a transformed tree or as a row-covariance whitening step that returns whitened `X` and `Y`.
- `OUM` is special because its mean structure is assembled by `.mvgls_oum_design()`.
- The current `mvSIMMAP()` code instead propagates segment-specific transitions and builds:
  - a dense covariance `V`
  - a dense mean design `D`
- That exact `mvSIMMAP()` parameterization is not directly compatible with `mvgls`, because `mvgls` assumes a separable covariance of the form `C \otimes R` with an `n x m` mean design.

## Proposed User API

- Tentative entry point:
  - `mvgls(formula, data, tree, model="SIMMAPmixed", process=..., process.groups=..., ...)`
- Required inputs:
  - `tree`: a `simmap` object with `mapped.edge` and `maps`
  - `process`: a named vector keyed by SIMMAP regime, with values in `c("BM", "OU", "OUM", "EB")`
- Likely optional inputs:
  - `process.groups`: a named vector keyed by SIMMAP regime when shared process blocks are desired
  - SIMMAP-specific starting values and bounds through the existing `start`, `low`, `up`, or `param` conventions
  - `root` and `root_std`, following the current OUM semantics as closely as possible

## First Implementation Target

The first `mvgls` implementation should preserve the existing separable GLS structure rather than trying to embed the full standalone `mvSIMMAP()` parameterization.

The practical target is:

- one shared trait covariance matrix `R`
- scalar regime-specific process parameters:
  - a rate multiplier for each regime
  - `alpha_r > 0` for OU regimes
  - `beta_r` for EB regimes
- a tip-level row covariance `C` built from SIMMAP segment recursion
- an OUM-like design matrix containing:
  - one root column when appropriate
  - one `theta.regime` column for each OU regime
  - the ordinary non-intercept regression predictors from the formula

This is the best fit to the current `mvgls` architecture because it lets the model remain separable:

- row covariance: `C(process, tree, map)`
- trait covariance: `R`
- full covariance: `C \otimes R`

It also preserves the core `mvgls` idea that regression effects and OU optima live in the mean structure, while the phylogenetic correlation lives in the row covariance.

## Proposed Statistical Formulation

For a branch segment in regime `r` of duration `t`, define scalar transition pieces:

- BM:
  - `Phi_r(t) = 1`
  - `Q_r(t) = rate_r * t`
- OU:
  - `Phi_r(t) = exp(-alpha_r * t)`
  - `Theta_r(t) = 1 - Phi_r(t)`
  - `Q_r(t) = rate_r * (1 - exp(-2 * alpha_r * t)) / (2 * alpha_r)`
- EB:
  - `Phi_r(t) = 1`
  - `Q_r(t) = rate_r * eb_factor(beta_r, t, start_age)`

Because these transition multipliers are scalar, the multivariate process stays compatible with a shared trait covariance `R`.

The resulting implementation can mirror the recursion already used in `mvSIMMAP()`, but with scalar node quantities rather than full `p x p` matrices:

- scalar ancestry multiplier `A`
- scalar OU-optimum weights `B`
- scalar accumulated variance term `V`

From those quantities, the code can build:

- a tip-level row covariance matrix `C`
- a SIMMAP-specific mean-design block for root and OU optima

That block can then be combined with the formula design matrix inside `mvgls`, matching the role that `.mvgls_oum_design()` already plays for `OUM`.

## Integration Points

### `R/mvgls.r`

- accept the new model name
- require a `simmap` tree
- carry `process` and any SIMMAP-specific precalculations into `corrModel`
- adjust the effective predictor count `m` when the design is expanded by OU regime columns

### `R/fun.r`

- extend `.prepModel()` to precompute:
  - mapped edge segments
  - node ordering
  - MRCA lookup
  - regime indices
- add a helper analogous to `.mvgls_oum_design()` for the mixed SIMMAP mean structure

### `R/penalized.r`

- extend `.corrStr()` to allow the new model
- extend `.transformTree()` with a `SIMMAPmixed` branch that:
  - builds the row covariance `C`
  - computes a whitening factor from `C`
  - builds and whitens the augmented design matrix
  - whitens `Y`
  - returns the determinant and weighted matrices
- extend `.setBounds()` and `.startGuess()` to handle block parameters for rates, OU `alpha`, and EB `beta`

### `R/mvSIMMAP.r`

- reuse or refactor internal SIMMAP helpers where practical, especially:
  - segment parsing
  - EB age accumulation across contiguous segments
  - branch-local transition logic

## Main Challenges

- Architectural fit:
  - the current standalone `mvSIMMAP()` parameterization is nonseparable, while `mvgls` is built around separable GLS
- Mean versus covariance split:
  - OU optima need to enter the `mvgls` design matrix, not a fully vectorized `D`
- Parameter bookkeeping:
  - the new model needs stable ordering, starts, bounds, printing, and summary behavior for a heterogeneous parameter vector
- EB continuity:
  - EB scaling must accumulate across contiguous segments of the same EB regime, as already handled in `mvSIMMAP()`
- Identifiability:
  - root handling, short regimes, duplicated or weakly informed OU optima, and flat EB surfaces will all need guardrails
- Performance:
  - even the restricted version likely needs a dense row covariance `C`, which is slower than the current pruning-based tree transforms
- Penalized methods:
  - `LL` should be the easiest first target
  - `LOOCV` and `EmpBayes` should work in principle once whitening is stable, but runtime and numerical behavior will need explicit validation
- Long-term extensibility:
  - if the package later wants regime-specific full `sigma_r` and matrix `alpha_r` inside `mvgls`, that likely requires a more invasive nonseparable likelihood path rather than a small extension

## Validation

Once implementation begins, the expected validation layers are:

- intercept-only tests that reduce to existing `BM`, `EB`, `OUM`, or `BMM` behavior where overlap exists
- regression tests with one or two predictors under mixed regime assignments
- equivalence checks against scalar or separable special cases of `mvSIMMAP()`
- numerical tests for:
  - non-ultrametric trees
  - short mapped regimes
  - measurement error
  - weakly identified OU and EB settings

Expected first-version limitations:

- no regime-specific full `sigma` matrices inside `mvgls`
- no regime-specific matrix `alpha`
- likely slower than the current `BM`, `OU`, and `OUM` `mvgls` paths

## Open Questions

- What should the public model name be: `SIMMAPmixed`, `SIMMAP`, or something closer to current naming?
- Should the first version include regime-specific rate scalars for every regime, or only OU and EB-specific parameters plus one shared overall rate?
- Should the first pass support only `LL`, with penalized methods added after the model stabilizes?
- How much of the current `mvSIMMAP()` helper layer should be moved into shared internals before implementation starts?
- Should the fitted object expose the SIMMAP-expanded design explicitly for diagnostics and prediction helpers?

## Next Steps

1. Freeze the restricted statistical model for version 1, especially the rate-scalar decision.
2. Refactor reusable SIMMAP helpers out of `R/mvSIMMAP.r` if that reduces duplicate logic.
3. Prototype the scalar row-covariance and mixed-design builders before touching the `mvgls()` front end.
4. Implement `LL` first, then decide whether `LOOCV` and `EmpBayes` are stable enough for the same branch.

## References

- `R/mvgls.r`
- `R/penalized.r`
- `R/fun.r`
- `R/mvSIMMAP.r`
- commit `efda128` (`Add experimental SIMMAP mixed-process fitter`)
- commit `bad77a5` (`Refine mvSIMMAP grouping semantics and summaries`)
- commit `9f7dfe2` (`Improve mvSIMMAP diagnostics and fit-object compatibility`)
