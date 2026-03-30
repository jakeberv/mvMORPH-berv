# mvSIMMAP Simulation Validation

## Purpose

This spec records what has already been learned from local mixed-process `mvSIMMAP()` simulation-and-recovery checks, so future threads do not have to reconstruct the benchmark setup or rerun the same first-pass analyses just to recover the basic conclusions.

## Status

- Status: active
- Last updated: 2026-03-30
- Applies to branch: `master`

## Scope

- In scope:
  - local exact-tree simulation checks for the standalone mixed-process `mvSIMMAP()` fitter
  - the scaffold-based simulation workflow now available through `simulate(mvSIMMAP_object, ...)`
  - practical recovery behavior for a background `BM` regime plus a derived `OU` regime
  - timing and method comparisons for the dense methods currently exposed by `mvSIMMAP()`
- Out of scope:
  - tree-uncertainty studies
  - mismatched-tree inference
  - full-scale simulation reports with hundreds or thousands of replicates
  - `mvgls(..., model="SIMMAPmixed")`, which is still design-only

## Background

The current mixed-process public API lives in standalone `mvSIMMAP()`. It can now also act as a simulation scaffold through:

- `mvSIMMAP(tree, data = NULL, optimization = "fixed", ...)`
- `simulate(mvSIMMAP_object, param = list(...))`

That simulation path reuses the exact tree, painted SIMMAP mapping, `process`, and `process.groups` structure retained inside the `mvSIMMAP` object.

The local checks summarized here were designed as an initial identifiability and workflow sanity pass, not as a final scientific claim. The emphasis was:

- exact matched-tree recovery first
- moderate tree size
- small numbers of traits
- strong derived-regime OU pull
- clear logging of convergence, Hessian status, and parameter bias

A reusable driver for the canonical exact-tree benchmark now lives in:

- `tests/experimental_mvsimmap_recovery_benchmark.R`

## Canonical Local Benchmark Setup

All matched-tree recovery runs summarized below used the same basic tree and regime layout:

- tree seed: `20260330`
- tree generator: `phytools::pbtree(n = 50, scale = 1)`
- painted shift node: `53`
- total tips: `50`
- derived `OU` clade size: `28` tips
- tree height: `1`
- derived clade root depth: `0.04175815`
- remaining derived clade height: `0.9582418`
- regime assignment:
  - background regime `A = BM`
  - derived regime `B = OU`
- inference tree:
  - exact matched SIMMAP tree used for simulation and fitting

Common fitting settings:

- `optimization = "L-BFGS-B"`
- `param = list(decomp = "diagonalPositive", decompSigma = "cholesky")`
- `control = list(maxit = 500, retry.unreliable = TRUE, retry.max = 2, retry.jitter = 0.05, retry.seed = 1)`

Common simulation pattern:

1. build a fixed scaffold with `mvSIMMAP(..., data = NULL, optimization = "fixed")`
2. generate one or more replicates with `simulate(scaffold, seed = ..., param = true_param)`
3. refit the same mixed model on the exact matched SIMMAP tree
4. summarize convergence, Hessian status, bias, and RMSE

## Current Validation Results

### 1D exact-tree recovery

Setup:

- traits: `1`
- method: `inverse`
- replicates: `10`
- true parameters:
  - `root = 0`
  - `theta.B = 1.5`
  - `sigma_A = 0.04`
  - `sigma_B = 0.04`
  - `alpha_B = 3`

Mean recovery across the 10 replicates:

- `root = 0.0052`
- `theta.B = 1.4855`
- `sigma_A = 0.0370`
- `sigma_B = 0.0567`
- `alpha_B = 5.5131`

RMSE:

- `root = 0.1087`
- `theta.B = 0.1153`
- `sigma_A = 0.0145`
- `sigma_B = 0.0315`
- `alpha_B = 4.8125`

Operational notes:

- all `10/10` fits converged
- all `10/10` Hessians were reported as `reliable`
- the main weakness was the expected `alpha_B` / `sigma_B` tradeoff
- in this setting, the mean and background-BM pieces were much easier to recover than the derived OU pull strength

### 2D exact-tree recovery

Setup:

- traits: `2`
- method: `inverse`
- replicates: `10`
- true `theta.B = (1.5, -1.0)`
- true `sigma_A = sigma_B =`

```text
[ 0.040  0.012
  0.012  0.030 ]
```

- true `alpha_B = diag(3.0, 2.5)`

Mean recovery across the 10 replicates:

- `theta.B = (1.5173, -0.9926)`
- mean `sigma_A =`

```text
[ 0.0370  0.0096
  0.0096  0.0283 ]
```

- mean `sigma_B =`

```text
[ 0.0579  0.0168
  0.0168  0.0367 ]
```

- mean `alpha_B = diag(5.7181, 3.6391)`

Selected RMSE values:

- `thetaB_y1 = 0.2139`
- `thetaB_y2 = 0.0804`
- `sigmaB_11 = 0.0332`
- `sigmaB_22 = 0.0120`
- `alphaB_11 = 5.1257`
- `alphaB_22 = 1.9035`

Operational notes:

- all `10/10` fits converged
- all `10/10` Hessians were `reliable`
- one OU pull component improved relative to the 1D case, but the first trait still showed a strong `alpha` / `sigma` ridge
- the second trait's `alpha` / `sigma` correlation was much weaker than the first trait's

### 3D exact-tree recovery with `inverse`

Setup:

- traits: `3`
- method: `inverse`
- replicates: `10`
- true `theta.B = (1.5, -1.0, 0.75)`
- true `sigma_A = sigma_B =`

```text
[ 0.040  0.012  0.008
  0.012  0.030  0.009
  0.008  0.009  0.025 ]
```

- true `alpha_B = diag(3.0, 2.5, 2.0)`
- elapsed time on 4 cores: `245.975` seconds

Mean recovery across the 10 replicates:

- `theta.B = (1.5612, -0.9919, 0.7458)`
- mean `sigma_A =`

```text
[ 0.0370  0.0096  0.0035
  0.0096  0.0283  0.0061
  0.0035  0.0061  0.0216 ]
```

- mean `sigma_B =`

```text
[ 0.0592  0.0170  0.0164
  0.0170  0.0369  0.0123
  0.0164  0.0123  0.0327 ]
```

- mean `alpha_B = diag(5.8610, 3.7337, 3.2046)`

Selected RMSE values:

- `thetaB_y1 = 0.3632`
- `thetaB_y2 = 0.0771`
- `thetaB_y3 = 0.1016`
- `sigmaB_11 = 0.0342`
- `sigmaB_22 = 0.0115`
- `sigmaB_33 = 0.0165`
- `alphaB_11 = 5.2219`
- `alphaB_22 = 2.0834`
- `alphaB_33 = 2.7371`

Operational notes:

- all `10/10` fits converged
- all `10/10` Hessians were `reliable`
- the third trait improved interpretability and made the compact summary output more useful
- it did not eliminate OU identifiability problems for the strongest-confounded dimensions

### 3D exact-tree recovery with `rpf`

The same 3D datasets were rerun with `method = "rpf"` instead of `method = "inverse"`.

Results:

- elapsed time on 4 cores: `205.31` seconds
- speedup relative to `inverse`: about `16.5%`
- recovered means were numerically almost identical to the `inverse` run
- all `10/10` fits converged
- all `10/10` Hessians were `reliable`

Practical conclusion:

- among the methods currently available to `mvSIMMAP()`, `rpf` is the best default performance choice when the fit is numerically well behaved
- however, the current implementation is still dense-only, so even `rpf` remains much slower than a true sparse or pruning-based implementation would be

### 3D exact-tree model comparison against global `BM1`

The same 3D `BM -> OU` benchmark datasets were also compared against a single-regime Brownian baseline:

- mixed model: `mvSIMMAP(process = c(A = "BM", B = "OU"))`
- baseline model: `mvBM(model = "BM1")`
- fitting method: `rpf`
- replicates: `10`

Results:

- mean `AIC_mixed = -361.0092`
- mean `AIC_BM1 = -257.1560`
- mean `delta_AIC = AIC_BM1 - AIC_mixed = 103.8532`
- median `delta_AIC = 104.8889`
- `delta_AIC` range: `85.4453` to `116.8239`
- mixed model wins: `10/10`
- mixed model wins by more than 10 AIC units: `10/10`
- mean two-model Akaike weight for the mixed model: effectively `1.0`

Interpretation:

- for this exact-tree mixed-process generating scenario, the mixed model is not just slightly better than global `BM1`; it is decisively preferred on every replicate
- the main signal is not a small penalty tradeoff but a large log-likelihood gain that overwhelms the mixed model's extra parameters
- the benchmark driver now reports this comparison directly alongside parameter-recovery summaries

## What We Have Learned So Far

### Estimation behavior

- Exact matched-tree recovery is clearly viable for the current standalone `mvSIMMAP()` implementation.
- In the canonical mixed `BM -> OU` simulation, the fitted mixed model is decisively favored over a global `BM1` alternative by AIC.
- The easiest parameters to recover in these local checks were:
  - `theta`
  - background `BM` covariance terms
- Derived-regime `OU` covariance terms were somewhat biased upward but still broadly informative.
- Derived-regime `OU` pull parameters (`alpha`) were consistently the hardest parameters to recover.
- Increasing trait dimension helped some `alpha` components, but did not remove the classic `alpha` / `sigma` confounding ridge.

### Practical workflow

- `mvSIMMAP(..., data = NULL, optimization = "fixed")` is good enough to serve as a simulation scaffold.
- `simulate(mvSIMMAP_object, param = list(...))` is sufficient for exact-tree validation loops.
- `tests/experimental_mvsimmap_recovery_benchmark.R` can now rerun the canonical 3D exact-tree benchmark directly, with env overrides for method, trait count, replicate count, and output directory.
- That benchmark driver now includes a built-in global `BM1` baseline and prints per-replicate and summary AIC comparisons.
- The new compact multivariate summary output is useful for debugging and for reading 2D/3D fits without digging through raw matrices.
- `parameters(fit)` is useful for flattening multivariate estimates into a simulation-summary table; `fit$sigma`, `fit$alpha`, and `coef(fit)` remain the right interfaces when full matrices are desired.

### Performance

- The dense implementation becomes expensive quickly as trait dimension rises.
- In the local matched-tree benchmark:
  - 3D + 10 replicates + `inverse` took about `246` seconds on 4 cores
  - the same run with `rpf` took about `205` seconds on 4 cores
- So method choice matters, but the larger bottleneck is still the dense likelihood strategy itself.

## Recommended Current Practice

For future local `mvSIMMAP()` simulation checks:

1. Start with exact matched-tree validation before adding tree misspecification.
2. Prefer `method = "rpf"` unless there is a specific numerical reason to use `inverse` or `pseudoinverse`.
3. Record all of:
   - convergence
   - Hessian status
   - per-parameter bias and RMSE
4. Pay particular attention to the derived-regime `alpha` / `sigma` relationship rather than interpreting `alpha` in isolation.
5. Use a reasonably large and deep derived clade when testing OU recovery; tiny derived clades are too easy to mistake for algorithmic failure when the issue is really weak information.

## Main Remaining Open Questions

- How much better does `alpha` recovery get with larger trees and more derived-regime tips?
- How sensitive are these results to shift depth and derived-clade placement?
- How much does recovery degrade when the painted SIMMAP tree is not exact?
- Should a reusable `tests/reports/` markdown summary be generated from these local matched-tree campaigns?
- Is it worth adding a lighter-weight scripted benchmark so future method changes can be checked without rerunning the full 3D batch?

## References

- `R/mvSIMMAP.r`
- `tests/testthat/test-mvSIMMAP.R`
- `tests/experimental_mvsimmap_recovery_benchmark.R`
- `man/mvSIMMAP.Rd`
- `man/parameters.Rd`
- `specs/020-validation-and-open-questions.md`
- `specs/050-mvgls-simmap-mixed-design.md`
