# mvMORPH
mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data    

## Experimental Fork Notice

This repository is an experimental fork of `mvMORPH` where I am testing new models and features built with Codex.

Nothing in this fork is guaranteed at this point: APIs may change, behavior may shift, numerical results may need additional validation, and features may be incomplete or removed without notice.

If you need the stable `mvMORPH` package, use the CRAN release or the upstream repository rather than relying on this fork.

## Experimental Additions In This Fork

This fork currently includes several experimental additions on top of upstream `mvMORPH`:

- `mvgls(..., model="OUM")` has been extended so user-supplied covariates can be retained alongside regime-specific optima rather than being collapsed away.
- Experimental correlation-aware `BMM` support is available through `mvgls(..., model="BMM", bmm.structure="corrpower")`, with a simpler nested comparison mode via `bmm.scale=FALSE`.
- A standalone experimental `mvSIMMAP()` fitter can fit painted `simmap` trees with mixed `BM`, `OU`, `OUM`, and `EB` regimes, including shared `process.groups`, richer summaries, and improved post-fit diagnostics.
- This fork also includes simulation harnesses, reports, and design notes that were added to validate and iterate on these experimental paths.

For fuller branch-local summaries of the implemented features, see `specs/010-fork-delta.md` and `specs/020-validation-and-open-questions.md`. For the draft design note covering planned `mvgls` integration of mixed-process SIMMAP support, see `specs/050-mvgls-simmap-mixed-design.md`.

This package allows the fitting of multivariate evolutionary models (Ornstein-Uhlenbeck, Brownian motion, Early burst, Shift models) on species trees and time series.
It also provides functions to compute log-likelihood of users specified models with fast methods (*e.g.*, for Bayesian approaches or customized comparative methods), simulates correlated traits under various models, constrain various parts of multivariate models...

The package implement now efficient methods for high-dimensional multivariate comparative methods (mvgls) based on Penalized likelihood as well as associated tests (Wilks, Pillai...)

The package is designed to handle ultrametric and non-ultrametric trees (*i.e.* with fossil species) and missing data in multivariate datasets (NA values), SIMMAP mapping of discrete traits, measurement error, etc...

See the package vignettes for details and examples: browseVignettes("mvMORPH").

## Current Status

This experimental fork currently reports package version `1.2.2` in `DESCRIPTION`.

For a stable release of `mvMORPH`, use the CRAN package or the upstream repository rather than this experimental fork.

## **Package Installation**

For stable use, install the CRAN or upstream `mvMORPH` release.

If you want to try this experimental fork, the safest approach is to install it into an isolated library or project rather than into your default R library. This fork still installs as package `mvMORPH`, so a normal install can overwrite an existing user-level `mvMORPH` installation.

### Stable upstream install

The command below installs the upstream GitHub package, not this experimental fork:

```
library(remotes)

install_github("JClavel/mvMORPH", build_vignettes = TRUE)

```

### Experimental fork safe install

The example below installs this fork into a dedicated library that you can load explicitly when you want the experimental version:

```
library(remotes)

exp_lib <- path.expand("~/Library/R/experimental-mvMORPH")
dir.create(exp_lib, recursive = TRUE, showWarnings = FALSE)

install_github("jakeberv/mvMORPH-berv", build_vignettes = TRUE, lib = exp_lib)

library(mvMORPH, lib.loc = exp_lib)
```

If you want to keep using your usual packages while preferring this experimental `mvMORPH`, prepend the experimental library to `.libPaths()` instead of replacing your normal library search path:

```
exp_lib <- path.expand("~/Library/R/experimental-mvMORPH")
.libPaths(c(exp_lib, .libPaths()))

library(mvMORPH)
library(ape)
library(phytools)
```

This setup makes R look in the experimental library first for `mvMORPH`, then fall back to your usual user and Homebrew libraries for everything else.

Only load one `mvMORPH` installation per R session. If you want to keep this fork fully isolated inside a project, using `renv` is also a good option.

### Note for macOS + Homebrew R

On macOS with Homebrew-managed R, Homebrew installs the R executable and system libraries, but user-installed packages are typically placed in your personal R library. If you install this fork without specifying `lib=`, R will usually install `mvMORPH` into the first writable path in `.libPaths()`, which can replace your usual user-level `mvMORPH` package without affecting the Homebrew R installation itself.

Using a dedicated `lib=` path, as shown above, keeps this experimental fork separate from the `mvMORPH` package you may already use in your default library.


(The installation may crash if your dependencies are not up to date. Note that you may also need to install Rtools to compile the C codes included in the package. For [Windows] (https://cran.r-project.org/bin/windows/Rtools/) and for [Mac] (https://mac.r-project.org/) (and [Tools] (https://mac.r-project.org/tools/) )

## **Report an issue**
Bugs encountered when using this experimental fork can be reported [here](https://github.com/jakeberv/mvMORPH-berv/issues)

## **Package citation**

**Clavel, J., Escarguel, G., Merceron, G. 2015.** mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.
