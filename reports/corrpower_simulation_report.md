# Experimental BMM Simulation Report

Date: 2026-03-24

## Executive Summary

The simulation evidence now supports a clear working conclusion:

- The original `corrstrength` parameterization was scientifically promising but practically unstable, especially in higher-dimensional and weak-signal settings.
- The new `corrpower` parameterization improved robustness by reducing catastrophic failures.
- A mild stabilization pass on `corrpower` made it clearly better than `corrstrength` on the hard comparison set.
- With stabilized `corrpower`, adding fossil terminals improved recovery of the true regime covariance structure in about 70% of fossil scenarios and modestly reduced both pathological-scale fits and boundary behavior.

At this point, the best-supported experimental option is stabilized `corrpower`, not `corrstrength`.

## Simulation Campaigns Covered

This report synthesizes the main simulation campaigns run during development:

| Campaign | Scope | Main purpose |
| --- | --- | --- |
| Pilot extant-only `corrstrength` screen | Local pilot screen, mainly `p = 2`, with attempts to expand to `p = 4` | Establish whether the new correlation-structure model was viable at all |
| `corrstrength` fossil benchmark | 384 scenarios, 3840 paired extant/fossil fits | Test whether fossil terminals improve identifiability under `corrstrength` |
| Unstabilized `corrpower` vs `corrstrength` comparison | 160 hard-case runs | Test whether `corrpower` is a better model family than `corrstrength` |
| Stabilized `corrpower` vs `corrstrength` comparison | 160 hard-case runs | Test whether weak regularization materially improves `corrpower` |
| Stabilized `corrpower` fossil recovery grid | 6720 paired fits, 5760 fossil scenarios | Test whether fossils improve truth recovery under the stabilized `corrpower` fitter |

## 1. Pilot Extant-Only `corrstrength` Screen

The first extant-only `corrstrength` pilot established that the model was not obviously unusable, but it was already showing practical fragility.

From the recorded local pilot summary:

- success rate: `1.000`
- warning rate: `0.000`
- convergence warning rate: `0.000`
- median derived-regime covariance error: `0.300`
- share with covariance error `< 0.25`: `0.431`
- pathological-scale flag rate: `0.139`

The important lesson was that the model could fit and recover some signal, but pathology remained common enough that it was not ready for broader use. Attempts to extend the screen into harder `p = 4` settings also bogged down or became unstable. That early result motivated a rethink of the parameterization rather than more feature work.

## 2. `corrstrength` Fossil Benchmark

The first major fossil benchmark showed that fossils were helping, but not enough to rescue `corrstrength`.

Summary from the archived analysis bundle:

- rows: `3840`
- unique scenarios: `384`
- extant success: `96.4%`
- fossil success: `96.3%`
- timeout rate: `3.6%`
- fossils better on covariance: `60.4%`
- fossils better on absolute correlation: `52.6%`
- fossils better on mean rate: `55.2%`
- median covariance error: `0.352 -> 0.320` for extant-only to fossil+extant
- median covariance improvement: `0.024`
- pathological-scale rate: `5.9% -> 9.7%`
- boundary rate: `52.9% -> 46.8%`

Interpretation:

- Fossils were already improving covariance recovery on average.
- More fossils helped more.
- Weak-signal scenarios benefited more than strong-signal scenarios.
- But `corrstrength` still remained too unstable, especially in higher-dimensional settings (`p = 4` was the bad corner).

This benchmark justified continuing to take fossils seriously, but it also showed that the model architecture itself still needed work.

## 3. Unstabilized `corrpower` vs `corrstrength`

The next major question was whether `corrpower` was genuinely better than `corrstrength` in the hardest part of the design space.

On the first targeted comparison set (`160` hard-case runs):

- `corrpower` better on regime-`B` covariance recovery: `65.6%`
- `corrpower` better on regime-summary recovery: `64.4%`
- `corrpower` better on log-likelihood: `57.5%`
- median regime-`B` covariance error:
  - `corrpower`: `1.324861`
  - `corrstrength`: `1.605179`
- median regime-summary error:
  - `corrpower`: `0.227341`
  - `corrstrength`: `0.279938`
- pathological rate:
  - `corrpower`: `0.0%`
  - `corrstrength`: `25.0%`
- boundary rate:
  - `corrpower`: `60.6%`
  - `corrstrength`: `60.0%`

Clean subset, excluding pathological runs in either model:

- `corrpower` covariance win rate: `54.2%`
- median regime-`B` covariance error:
  - `corrpower`: `1.212830`
  - `corrstrength`: `1.206407`

Interpretation:

- `corrpower` was clearly better at avoiding catastrophic failures.
- However, it was not yet a decisive scientific win once both models were behaving.
- The main remaining defect was still weak identifiability, showing up as high boundary behavior.

## 4. Stabilized `corrpower` vs `corrstrength`

The stabilization pass added two important changes to `corrpower`:

- optimize non-reference `scale` on a log scale rather than as a squared parameter
- add weak shrinkage toward `scale = 1` and `corr_power = 1`

On the same targeted hard-case comparison set (`160` runs), the stabilized `corrpower` model performed much better:

- `corrpower` better on regime-`B` covariance recovery: `75.0%`
- `corrpower` better on regime-summary recovery: `81.9%`
- `corrpower` better on log-likelihood: `39.4%`
- median regime-`B` covariance error:
  - `corrpower`: `1.220909`
  - `corrstrength`: `1.605179`
- median regime-summary error:
  - `corrpower`: `0.214712`
  - `corrstrength`: `0.279938`
- pathological rate:
  - `corrpower`: `0.0%`
  - `corrstrength`: `25.0%`
- boundary rate:
  - `corrpower`: `30.6%`
  - `corrstrength`: `60.0%`

Clean subset:

- `corrpower` covariance win rate: `66.7%`
- median regime-`B` covariance error:
  - `corrpower`: `1.144657`
  - `corrstrength`: `1.206407`

Interpretation:

- The stabilization pass materially improved `corrpower`.
- The model no longer just avoided blowups; it also recovered the target better more often even when both models were otherwise well-behaved.
- The drop in unpenalized log-likelihood wins was expected and acceptable because the model was now trading a little fit for much better truth recovery and stability.

This was the decisive result that made stabilized `corrpower` the preferred experimental option.

## 5. Stabilized `corrpower` Fossil Recovery Grid

The final and most informative campaign asked the practical question: once the fitter is stabilized, do fossil terminals genuinely help recover the truth?

The finished fossil-recovery campaign produced:

- total rows: `6720`
- control rows (`fossil_fraction = 0`): `960`
- actual fossil scenarios: `5760`
- fit success overall: `100%`
- fit success on fossil scenarios: `100%`

### Overall fossil effect

Across the fossil scenarios only:

- fossils better on regime-`B` covariance recovery: `70.5%`
- fossils better on regime-summary recovery: `70.2%`
- median regime-`B` covariance error:
  - fossil+extant: `0.179`
  - extant-only: `0.215`
- median covariance improvement (`extant - fossil`): `0.029`
- mean covariance improvement: `0.042`
- median regime-summary RMSE:
  - fossil+extant: `0.120`
  - extant-only: `0.148`
- median summary RMSE improvement: `0.022`
- median regime-`B` mean-rate absolute error:
  - fossil+extant: `0.134`
  - extant-only: `0.163`
- median regime-`B` mean absolute correlation error:
  - fossil+extant: `0.050`
  - extant-only: `0.057`
- pathological-scale rate:
  - fossil+extant: `0.1%`
  - extant-only: `0.3%`
- boundary rate:
  - fossil+extant: `15.9%`
  - extant-only: `18.7%`

Sanity check:

- the `fossil_fraction = 0` control rows showed essentially zero median difference between the “full” and extant-only fits

That means the harness is behaving as intended and the fossil advantage is not an artifact of the bookkeeping.

### Effect of fossil fraction

More fossils helped more, very cleanly:

| Fossil Fraction | Fossils Better on Covariance | Fossils Better on Summary | Median Covariance Improvement |
| --- | ---: | ---: | ---: |
| `0.25` | `63.6%` | `63.7%` | `0.015` |
| `0.50` | `71.4%` | `70.3%` | `0.031` |
| `0.75` | `76.5%` | `76.6%` | `0.045` |

### Effect of fossil depth

Depth mattered less than expected in the current design:

| Fossil Depth | Fossils Better on Covariance | Fossils Better on Summary | Median Covariance Improvement |
| --- | ---: | ---: | ---: |
| `deep` | `70.3%` | `71.3%` | `0.025` |
| `shallow` | `70.7%` | `69.1%` | `0.033` |

The current simulation design therefore supports a general fossil effect more strongly than a strong depth-specific effect.

### Effect of extant sample size

| `n_extant` | Fossils Better on Covariance | Median Covariance Improvement | Full Boundary Rate | Extant Boundary Rate |
| --- | ---: | ---: | ---: | ---: |
| `50` | `69.4%` | `0.043` | `20.8%` | `24.2%` |
| `100` | `70.2%` | `0.029` | `15.7%` | `18.2%` |
| `200` | `72.0%` | `0.023` | `11.4%` | `13.8%` |

Interpretation:

- fossils remained useful at all three sample sizes
- absolute improvement was largest at smaller extant sample sizes
- as expected, the overall problem gets easier as `n_extant` increases

### Effect of dimensionality

| Traits (`p`) | Fossils Better on Covariance | Fossils Better on Summary | Median Covariance Improvement |
| --- | ---: | ---: | ---: |
| `2` | `66.3%` | `70.0%` | `0.023` |
| `4` | `74.7%` | `70.3%` | `0.033` |

This is encouraging. The covariance benefit from fossils was actually larger in the harder `p = 4` setting.

### Other design factors

The fossil advantage was stable across the remaining axes:

- base signal `0.2`: fossils better on covariance `70.3%`
- base signal `0.5`: fossils better on covariance `70.7%`
- mapped branch fraction target `0.2`: fossils better on covariance `67.4%`
- mapped branch fraction target `0.8`: fossils better on covariance `73.6%`
- `relation = stronger`: fossils better on covariance `70.7%`
- `relation = weaker`: fossils better on covariance `70.3%`

The low-signal runs still had more boundary behavior overall, but fossils helped there too.

## Cross-Campaign Synthesis

Taken together, the simulation record tells a coherent story:

1. The original `corrstrength` model captured a real signal but remained too unstable to trust broadly.
2. Fossils already helped under `corrstrength`, which argued that the identifiability problem was not just conceptual; it was partly a data-geometry problem.
3. The `corrpower` reparameterization reduced catastrophic failures but was initially still too boundary-prone.
4. A mild stabilization pass turned `corrpower` into a clearly better model than `corrstrength` on the hard comparison set.
5. Under stabilized `corrpower`, fossil terminals provide consistent and meaningful gains in recovery of the true regime covariance structure.

The strongest global conclusion is therefore:

> Stabilized `corrpower` plus fossil terminals is the first version of this experimental BMM correlation model that looks scientifically promising rather than merely technically interesting.

## Practical Recommendations

### Recommended current experimental default

If the experimental BMM correlation structure is to be used going forward, it should be:

- model family: `corrpower`
- fitter: stabilized `corrpower`
- interpretive target: regime covariance matrices and derived regime summaries, not raw `corr_power` alone

### What to do next

The best next step is not a new model family yet. It is to consolidate and communicate what now works:

1. Treat stabilized `corrpower` as the preferred experimental option.
2. Preserve the full simulation record and attach this report to the branch.
3. Produce a compact figure set from the fossil-recovery grid.
4. If more simulation work is needed, make it targeted:
   - test where the fossil benefit saturates
   - test robustness to incomplete fossil sampling
   - test whether fossil placement can be optimized

### What not to overstate

- `corrpower` is improved, not finished.
- Boundary behavior is still present, even though it is much better than before.
- Raw likelihood is not the right score when comparing fossil+extant fits to extant-only fits, because those are different datasets.
- These runs used complete simulated fossil data; missing data and fossil measurement error remain a separate challenge.

## Data Sources

This report was assembled from the following simulation outputs and summaries:

- Archon analysis summary for the targeted corrpower comparison campaign
- Combined CSV for targeted corrpower comparison campaign, run 1
- Combined CSV for targeted corrpower comparison campaign, stabilized run 2
- Combined CSV for corrpower fossil-recovery campaign, run 3
