# Newton polish follow-up results

Historical follow-up report. Certification semantics and hardware validation
are superseded by [the bug-fix report](newton_polish_bugfix_results.md).

Implemented 2026-11-02. `polish = "auto"` remains the default. The phase-1
stopping rules and `R/cotCCPTrain.R` are unchanged.

## Changes

- Shared FW-gap certificate and five fit statuses, including certification with
  refinement disabled. Failed fits store NA weights, are excluded from every
  selector, and cannot poison later warm starts. Optimization errors propagate.
- Float32 measures retain pre-cast CPU float64 points. Default and supported R
  costs are rebuilt in bounded blocks. Unsupported R costs fall back with a
  reason; custom KeOps strings remain outside scope. Final rounded weights are
  certified against the same float64 costs; device loss is reported separately.
- Cost normalization uses a dtype-matched scalar. This also fixes an existing
  p = 3 discrepancy caused by float32 rounding of the R scalar `1/p`.
- Cached cost diameters and target self terms, skipped float64 dtype rechecks,
  named stage and probe timings, and marginal-residual warm line searches.
  Warm self solves use a gauge-corrected update in the high-epsilon regime;
  low-epsilon solves retain damping. Near evaluator resolution, probes switch
  to cold evaluation. Certification is cold and existing cold results are reused.
- A.4 screening was not added: the Hessian was not the dominant cost after A.3.

## Focused validation

The n = 256 default regression passed for both solver dtypes:

| multiplier | float64 support | float32 support | divergence (approximately) |
|---|---:|---:|---:|
| 1 | 39 | 39 | 1.694048e-4 |
| 5 | 34 | 34 | 9.457221e-7 |
| 25 | 31 | 31 | 6.91486e-9 |

All auto fits certify (strict or resolution); all six unpolished fits retain
134 positive weights and report `uncertified`. Returned-weight gaps were
checked against the actual certification costs.

All 30 saved prototype fixtures (six starts at each multiplier 0.05, 0.25, 1,
5, 25) produced identical supports between warm and cold runs. The largest
absolute loss difference was 1.398881e-14, within
`1e-12 * max(1, abs(loss))`. These difficult fixtures were validated with
`sinkhorn_niter = 100000`; the default remains 10000. Some probes on these
fixtures exhaust the default budget and require a cold fallback.

At multiplier 0.05, cold probes used 420,981 update rounds in total and warm
probes used 417,302: 3,679 fewer (including discarded warm work and cold
fallbacks, but excluding final certification solves). This did not translate to lower total time
for those fixtures: 128.988 s cold versus 158.221 s warm. Warm starts do not
uniformly improve time on nearly disconnected problems.

Paired one-iteration scaling, with diameter preparation outside the timed core:

| n | cold probes: total / line search (s) | warm probes: total / line search (s) |
|---|---:|---:|
| 256 | 0.128 / 0.025 | 0.142 / 0.005 |
| 1000 | 1.254 / 0.360 | 1.439 / 0.122 |
| 2000 | 5.584 / 1.300 | 5.994 / 0.360 |

At n = 2000, line-search time fell by 72%, exceeding the 50% target, with
45 versus 8 probe update rounds. Total time for a single iteration increased:
the warm run needs a separate cold certification (1.368 s), while the cold
run reuses its accepted probe. Both arms have the same support and agree in
loss within the requested tolerance. Stage timings account for the total.
These timings are indicative single runs, not a claim of uniform end-to-end
speedup.

Tests cover certificate tiers, invalid weights/losses, NaN path exclusion and
all-failed errors, selector fallback/warnings, unsupported configurations,
cost rebuilding and custom-cost fallbacks, p = 1/2/3 block equivalence,
warm/cold agreement, cache reuse, timing accounting, and redundant solve removal.
The newton, selector, CCP, OOP, COT problem/class, grid-search, and cost suites
were exercised. KeOps-dependent tests are skipped because rkeops is unavailable.
The separate dual-options suite has an existing setup-argument signature
mismatch (10 versus 12 arguments); those signatures already differ in HEAD.
A direct dual fit confirms `not_applicable` is reported.

## Reproduction and deferred runs

Focused regression and scaling:

```sh
Rscript dev/newton_polish_followup_validate.R
Rscript dev/newton_polish_followup_fixtures.R
```

Tests (including tests skipped on CRAN):

```sh
NOT_CRAN=true Rscript -e 'pkgload::load_all(".", quiet=TRUE); for (f in c("newton-followup", "newton-polish", "lambda-selector", "ccp-solver", "cotOOP", "cotProblem", "cotClass", "gridSearch", "costMethod")) testthat::test_file(paste0("tests/testthat/test-", f, ".R"), reporter="summary")'
```

As requested, n = 5000 and the full paired convergence experiment were not run.
To run them with separate none/auto output files and checkpoints:

```sh
bash dev/newton_polish_followup_full.sh
```

The full script uses the existing local `dev/convergence_lambda_grid.R`,
`dev/polish_default_validate.R`, and saved prototype fixtures. These development
inputs remain ignored by git, as before. The new validation scripts and this
report are explicitly included in version control visibility.
