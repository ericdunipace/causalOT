# Convergence simplification validation

Date: 2026-11-02 (America/Los_Angeles). Branch: `convergence-simple`, based on `74fa9c9`.
Implementation and validation used `/private/tmp/causalOT-convergence-simple`. The original checkout and the `phase1-convergence` worktree were not edited. No push was made.

## Result

Final outcomes are `converged` and `reason`. Eligible fits use `gap <= atol + tol * max(0, D)` with the user's solver `tol` in both phases. Unpolished fits also require normal stopping and use `atol = ot_tol`. Cold polished fits use `64 * eps64 * scale` or `eps32 * scale`, depending on returned weight precision. Rejected polish judges the retained starting point using its cold gap and divergence. Numerical failures still store NA weights and restart the path from the previous usable point.

Stopping rules, the line search, `polish.control$warm_start`, and polish internals were not changed. The internal polish classifier and its precision floors remain in place. No phase-1 confirmation solve, `phase1_certify`, `cot_sinkhorn_confirm`, new `marginal_tol`, or `phase1_y` cache was added.

## Files changed

- `R/cotNewtonPolish.R`: add the finite-input `cot_converged` helper; retain internal polish stopping.
- `R/cotOOP.R`: apply the common test, replace fit statuses with reasons, handle cold rejected polish, and update selection messages.
- `R/cotCCPTrain.R`: propagate reasons and rename the legacy stall predicate without changing its logic.
- `R/cotClass.R`, `man/cotOptions.Rd`, `man/OTProblem.Rd`, `NEWS.md`: document tolerances, outcomes, and numeric diagnostics.
- `tests/testthat/test-convergence-simple.R`: cover absolute and relative tolerance, zero/negative divergence, non-finite inputs, iteration caps, zero-divergence solves, excessive polished gaps, rejected polish, polish errors, and out-of-scope flags.
- `tests/testthat/test-newton-followup.R`, `test-newton-polish.R`, `test-lambda-selector.R`, `test-cotOOP.R`: update outcome and threshold assertions.
- `.gitignore`, `dev/convergence_simple_validate.R`, this report, and four `dev/convergence_simple_summary_*.csv` files: reproducible focused validation and recorded results.

## Focused tests

Each file was run separately through `devtools::test(filter = ...)`, outside the sandbox with MPS available and `NOT_CRAN=true`. Only the listed files were run; the full suite and the known failing `test-dual_opts.R` were not run. Final results supersede intermediate runs while updating old assertions.

| File | Passed assertions | Failures | Errors | Warnings | Skips |
|---|---:|---:|---:|---:|---:|
| ccp-solver | 36 | 0 | 0 | 0 | 0 |
| convergence-simple | 75 | 0 | 0 | 0 | 0 |
| cotOOP | 116 | 0 | 0 | 0 | 0 |
| lambda-selector | 135 | 0 | 0 | 0 | 0 |
| newton-followup | 113 | 0 | 0 | 0 | 0 |
| newton-polish | 154 | 0 | 0 | 0 | 0 |

Total: **629 passing assertions**, no failures, errors, warnings, or skips.

`devtools::document()` regenerated the docs. The installed newer roxygen also rewrote unrelated topics and metadata; those unrelated rewrites were discarded. It reported pre-existing S3 export and Measure documentation warnings. Both retained Rd files parse successfully. `git diff --check` passes.

## n = 256 regression: actual measurements

Seed 1; 256 observations, 3 covariates; 134 source observations; lambda is the multiplier times the median squared pairwise distance. Every phase-1 fit stopped normally after **3 iterations**. No-polish cases retain 134 positive weights and do not converge. Auto-polished cases converge and retain **39 / 34 / 31** positive weights at multipliers **1 / 5 / 25**.

The measured cost scale is about **46.9974**, rather than the illustrative 31 in the plan. Accordingly, float32 `atol` is about **5.60253e-6** and float64 `atol` about **6.67873e-13**. Float32 calibration passes on CPU and MPS.

| Device / dtype | Lambda multiplier | Polish | Gap | D | atol | atol + tol * max(0,D) | Converged | Positive weights |
|---|---:|---|---:|---:|---:|---:|---|---:|
| CPU float64 | 1 | none | 0.0575325784 | 0.00320381652 | 0.001 | 0.00100032038 | FALSE | 134 |
| CPU float64 | 1 | auto | 1.99396055e-13 | 0.000169404824 | 6.67873281e-13 | 1.69411502e-08 | TRUE | 39 |
| CPU float64 | 5 | none | 0.0295151928 | 0.000334235732 | 0.001 | 0.00100003342 | FALSE | 134 |
| CPU float64 | 5 | auto | 4.05169231e-11 | 9.457221e-07 | 6.67873281e-13 | 9.52400833e-11 | TRUE | 34 |
| CPU float64 | 25 | none | 0.0176285751 | 0.000203875357 | 0.001 | 0.00100002039 | FALSE | 134 |
| CPU float64 | 25 | auto | 5.41122702e-13 | 6.91485047e-09 | 6.67873281e-13 | 1.35935833e-12 | TRUE | 31 |
| CPU float32 | 1 | none | 0.0575318336 | 0.00320398808 | 0.001 | 0.0010003204 | FALSE | 134 |
| CPU float32 | 1 | auto | 4.95067718e-08 | 0.000169404824 | 5.60252738e-06 | 5.61946786e-06 | TRUE | 39 |
| CPU float32 | 5 | none | 0.0295163393 | 0.000334024429 | 0.001 | 0.0010000334 | FALSE | 134 |
| CPU float32 | 5 | auto | 4.33319205e-08 | 9.45722098e-07 | 5.60252738e-06 | 5.60262195e-06 | TRUE | 34 |
| CPU float32 | 25 | none | 0.017649889 | 0.000196814537 | 0.001 | 0.00100001968 | FALSE | 134 |
| CPU float32 | 25 | auto | 8.72450994e-08 | 6.91485758e-09 | 5.60252738e-06 | 5.60252807e-06 | TRUE | 31 |
| MPS float32 | 1 | none | 0.0575327873 | 0.00320374966 | 0.001 | 0.00100032037 | FALSE | 134 |
| MPS float32 | 1 | auto | 8.24346633e-08 | 0.000169404824 | 5.60252738e-06 | 5.61946786e-06 | TRUE | 39 |
| MPS float32 | 5 | none | 0.0295159817 | 0.000334620476 | 0.001 | 0.00100003346 | FALSE | 134 |
| MPS float32 | 5 | auto | 7.99563096e-08 | 9.45722098e-07 | 5.60252738e-06 | 5.60262195e-06 | TRUE | 34 |
| MPS float32 | 25 | none | 0.0176448822 | 0.000193834305 | 0.001 | 0.00100001938 | FALSE | 134 |
| MPS float32 | 25 | auto | 7.03276763e-08 | 6.91485624e-09 | 5.60252738e-06 | 5.60252807e-06 | TRUE | 31 |

## Base equivalence and numerical limits

- The six small cases (mirror, LBFGS, and accelerated CCP, each on default MPS and CPU float64) match `74fa9c9` bit for bit in iterations, losses, weights, evaluation outputs, and potentials.
- All 18 regression cases match `74fa9c9` bit for bit in the full phase-1 result, including iteration count, loss, stopping flags, FW gap, and divergence. All nine unpolished weight vectors are bit-for-bit identical.
- All polished supports match. Seven of nine polished weight vectors are also bit-for-bit identical. At multiplier 1, CPU float64 differs by at most `3.381492e-10` and MPS float32 by at most `6.984919e-10`; each loss differs by only `2.220446e-16`, below the float64 evaluator floor (`6.678733e-13`).
- A fresh repeat of the branch's CPU float64 multiplier-1 polished case differs from the first branch run by `4.570406e-10` in weights and `6.661338e-16` in loss, with identical support. Its measured gap is `5.140790e-10`, still below its final threshold `1.694115e-8`. This demonstrates run-to-run numerical variation; polished bitwise identity is not claimed.
- The comparison harness requires exact phase-1 results and exact unpolished weights, and verifies polished supports and loss agreement within evaluator roundoff. It also reports every polished weight difference. The initial overly strict comparison of polished weights failed, leading to the repeat check above.

## Reproduction

Create a separate base snapshot using `git archive 74fa9c9`, extracted into `/private/tmp/causalOT-convergence-base`. Do not load the user's live checkout as a baseline. From this branch worktree, run the following outside the sandbox for MPS access. The `equivalence` command includes the n = 256 regression and requires the user's approval before execution.

```sh
Rscript dev/convergence_simple_validate.R tests default convergence-simple newton-followup newton-polish lambda-selector ccp-solver cotOOP
Rscript dev/convergence_simple_validate.R equivalence /private/tmp/causalOT-convergence-base
```

The final comparison was rerun directly after correcting its polished-point criterion; both the phase-1 and regression comparisons pass. Raw logs and RDS snapshots remain in the worktree under ignored `dev/convergence_simple_*` paths. The four tracked summary CSVs contain the final test counts, all 18 measured regression rows, six small equivalence cases, and full regression equivalence diagnostics.
