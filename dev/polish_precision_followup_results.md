# Polish precision and diagnostics follow-up

Date: 2026-11-03 (America/Los_Angeles). Branch: `convergence-simple`; follows `e9b930c`.
The fixes and validation are in `/private/tmp/causalOT-convergence-simple`. The original checkout, saved evaluation checkpoints, and `phase1-convergence` worktree were not edited. No simulation grid or full test suite was run; no push was made.

## Fixes

- `R/cotNewtonPolish.R`: use tensor-left scalar multiplication in the blocked softmin, cold cross/source-self updates, target-self updates, warm weighted log terms, Hessian scaling, and diagonal ridge scaling. This avoids R/torch numeric-left dispatch rounding non-float32-representable scalars before multiplication. Embedded TorchScript scalar arithmetic is unchanged. Tolerances, iteration budgets, warm-attempt policy, stopping rules, and line-search acceptance rules are unchanged.
- `R/cotOOP.R`: when no polished weights are installed, keep the complete `info` record and append the convergence decision using `fw_gap_start` and `divergence_start`. `already_optimal` remains visible when a converged polish accepts zero steps; timings, `n_accepted`, `nnz`, and other diagnostics remain available.
- A rejected dtype round-trip restores the original weights and uses the ordinary cold-start tolerance outcome (`ok` or `gap above tolerance`). The rejection detail is stored in `rejection_reason`. It is not converted to `polish error: ...`. Genuine evaluator exceptions still fall back to the phase-1 rule, retain any available refinement diagnostics, and report `polish error: ...` only when that rule fails. A separate regression test distinguishes these errors from ordinary rejection.
- `tests/testthat/test-polish-precision.R` forces blocked costs at lambda `0.22916026113137714`, compares against independent R softmin/residual/Hessian formulas and dense evaluation, and compares full warm, cold-only, and dense polish results.
- `tests/testthat/test-convergence-simple.R` now asserts full diagnostic retention for a zero-step optimal polish, preservation of the original non-installed polish status, and both passing and failing cold-start outcomes after an actual float32 installation/rejection/restore path. `test-newton-polish.R` checks `rejection_reason` for rejected installation.
- `NEWS.md` describes these fixes. `.gitignore` admits the bounded replay harness, this report, and its two small summary CSVs.

## Focused test results

| File | Assertions passed | Failures | Errors | Warnings | Skips |
|---|---:|---:|---:|---:|---:|
| polish-precision | 30 | 0 | 0 | 0 | 0 |
| convergence-simple | 126 | 0 | 0 | 0 | 0 |
| newton-followup | 113 | 0 | 0 | 0 | 0 |
| newton-polish | 154 | 0 | 0 | 0 | 0 |

**423 passing assertions.** All four files were run individually using `devtools::test(filter = ...)`. Precision and new outcome tests ran on CPU; the existing follow-up and polish files ran outside the sandbox with MPS available. The existing n=256 CPU float64/float32 checks inside `test-newton-polish.R` also pass. These counts include only files rerun for this follow-up; historical counts for other files are not presented as fresh validation.

An independent pre/post check loaded the original softmin from `e9b930c` into an isolated R environment: its maximum error against an R double-precision reference was `2.141803e-8`, while the corrected function had zero error for that fixture. The new regression requires error below `1e-14`, so it fails on the original code.

## Saved n=1000 replay

Read-only inputs were the existing CPU float64 n=1000, d=3, lambda-multiplier-0.05 `phase1.rds` and `dev/newton_warm_diagnostic_results/diagnostic.rds` from the original checkout. The replay uses the saved candidate and incoming potentials directly; no phase-1 refit or full n=1000 optimization is performed. Corrected evaluations are written only into the new directory `dev/polish_precision_followup_results` in this worktree, separate from all original results.

Warm evaluations use a 256-iteration cap; the fresh cold/target-self evaluations use a 1000-iteration cap. The marginal and self tolerances remain `1e-13`.

| Initialization / weights | Iterations | Cross marginal residual | Self residual | Absolute loss difference from corrected cold candidate |
|---|---:|---:|---:|---:|
| previous_new_weights | 128 | 9.39574879e-14 | 2.22044605e-16 | 0 |
| old_cold_new_weights | 44 | 8.3751925e-14 | 2.22044605e-16 | 0 |
| corrected_cold_new_weights | 1 | 6.8704561e-14 | 3.86357613e-14 | 5.55111512e-17 |
| previous_unchanged_weights | 53 | 9.06397955e-14 | 2.22044605e-16 | — |

Both saved initializations converge to candidate loss **0.440211071933707**. Their losses agree with the fresh corrected cold evaluation to the printed precision, and their gradients differ from it by at most `8.442136e-13`. Starting from the corrected cold potentials takes one iteration. Unchanged weights also converge under the same tolerance; that row has its own objective, so no candidate-loss comparison is made. The original diagnostic reported exhaustion at 10,000 warm iterations under the affected arithmetic.

This verifies the precision correction on the saved case and full polish agreement on the small blocked fixture. It does not measure a full corrected simulation grid or establish overall benchmark speedups.

## Reproduction

From the branch worktree:

```sh
Rscript dev/convergence_simple_validate.R tests default polish-precision convergence-simple newton-followup newton-polish
Rscript dev/polish_precision_followup_validate.R "/Users/eifer/GoogleDrive/Research/Jose Zubizarreta Research/causalOT" /private/tmp/polish-precision-new-results
```

The test command needs access outside the sandbox to exercise MPS. The replay is CPU-only and refuses to reuse an existing output directory. `dev/polish_precision_followup_summary.csv` and `dev/polish_precision_followup_tests.csv` contain the recorded measurements and fresh focused test counts. `git diff --check` passes.
