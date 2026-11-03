# Newton polish bug-fix results

2026-11-02. Implements `dev/newton_polish_bugfix_plan.md`. Both `polish = "auto"` and `warm_start = TRUE` remain the defaults.

## Changes

- Float64 conversion now transfers tensors to CPU before changing dtype, avoiding unsupported MPS float64 conversion.
- Unpolished and error-fallback fits use only phase 1’s stored FW gap, full debiased divergence, and stopping flag. No extra Sinkhorn evaluation runs. Certification requires the stopping rule and the strict solver-relative threshold, with only the double-arithmetic floor. There is no loose resolution tier.
- Evaluator/cost-building errors preserve valid phase-1 weights and their diagnostics, with the error recorded. Only invalid weights or a non-finite phase-1 loss mark a fit failed. Failed fits retain NA storage, path restoration, and selector exclusion.
- Removed `polish.control$certify`. Bootstrap and delta selection continue to include all valid fits, including uncertified ones. Selected-fit warnings use the recorded reason when gap/D is unavailable.
- Polished fits compute one final certificate in `info()$polish[[lambda]][[delta]]$certification`; `status` records the refinement outcome separately. Unpolished certificates have `source = "phase1"`.
- Centralized the all-failed path check in a shared helper called by each solve entry point, retaining the selector guard. Removed the delta-level lambda comparisons. Timing tests now bound unaccounted `other` time rather than asserting a remainder-based identity.
- The deferred full-run script now explicitly requires comparing warm/cold multi-iteration total time at n = 5000 before reconsidering the warm-start default. This comparison was not run.

## Hardware and scope

The same R 4.6.1 / torch 0.17.0 installation reports MPS unavailable **inside the sandbox**, but MPS available **outside the sandbox**. The earlier CPU-only conclusion was caused by that restriction, not unavailable hardware. No runtime installation was changed.

Final validations were run outside the sandbox: `backends_mps_is_available() = TRUE`, `cuda_is_available() = FALSE`, default device **MPS**, default dtype **float32**. The requested suite set was also rerun with CPU/float64 defaults in the same unrestricted runtime. Explicitly requested dtypes/devices in individual tests were preserved. The default-accelerator test ran on MPS and passed.

## Test suites

Counts below refer to reported failures/errors and unexpected warnings; intentional warning assertions passed.

| Suite | Default device: MPS | CPU float64 defaults | Skips: MPS / CPU |
|---|---|---|---:|
| newton-followup | passed; 0 warnings | passed; 0 warnings | 0 / 1 |
| newton-polish | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| lambda-selector | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| ccp-solver | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| cotOOP | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| cotProblem | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| cotClass | passed; 0 warnings | passed; 0 warnings | 3 / 3 |
| gridSearch | passed; 0 warnings | passed; 0 warnings | 0 / 0 |
| costMethod | passed; 0 warnings | passed; 0 warnings | 2 / 2 |
| dual_opts | FAILED: 1 failure; 0 warnings | FAILED: 1 failure; 0 warnings | 3 / 3 |

Every final failure and warning:

- **Default device (MPS):** `test-dual_opts.R`, test **“training function works for dual optimizer”**, line 323: `cot_names` and `otp_names` differ in length, 10 versus 12. One failure; zero unexpected warnings.
- **CPU float64 defaults:** the same test and failure, 10 versus 12. One failure; zero unexpected warnings.
- No other final failures, errors, or unexpected warnings occurred in either run. The same dual-options failure was reproduced against an isolated HEAD checkout.

Corrected during MPS validation: `test-newton-followup.R`, **“default accelerator fits actually refine using float64 data costs”**, initially failed because it used object identity to compare two equal torch dtype objects. The assertion now uses dtype value equality. Its targeted rerun and the final full default-MPS run both passed, with no hardware skip.

Skips: three `cotClass` tests, two `costMethod` tests, and three `dual_opts` tests requiring unavailable rkeops. The CPU-default run also intentionally skips the default-accelerator test because its selected default device is CPU. The MPS run does not skip that test.

Additional assertions cover zero post-phase-1 calls to `cot_polish_eval`, `cot_polish_self`, and `cot_ccp_y_self_term` in the unpolished hook; negligible hook time; strict gap/D rejection; niter reasons; preserved weights on evaluator error; NaN exclusion/all-failed errors; bootstrap selection of an uncertified fit over a certified fit; timing overhead bounds; and removal of the certify control.

## n = 256 regression

All 18 runs passed the regression assertions. All auto fits refined and certified (strict or resolution tier). Their supports are 39/34/31. All none fits retain 134 positive weights and report uncertified, using phase-1 diagnostics without an extra solve.

`default_float32` below means device omitted and dtype explicitly float32; its actual device was **MPS**. Both CPU configurations were explicitly set to CPU.

| Configuration | Actual device | Multiplier | Polish | Certificate | Gap | Divergence | Support | Cost source |
|---|---|---:|---|---|---:|---:|---:|---|
| cpu_float64 | cpu | 1 | none | uncertified | 0.0575326 | 0.00320382 | 134 | solver_rounded |
| cpu_float64 | cpu | 1 | auto | certified | 1.97176e-13 | 0.000169405 | 39 | solver_rounded |
| cpu_float64 | cpu | 5 | none | uncertified | 0.0295152 | 0.000334236 | 134 | solver_rounded |
| cpu_float64 | cpu | 5 | auto | certified_resolution | 4.05169e-11 | 9.45722e-07 | 34 | solver_rounded |
| cpu_float64 | cpu | 25 | none | uncertified | 0.0176286 | 0.000203875 | 134 | solver_rounded |
| cpu_float64 | cpu | 25 | auto | certified_resolution | 5.41123e-13 | 6.91485e-09 | 31 | solver_rounded |
| cpu_float32 | cpu | 1 | none | uncertified | 0.0575318 | 0.00320399 | 134 | solver_rounded |
| cpu_float32 | cpu | 1 | auto | certified_resolution | 4.95068e-08 | 0.000169405 | 39 | data_float64 |
| cpu_float32 | cpu | 5 | none | uncertified | 0.0295163 | 0.000334024 | 134 | solver_rounded |
| cpu_float32 | cpu | 5 | auto | certified_resolution | 4.33319e-08 | 9.45722e-07 | 34 | data_float64 |
| cpu_float32 | cpu | 25 | none | uncertified | 0.0176499 | 0.000196815 | 134 | solver_rounded |
| cpu_float32 | cpu | 25 | auto | certified_resolution | 8.72451e-08 | 6.91486e-09 | 31 | data_float64 |
| default_float32 | mps | 1 | none | uncertified | 0.0575328 | 0.00320375 | 134 | solver_rounded |
| default_float32 | mps | 1 | auto | certified_resolution | 8.24347e-08 | 0.000169405 | 39 | data_float64 |
| default_float32 | mps | 5 | none | uncertified | 0.029516 | 0.00033462 | 134 | solver_rounded |
| default_float32 | mps | 5 | auto | certified_resolution | 7.99563e-08 | 9.45722e-07 | 34 | data_float64 |
| default_float32 | mps | 25 | none | uncertified | 0.0176449 | 0.000193834 | 134 | solver_rounded |
| default_float32 | mps | 25 | auto | certified_resolution | 7.03277e-08 | 6.91486e-09 | 31 | data_float64 |

## Phase-1 stopping equivalence

An isolated `git archive HEAD` checkout and the working tree were loaded in separate R processes. The hook was bypassed in both to compare the phase-1 solvers directly. With seed 804, n_source = 32, n_target = 48, lambda = 5, niter = 100, and tol = 1e-4, the resulting data frames of iteration counts and phase-1 losses are **bit-for-bit identical** (`identical(head, fixed)`), in both default-MPS and CPU/float64 configurations. These comparisons were rerun outside the sandbox with MPS available.

| Solver configuration | Actual device | Iterations | Phase-1 loss (rounded for display) |
|---|---|---:|---:|
| mirror | MPS | 100 | 1.50528109 |
| lbfgs | MPS | 7 | 0.0117753744 |
| ccp_accelerated | MPS | 100 | 4.45953178 |
| mirror | CPU | 100 | 1.50528132 |
| lbfgs | CPU | 11 | 0.0115634542 |
| ccp_accelerated | CPU | 100 | 3.56006475 |

The stopping conditions themselves are unchanged. LBFGS stores the existing `self$loss` used by `optimality_ok`, whose `sinkhorn_dist` subtracts both source and target self terms. CCP stores `ccp_loss - ccp_loss_offset`, where the offset is the target self term. Both reported divergences therefore include the y self term. Raw diagnostics reset at each solve to prevent stale certificates.

## Reproduction

```sh
NOT_CRAN=true Rscript dev/newton_polish_bugfix_validate.R default
NOT_CRAN=true Rscript dev/newton_polish_bugfix_validate.R cpu_float64
Rscript dev/newton_polish_bugfix_validate.R regression
```

Run these commands with Metal/MPS access, outside an OS sandbox that hides the GPU. The scripts print the actual device and backend availability and retain per-test failure/warning details in CSV/RDS files. Default-device validation must print `mps` to reproduce the MPS results.

Phase-1 comparison (HEAD checkout requires the usual built package library):

```sh
Rscript dev/newton_polish_bugfix_validate.R phase1 /path/to/HEAD /tmp/phase1-head.rds
Rscript dev/newton_polish_bugfix_validate.R phase1 . /tmp/phase1-fixed.rds
Rscript -e 'stopifnot(identical(readRDS("/tmp/phase1-head.rds"), readRDS("/tmp/phase1-fixed.rds")))'
```

No n = 5000 timing or full convergence-grid rerun was performed.
