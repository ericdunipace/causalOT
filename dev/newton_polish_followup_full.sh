#!/usr/bin/env bash
# Run from the repository root. Long-running validations, intentionally deferred.
# Before reconsidering the evidence-based warm_start default, compare warm and
# cold MULTI-ITERATION total time at n = 5000. The one-iteration timing below
# cannot establish an end-to-end speedup. That additional comparison is deferred;
# it is intentionally not launched by the present bug-fix validation.
set -euo pipefail
POLISH_FOLLOWUP_NS=256,1000,2000,5000 Rscript dev/newton_polish_followup_validate.R
Rscript dev/newton_polish_followup_fixtures.R
Rscript dev/polish_default_validate.R
for mode in none auto; do
  CONV_TEST_POLISH="$mode" \
  CONV_TEST_OUTPUT="dev/newton_followup_full_${mode}.rds" \
  CONV_TEST_CHECKPOINT_DIR="dev/newton_followup_full_${mode}_batches" \
  CONV_TEST_PLOT_PREFIX="dev/newton_followup_full_${mode}" \
  Rscript dev/convergence_lambda_grid.R
done
