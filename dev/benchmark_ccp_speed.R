library(causalOT)

set.seed(20260708)
n <- 200L
m <- 200L
d <- 2L
lambda <- 100.0
x <- matrix(stats::rnorm(n * d), n, d)
y <- matrix(stats::rnorm(m * d), m, d)

run_cot <- function(opt) {
  cot <- causalOT:::COT$new(source = x, target = y, options = opt)
  elapsed <- system.time(cot$solve())[["elapsed"]]
  list(elapsed = elapsed, weights = cot$weights)
}

base_ccp <- cotOptions(
  lambda = lambda,
  opt.direction = "primal",
  solver = "mirror",
  cost.online = "tensorized",
  niter = 50L,
  ot_niter = 50L,
  ot_tol = 1e-5,
  tol = 1e-6
)

dual <- cotOptions(
  lambda = lambda,
  opt.direction = "dual",
  cost.online = "tensorized",
  niter = 50L,
  ot_niter = 50L,
  ot_tol = 1e-5,
  tol = 1e-6
)

res <- list(
  ccp_single_tau = run_cot(base_ccp),
  dual = run_cot(dual)
)

ccp_sparse <- causalOT:::COT$new(source = x, target = y, options = base_ccp)
ccp_sparse_priv <- ccp_sparse$.__enclos_env__$private$optimizer$.__enclos_env__$private
ccp_sparse_priv$ccp_diagnostics_every <- 25L
res$ccp_sparse_diagnostics <- list(
  elapsed = system.time(ccp_sparse$solve())[["elapsed"]],
  weights = ccp_sparse$weights
)

ccp_backtracking <- causalOT:::COT$new(
  source = x,
  target = y,
  options = base_ccp
)
ccp_backtracking_priv <- ccp_backtracking$.__enclos_env__$private$optimizer$.__enclos_env__$private
ccp_backtracking_priv$ccp_backtracking <- TRUE
res$ccp_backtracking <- list(
  elapsed = system.time(ccp_backtracking$solve())[["elapsed"]],
  weights = ccp_backtracking$weights
)

timings <- vapply(res, `[[`, numeric(1L), "elapsed")
print(timings)

reference <- res$dual$weights
weight_error <- vapply(
  res,
  function(x) max(abs(x$weights - reference)),
  numeric(1L)
)
print(weight_error)
