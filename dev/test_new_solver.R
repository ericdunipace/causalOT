devtools::load_all(".")

set.seed(2141)
n <- 15
m <- 13
d <- 2
lambda <- 100

x <- matrix(rnorm(n * d), n, d)
y <- matrix(rnorm(m * d), m, d)

opt <- cotOptions(
  lambda = lambda,
  opt.direction = "primal",
  solver = "mirror",
  niter = 1000L,
  ot_niter = 5L,
  ot_tol = 1e-8,
  tol = 1e-10,
  device = torch::torch_device("cpu"),
  dtype = torch::torch_double()
)

cot <- causalOT:::COT$new(source = x, target = y, options = opt)
cot$solve()

sum(cot$weights)
range(cot$weights)
head(cot$weights)

ots <- cot$.__enclos_env__$private$optimizer$.__enclos_env__$private$ot_objects
ot <- ots[[ls(ots)[[1L]]]]
causalOT:::cot_ccp_diagnostics(ot)


opt2 <- cotOptions(
  lambda = lambda,
  opt.direction = "primal",
  cost.online = "tensorized",
  niter = 1000,
  ot_niter = 5,
  ot_tol = 1e-8,
  tol = 1e-10,
  line_search_fn = "strong_wolfe",
  device = torch::torch_device("cpu"),
  dtype = torch::torch_double()
)

cot2 <- causalOT:::COT$new(source = x, target = y, options = opt2)
suppressWarnings(cot2$solve())

pkgload::load_all(".")
opt3 <- cotOptions(
  lambda = lambda,
  opt.direction = "dual",
  cost.online = "tensorized",
  niter = 1000,
  ot_niter = 5,
  ot_tol = 1e-8,
  tol = 1e-8,
  line_search_fn = "strong_wolfe",
  device = torch::torch_device("cpu"),
  dtype = torch::torch_double(),
  dual.solver = "fixed_point_joint"
  # dual.solver = "joint"
)

cot3 <- causalOT:::COT$new(source = x, target = y, options = opt3)
# debugonce(
#   cot3$.__enclos_env__$private$optimizer$.__enclos_env__$private$optimization_loop
# )
cot3$solve(
  debug_joint = TRUE,
  fixed_point_eta = 2,
  fixed_point_beta = 0.5,
  fixed_point_weight_floor = 1e-6,
  fixed_point_objective_check_every = 50,
  fixed_point_eta_anneal_every = 5,
  fixed_point_eta_min = 1e-8
)

cbind(cot$weights, cot2$weights, cot3$weights)

c(
  ccp = ot_distance(
    x1 = x,
    x2 = y,
    a = cot$weights,
    b = rep(1 / m, m),
    p = 2,
    penalty = lambda
  ),
  primal = ot_distance(
    x1 = x,
    x2 = y,
    a = cot2$weights,
    b = rep(1 / m, m),
    p = 2,
    penalty = lambda
  ),
  dual = ot_distance(
    x1 = x,
    x2 = y,
    a = cot3$weights,
    b = rep(1 / m, m),
    p = 2,
    penalty = lambda
  )
)


system.time(causalOT:::COT$new(source = x, target = y, options = opt)$solve())
system.time(causalOT:::COT$new(source = x, target = y, options = opt2)$solve())
system.time(causalOT:::COT$new(source = x, target = y, options = opt3)$solve())


microbenchmark::microbenchmark(
  ccp = causalOT:::COT$new(source = x, target = y, options = opt)$solve(),
  dual = causalOT:::COT$new(source = x, target = y, options = opt3)$solve(),
  times = 100L
)
