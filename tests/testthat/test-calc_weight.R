arg.names <- c("w0",  "w1", "gamma","args","estimand", "method")


testthat::test_that("works for Wass", {
  testthat::skip_on_cran()
 
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 2, 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = "osqp"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # testthat::skip_if_not_installed("gurobi")
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                          constraint = 2, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "gurobi"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 2, 
                                                           estimand = e, 
                                                           method = "Wasserstein",
                                                           solver = "mosek"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass, grid/formula", {
  testthat::skip_on_cran()
  
  set.seed(23483)
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  
  for (e in estimates) {
    # print(e)
    testthat::expect_silent(
      weight.check[[e]] <- calc_weight(data = data, 
                                       constraint = NULL,
                                       grid.search = TRUE,
                                       estimand = e, 
                                       p = 4,
                                       formula = "~.+0",
                                       balance.constraints = 0.5,
                                       method = "Wasserstein",
                                       solver = "osqp",
                                       wass.method = "greenkhorn",
                                       iter = 10, eval.method = "bootstrap")
    )
  }
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  for (e in estimates) {
    # print(e)
    testthat::expect_silent(
      weight.check[[e]] <- calc_weight(data = data, 
             constraint = NULL,
             grid.search = TRUE,
             estimand = e, 
             p = 4,
             formula = "~.+0",
             balance.constraints = 0.5,
             method = "Wasserstein",
             solver = "mosek",
             wass.method = "greenkhorn",
             iter = 10, eval.method = "bootstrap")
      )
  }
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
})

testthat::test_that("works for Wass, lbfgs", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "lbfgs"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = list(penalty = 2), 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = solver, penalty = "entropy"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                          constraint = 2, 
  #                                                          estimand = e, 
  #                                                          method = "Wasserstein",
  #                                                          solver = "cplex"))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n0,n0), 
  # #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # # testthat::expect_match(all.equal(rep(1/n1,n1), 
  # #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 2, 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = solver,
                                                       penalty = "L2"))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass, lbfgs formula", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "lbfgs"
  estimates <- c("ATT", "ATC", "ATE")
  formula <- "~.+0"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = list(penalty = 2), 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = solver, penalty = "entropy",
                                                       formula = formula,
                                                       balance.constraints = 0.1))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = list(penalty = 2), 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = solver, penalty = "L2",
                                                       formula = formula,
                                                       balance.constraints = 0.1))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass divergence", {
  testthat::skip_on_cran()
  causalOT:::skip_if_no_geomloss()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = list(penalty = 10000), 
                                                       estimand = e, p = 2,
                                                       method = "Wasserstein",
                                                       solver = "lbfgs",
                                                       search = "LBFGS",
                                                       penalty = "entropy",
                                                       add.divergence = TRUE,
                                                       stepsize = 1e-2))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                      constraint = list(penalty = 10000), 
  #                                                      estimand = e, p = 2,
  #                                                      method = "Wasserstein",
  #                                                      solver = "mosek",
  #                                                      search = "LBFGS",
  #                                                      penalty = "L2",
  #                                                      add.divergence = TRUE,
  #                                                      stepsize = 1e-3))
  
})

testthat::test_that("works for Wass divergence, grid/formula", {
  testthat::skip_on_cran()
  causalOT:::skip_if_no_geomloss()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  # no grid search
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = list(penalty = 10000), 
                                                       estimand = e, p = 2,
                                                       method = "Wasserstein",
                                                       solver = "mosek",
                                                       search = "pgd",
                                                       penalty = "entropy",
                                                       add.divergence = TRUE,
                                                       stepsize = 1e-2,
                                                       formula = "~.+0",
                                                       niter = 5L,
                                                       balance.constraints = 0.1))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  testthat::skip("Interactive only")
  # grid search
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       grid.search = TRUE,
                                                       estimand = e, p = 2,
                                                       method = "Wasserstein",
                                                       solver = "mosek",
                                                       search = "pgd",
                                                       penalty = "entropy",
                                                       add.divergence = TRUE,
                                                       stepsize = 1e-2,
                                                       formula = "~.+0",
                                                       balance.constraints = 0.1,
                                                       n.boot = 10))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n0,n0), 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n0,n0), 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(rep(1/n1,n1), 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(rep(1/n1,n1), 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for SBW", {
  
  set.seed(23483)
  n <- 2^8
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 3, 
                                                       estimate = e, 
                                                       method = "SBW",
                                                       solver = "osqp"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  # testthat::skip_if_not_installed("gurobi")
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimate = e, 
  #                                                          method = "SBW",
  #                                                          solver = "gurobi"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, 
                                                           estimate = e, 
                                                           method = "SBW",
                                                           solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimate = e, 
  #                                                          method = "SBW",
  #                                                          solver = "cplex"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  
  
})

testthat::test_that("works for SBW grid", {
  testthat::skip_on_cran()
  
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint <- seq(0,1, 10)
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = constraint, 
                                                       grid.search = TRUE,
                                                       grid = 10,
                                                       estimand = e, 
                                                       method = "SBW",
                                                       solver = "osqp"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  sapply(weights[1:3], function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  sapply(weights[4], function(w) testthat::expect_equal(w$args$standardized.mean.difference , c(10,10)))
  
  # testthat::skip_if_not_installed("gurobi")
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                      constraint = constraint, 
  #                                                      grid.search = TRUE,
  #                                                      grid = 10,
  #                                                      estimand = e, 
  #                                                      method = "SBW",
  #                                                      solver = "gurobi"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  # sapply(weights[1:3], function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  # sapply(weights[4], function(w) testthat::expect_equal(w$args$standardized.mean.difference , c(10,10)))
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = constraint,
                                                       grid.search = TRUE,
                                                       grid = 10,
                                                       estimand = e, 
                                                       method = "SBW",
                                                       solver = "mosek"))
  sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  sapply(weights[1:3], function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  sapply(weights[4], function(w) testthat::expect_equal(w$args$standardized.mean.difference , c(10,10)))
  
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                      constraint = constraint,
  #                                                      grid.search = TRUE,
  #                                                      grid = 10,
  #                                                      estimand = e, 
  #                                                      method = "SBW",
  #                                                      solver = "cplex"))
  # sapply(weights, function(w) testthat::expect_equal(names(w), arg.names))
  # sapply(weights[1:3], function(w) testthat::expect_lte(w$args$standardized.mean.difference , 10))
  # sapply(weights[4], function(w) testthat::expect_equal(w$args$standardized.mean.difference , c(10,10)))
  # 
  
  
})

testthat::test_that("works for Wass, sw", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 2, p = 4,
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = "osqp",
                                                       sample_weight = sample_weights))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  
  # testthat::skip_if_not_installed("gurobi")
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                      constraint = 2, 
  #                                                      estimand = e, p = 4,
  #                                                      method = "Wasserstein",
  #                                                      solver = "gurobi",
  #                                                      sample_weight = sample_weights))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 2, p = 4,
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = "mosek",
                                                       sample_weight = sample_weights))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[3]]$w0, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[4]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[3]]$w1, check.attributes = FALSE), "Mean relative difference")
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[4]]$w1, check.attributes = FALSE), "Mean relative difference")
})

testthat::test_that("works for Wass, grid/formula, sw", {
  testthat::skip_on_cran()
  causalOT:::skip_if_no_geomloss()
  set.seed(23483)
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  for (e in estimates) {
    testthat::expect_silent(weight.check[[e]] <- calc_weight(data = data, 
                                                             constraint = NULL,
                                                             grid.search = TRUE,
                                                             estimand = e, p = 4,
                                                             formula = "~.+0",
                                                             balance.constraints = 1,
                                                             method = "Wasserstein",
                                                             solver = "osqp",
                                                             wass.method = "greenkhorn",
                                                             iter = 10,
                                                             sample_weight = sample_weights, eval.method = "bootstrap"))
  }
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  for (e in estimates) {
    testthat::expect_silent(weight.check[[e]] <- calc_weight(data = data, 
                                                              constraint = NULL,
                                                              grid.search = TRUE,
                                                              estimand = e, p = 4,
                                                              formula = "~.+0",
                                                              balance.constraints = 1,
                                                              method = "Wasserstein",
                                                              solver = "mosek",
                                                              wass.method = "greenkhorn",
                                                              iter = 10,
                                                              sample_weight = sample_weights, eval.method = "bootstrap"))
  }
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
})

testthat::test_that("works for NNM, sample weight", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w0[seq(1,n0,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 2.5, p = 4,
                                                           estimand = e, 
                                                           method = "NNM",
                                                           solver = "osqp",
                                                           sample_weight = sample_weights))
  for (w in weights) testthat::expect_equal(names(w), arg.names)
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  # weights <- lapply(estimates, function(e) calc_weight(data = data, 
  #                                                          constraint = 3, 
  #                                                          estimand = e, 
  #                                                          method = "NNM",
  #                                                          solver = "cplex",
  #                                                          sample_weight = sample_weights))
  # for(w in weights) testthat::expect_equal(names(w), arg.names)
  # testthat::expect_match(all.equal(sample_weights$w0, 
  #                                  weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  # 
  # testthat::expect_match(all.equal(sample_weights$w1, 
  #                                  weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                           constraint = 3, p = 4,
                                                           estimand = e, 
                                                           method = "NNM",
                                                           solver = "mosek",
                                                           sample_weight = sample_weights))
  for(w in weights) testthat::expect_equal(names(w), arg.names)
  
  testthat::expect_match(all.equal(sample_weights$w0, 
                                   weights[[1]]$w0, check.attributes = FALSE), "Mean relative difference")
  
  testthat::expect_match(all.equal(sample_weights$w1, 
                                   weights[[2]]$w1, check.attributes = FALSE), "Mean relative difference")
  
})

testthat::test_that("works for SCM grid", {
  testthat::skip_on_cran()
  causalOT:::skip_if_no_geomloss()
  skip_if
  set.seed(23483)
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  # testthat::expect_warning(
    weight.check[[estimates[1]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                grid.search = TRUE,
                                                estimand = estimates[1], 
                                                method = "SCM",
                                                solver = "mosek",
                                                wass.method = "sinkhorn",
                                                iter = 10, eval.method = "bootstrap")
  # )
  testthat::expect_silent(
    weight.check[[estimates[2]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                grid.search = TRUE,
                                                estimand = estimates[2], 
                                                method = "SCM",
                                                solver = "mosek",
                                                wass.method = "sinkhorn",
                                                iter = 10, eval.method = "bootstrap")
  )
  testthat::expect_silent(
    weight.check[[estimates[3]]] <- calc_weight(data = data, 
                                     constraint = NULL,
                                     grid.search = TRUE,
                                     estimand = estimates[3], 
                                     method = "SCM",
                                     solver = "mosek",
                                     wass.method = "sinkhorn",
                                     iter = 10, eval.method = "bootstrap")
  )
  
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
})

testthat::test_that("works for probit vs logit", {
  
  set.seed(23483)
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  # testthat::expect_warning(
  weight.check[[estimates[1]]] <- calc_weight(data = data, 
                                              constraint = NULL,
                                              method = "Logistic",
                                              estimand = estimates[1])
  # )
  testthat::expect_silent(
    weight.check[[estimates[2]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                method = "Logistic",
                                                estimand = estimates[2]
                                                )
  )
  testthat::expect_silent(
    weight.check[[estimates[3]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                method = "Logistic",
                                                estimand = estimates[3])
  )
  
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
  
  
  weight.check <- vector("list", length(estimates))
  names(weight.check) <- estimates
  
  # testthat::expect_warning(
  weight.check[[estimates[1]]] <- calc_weight(data = data, 
                                              constraint = NULL,
                                              method = "Probit",
                                              estimand = estimates[1])
  # )
  testthat::expect_silent(
    weight.check[[estimates[2]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                method = "Probit",
                                                estimand = estimates[2]
    )
  )
  testthat::expect_silent(
    weight.check[[estimates[3]]] <- calc_weight(data = data, 
                                                constraint = NULL,
                                                method = "Probit",
                                                estimand = estimates[3])
  )
  
  for (w in weight.check) testthat::expect_equal(names(w), arg.names)
})