testthat::skip("RKHS method deprecated")

testthat::test_that("check RKHS kernel works", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(23483)
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimands <- c("ATT", "ATC", "cATE","ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  #### run sims  ####
  # debugonce(calc_weight_RKHS)
  # debugonce(quadprog.DataSim)
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimand = "ATE",  opt.hyperparam = FALSE,
                   solver = c("mosek"), theta = c(1,1), gamma = c(1,1), 
                   p = power[1], metric = "mahalanobis"))
  
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimand = "ATT",  opt.hyperparam = FALSE,
                                                  solver = c("mosek"), theta = c(1,1), gamma = c(1,1), 
                                                  p = power[1], metric = "mahalanobis"))
  
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimand = "ATC",  opt.hyperparam = FALSE,
                                                  solver = c("mosek"), theta = c(1,1), gamma = c(1,1), 
                                                  p = power[1], metric = "mahalanobis"))
})

testthat::test_that("check RKHS kernel works with hyper param opt", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="Rmosek")
  testthat::skip_on_ci()
  set.seed(23483)
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  metric <- c("Lp")
  power <- c(1,2)
  solver <- "mosek"
  estimands <- c("ATT", "ATC", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  #### run sims  ####
  # debugonce(calc_weight_RKHS)
  # debugonce(quadprog.DataSim)
  # for(e in estimands) {
  #   testthat::expect_silent(out <- calc_weight_RKHS(data, estimand = e, opt.hyperparam = TRUE,
  #                                                 solver = c("gurobi"), theta = c(1,1), gamma = c(1,1), 
  #                                                 p = power[1], metric = "mahalanobis", iter = 10))
  # }
  
})

testthat::test_that("check RKHS kernel optimal for hainmueller", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "B"
  metric <- c("mahalanobis")
  kernel <- "polynomial"
  power <- c(1,2)
  solver <- "mosek"
  estimands <- c("ATT", "ATC", "cATE","ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  #### run sims  ####
  # debugonce(calc_weight_RKHS)
  # debugonce(RKHS_param_opt)
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimand = "ATE",  opt.hyperparam = TRUE,
                                                  p = 2:3, metric = metric, kernel = kernel,
                                                  opt.method = "stan", verbose = FALSE,
                                                  algorithm = "LBFGS"))
  
  est <- estimate_effect(data = data, weights = out, estimand = "ATE")
  
  # testthat::expect_lte( abs(est$estimate - 0), 1)
  
})