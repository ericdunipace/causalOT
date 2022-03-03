testthat::skip("RKHS method deprecated")

testthat::test_that("RKHS grid works", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(9870)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  ground_power <- 2
  solver <- "mosek"
  estimands <- c("ATT", "ATC", "CATE","feasible")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  # testthat::expect_silent(cplex.check  <- RKHS_grid_search(data = data, grid = NULL, estimand = "ATE", n.boot = 10, opt.hyperparam = FALSE, solver = "cplex"))
  # testthat::expect_silent(gurobi.check <- RKHS_grid_search(data = data, grid = NULL, estimand = "ATE", n.boot = 10, opt.hyperparam = FALSE, solver = "gurobi"))
  # testthat::expect_error(
  #   testthat::expect_warning(
      mosek.check  <- RKHS_grid_search(data = data, grid = NULL, 
                                       estimand = "ATE", n.boot = 10, 
                                       opt.hyperparam = FALSE, 
                                       solver = "mosek")
  #     )
  # )
  
})

testthat::test_that("RKHS grid works, opt.hyperparam", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(9870)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  ground_power <- 2
  solver <- "mosek"
  estimands <- c("ATT", "ATC","feasible")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  # testthat::expect_silent(cplex.check <- RKHS_grid_search(data = data, grid = NULL, estimand = "ATE", n.boot = 10, opt.hyperparam = TRUE, solver = "cplex", iter = 10))
  # testthat::expect_silent(gurobi.check <- RKHS_grid_search(data = data, grid = NULL, estimand = "ATE", n.boot = 10, opt.hyperparam = TRUE, solver = "gurobi", iter = 10))
  testthat::expect_condition(mosek.check <- RKHS_grid_search(data = data, grid = NULL, 
                                                             estimand = "ATE", n.boot = 10, 
                                                             opt.hyperparam = TRUE, solver = "mosek", iter = 10),
                             regexp = c("Algorithm|All")) # warning and error
  
})

testthat::test_that("RKHS grid gives expected value", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  set.seed(9870)
  library(causalOT)
  
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  solver <- "mosek"
  
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  weight <- RKHS_grid_search(data, grid = seq(0,100, length.out = 11), 
                             method = "RKHS.dose",
                            n.boot = 10,
                            solver = solver,
                            theta = c(1,1),
                            gamma = c(1,1),
                            p = 1,
                            opt.hyperparam = FALSE,
                            dist = distance,
                            estimand = "ATE")
  test <- do.call("calc_weight_RKHS", 
                  list(data = data, 
                       method = "RKHS.dose",
                       estimand = "ATE",
                       solver = solver,
                       opt.hyperparam = FALSE,
                       lambda = weight$lambda,
                       theta = c(1,1),
                       gamma = c(1,1),
                       p = 1,
                       dist = distance))
  test.fun <- function(n1,n2){testthat::expect_equal(n1, n2)}
  mapply(test.fun, n1 = weight[1:3], n2 = test[1:3])
})

