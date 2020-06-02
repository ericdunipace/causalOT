testthat::test_that("check RKHS kernel works", {
  set.seed(23483)
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE","ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  #### run sims  ####
  # debugonce(calc_weight_RKHS)
  # debugonce(quadprog.DataSim)
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimate = "ATE",  opt.hyperparam = FALSE,
                   solver = c("gurobi"), theta = c(1,1), gamma = c(1,1), 
                   p = power[1], metric = "mahalanobis"))
  
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimate = "ATT",  opt.hyperparam = FALSE,
                                                  solver = c("gurobi"), theta = c(1,1), gamma = c(1,1), 
                                                  p = power[1], metric = "mahalanobis"))
  
  testthat::expect_silent(out <- calc_weight_RKHS(data, estimate = "ATC",  opt.hyperparam = FALSE,
                                                  solver = c("gurobi"), theta = c(1,1), gamma = c(1,1), 
                                                  p = power[1], metric = "mahalanobis"))
})

testthat::test_that("check RKHS kernel works with hyper param opt", {
  set.seed(23483)
  n <- 2^9
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE","ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  #### run sims  ####
  # debugonce(calc_weight_RKHS)
  # debugonce(quadprog.DataSim)
  for(e in estimates) {
    testthat::expect_silent(out <- calc_weight_RKHS(data, estimate = e, opt.hyperparam = TRUE,
                                                  solver = c("gurobi"), theta = c(1,1), gamma = c(1,1), 
                                                  p = power[1], metric = "mahalanobis", iter = 10))
  }
  
})
