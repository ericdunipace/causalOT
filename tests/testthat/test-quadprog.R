testthat::test_that("quadprog.DataSim uses right wasserstein qp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Constrained Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method, cost = cost)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATT", method = method, cost = cost)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATE", method = method)
  
  cost_n0 <- causalOT::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(mean(cost_n1^power), 1, rep(1/n,n)))
  
  
})

testthat::test_that("quadprog.data.frame uses right wasserstein qp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Constrained Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  df <- data.frame(y = rep(0,n), z = z, x)
  
  bal.cov <- paste0("X", 1:p)
  
  cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATC", method = method, 
                            cost = cost, balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y")
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                         estimand = "ATT", 
                         method = method, cost = cost,
                         balance.covariates = bal.cov,
                         treatment.indicator = "z",
                         outcome = "y")
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                         estimand = "ATE", method = method,
                         balance.covariates = bal.cov,
                         treatment.indicator = "z",
                         outcome = "y")
  
  cost_n0 <- causalOT::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(mean(cost_n1^power), 1, rep(1/n,n)))
  
  
})

testthat::test_that("quadprog.DataSim uses right wasserstein qp with sample_weights", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Constrained Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w1[seq(1,n1,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.DataSim(data, constraint = constraint, 
                         estimand = "ATC", method = method, cost = cost,
                         sample_weight = sample_weights)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATT", method = method, cost = cost)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATE", 
                         method = method,
                         sample_weight = sample_weights)
  
  cost_n0 <- causalOT::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(mean(cost_n1^power), 1, rep(1/n,n)))
  
  
})

testthat::test_that("quadprog.data.frame uses right wasserstein qp with sample weights", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Constrained Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w1[seq(1,n1,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  df <- data.frame(y = rep(0,n), z = z, x)
  
  bal.cov <- paste0("X", 1:p)
  
  cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATC", method = method, 
                            cost = cost, balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATT", 
                            method = method, cost = cost,
                            balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATE", method = method,
                            balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  
  cost_n0 <- causalOT::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(mean(cost_n1^power), 1, rep(1/n,n)))
  
  
})

