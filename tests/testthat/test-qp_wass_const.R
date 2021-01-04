testthat::test_that("qp_wass_const correct", {
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
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATC", cost = cost)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                         c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATT", cost = cost)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                         c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATE")
  
  cost_n0 <- causalOT::cost_calc_lp(x[z==0,], x)
  cost_n1 <- causalOT::cost_calc_lp(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                         c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                         c(mean(cost_n1^power), 1, rep(1/n,n)))
  
})

testthat::test_that("qp_wass_const gives error if power has length > 1 ", {
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
  power <- 1:2
  
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
  
  testthat::expect_error(qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATC", cost = cost))
  
})

testthat::test_that("qp_wass_const correct", {
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
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATC", cost = cost,
                      sample_weight = sample_weights)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n0,n0)))
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATT", cost = cost,
                      sample_weight = sample_weights)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(mean(cost^power), 1, rep(1/n1,n1)))
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, estimand = "ATE",
                      sample_weight = sample_weights)
  
  cost_n0 <- causalOT::cost_calc_lp(x[z==0,], x)
  cost_n1 <- causalOT::cost_calc_lp(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(mean(cost_n0^power), 1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(mean(cost_n1^power), 1, rep(1/n,n)))
  
})
