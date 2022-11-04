testthat::test_that("qp_wass gives error if power has length > 1 ", {
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
  constraint = 1
  power <- 4:5
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  # cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  
  testthat::expect_error(qp <- qp_wass(x,z,K = constraint, p = power, estimand = "ATC"))
  
})

testthat::test_that("qp_wass has correct dimensions, no marginals", {
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
  constraint = 1
  power <- 4
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  # cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  #debugonce(qp_wass)
  qp <- qp_wass(x,z,K = constraint, p = power)
  
  testthat::expect_equal(prod(dim(qp$obj$Q)), (n0*n1)^2, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[1], 1 + n1, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[2], n0*n1, check.attributes=FALSE)
  
  
})

testthat::test_that("qp_wass has correct dimensions, marginals", {
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
  constraint = 1
  power <- 4
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  # cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  #debugonce(qp_wass)
  qp <- qp_wass(x,z,K = constraint, p = power,
                add.margins = TRUE)
  
  testthat::expect_equal(prod(dim(qp$obj$Q)), (n0*n1)^2, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[1], 1 + n1 + p, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[2], n0*n1, check.attributes=FALSE)
  
  
})

testthat::test_that("qp_wass works with sample weights", {
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
  constraint = 1
  power <- 4
  
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
  # cost <- causalOT::cost_calc_lp(x[z==0,], x[z==1,])
  #debugonce(qp_wass)
  qp <- qp_wass(x,z,K = constraint, p = power, 
                sample_weight = sample_weights, add.margins = TRUE)
  
  testthat::expect_equal(prod(dim(qp$obj$Q)), (n0*n1)^2, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[1], 1 + n1 + p, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[2], n0*n1, check.attributes=FALSE)
  
  
  qp <- qp_wass(x,z,K = constraint, p = power,
                sample_weight = sample_weights)
  
  testthat::expect_equal(prod(dim(qp$obj$Q)), (n0*n1)^2, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[1], 1 + n1, check.attributes=FALSE)
  testthat::expect_equal( dim(qp$LC$A)[2], n0*n1, check.attributes=FALSE)
  
  
})
