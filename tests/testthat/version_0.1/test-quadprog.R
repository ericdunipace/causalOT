testthat::test_that("quadprog.DataSim uses right wasserstein qp", {
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
  power <- 2
  method <- "Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method,cost = cost)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATT", method = method, cost = cost)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATE", method = method)
  
  cost_n0 <- causalOT:::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT:::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(1, rep(1/n,n)))
  
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Wasserstein"
  
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
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATC", method = method, 
                            cost = cost, balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y")
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                         estimand = "ATT", 
                         method = method, cost = cost,
                         balance.covariates = bal.cov,
                         treatment.indicator = "z",
                         outcome = "y")
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                         estimand = "ATE", method = method,
                         balance.covariates = bal.cov,
                         treatment.indicator = "z",
                         outcome = "y")
  
  cost_n0 <- causalOT:::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT:::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(1, rep(1/n,n)))
  
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Wasserstein"
  
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
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.DataSim(data, constraint = constraint, 
                         estimand = "ATC", method = method, cost = cost,
                         sample_weight = sample_weights)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATT", method = method, cost = cost)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATE", 
                         method = method,
                         sample_weight = sample_weights)
  
  cost_n0 <- causalOT:::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT:::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(1, rep(1/n,n)))
  
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 2
  method <- "Wasserstein"
  
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
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATC", method = method, 
                            cost = cost, balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATT", 
                            method = method, cost = cost,
                            balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1)))
  
  qp <- quadprog.data.frame(df, constraint = constraint, 
                            estimand = "ATE", method = method,
                            balance.covariates = bal.cov,
                            treatment.indicator = "z",
                            outcome = "y",
                            sample_weight = sample_weights)
  
  cost_n0 <- causalOT:::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT:::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n), n0*n))@x,
                              c(1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A %*% rep(1/(n*n1), n*n1))@x,
                              c(1, rep(1/n,n)))
  
  
})

testthat::test_that("quadprog handles all square terms", {
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
  power <- 2
  method <- "SBW"
  
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
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  # debugonce(quadprog.DataSim)
  
  mm <- model.matrix(as.formula("~ 0 + . + I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2) + I(X5^2) + I(X6^2)"),
                     data = data.frame(data$get_x0()))
  
  qp <- quadprog.DataSim(data, constraint = constraint, 
                         estimand = "ATT", method = "SBW",
                         sample_weight = sample_weights,
                         formula = "~. + I(.^2)")
  
  testthat::expect_equivalent(as.matrix(qp[[1]]$LC$A[2:13,]), t(mm))
  
  
  df <- data.frame(data$get_x(), z = data$get_z())
  sqp <- quadprog.data.frame(data = df, constraint = constraint, 
                         estimand = "ATE", method = method, cost = cost,
                         sample_weight = sample_weights,
                         formula = "~ . + I(.^2)",
                         balance.covariates = colnames(data$get_x()),
                         treatment.indicator = "z")
  
  testthat::expect_equivalent(as.matrix(sqp[[1]]$LC$A[2:13,]), t(mm))
  
  
})

testthat::test_that("quadprog.DataSim does mapping", {
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
  power <- 2
  method <- "Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method,cost = cost,
                         joint.mapping = TRUE)
  
  testthat::expect_equivalent((qp[[1]]$LC$A[,1:(n0*n1)] %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATT", method = method, 
                         cost = cost, joint.mapping = TRUE)
  testthat::expect_equivalent((qp[[1]]$LC$A[,1:(n0*n1)] %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1)))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATE", method = method
                         , joint.mapping = TRUE)
  
  cost_n0 <- causalOT:::cost_mahalanobis(x[z==0,], x)
  cost_n1 <- causalOT:::cost_mahalanobis(x[z==1,], x)
  testthat::expect_equivalent((qp[[1]]$LC$A[,1:(n0*n)] %*% rep(1/(n0*n), n0*n))@x,
                              c(1, rep(1/n,n)))
  testthat::expect_equivalent((qp[[2]]$LC$A[,1:(n*n1)] %*% rep(1/(n*n1), n*n1))@x,
                              c(1, rep(1/n,n)))
  
  
})

testthat::test_that("quadprog.DataSim various penalties", {
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
  constraint = 10
  power <- 2
  method <- "Wasserstein"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  cost <- causalOT:::cost_calc_lp(x[z == 0,], x[z == 1,])
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method,cost = cost,
                         penalty = "L2")
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)
                                ))
  
  testthat::expect_equivalent(qp[[1]]$obj$Q, Matrix::Diagonal(n0*n1,1))
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method,cost = cost,
                         penalty = "entropy")
  
  testthat::expect_equivalent((qp[[1]]$LC$A %*% c(rep(1/(n0*n1), n0*n1), rep(0, n0*n1)))@x,
                              c(1, rep(1/n0,n0)
                                ))
  testthat::expect_equivalent(qp[[1]]$obj$Q, NULL)
  testthat::expect_equivalent(length(qp[[1]]$obj$L), 2 * n0 * n1)
  
  
  qp <- quadprog.DataSim(data, constraint = constraint, estimand = "ATC", method = method,cost = cost,
                         penalty = "variance")
  
  testthat::expect_equivalent((qp[[1]]$LC$A[,1:(n0*n1)] %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n0,n0)
                                ))
  testthat::expect_equivalent(qp[[1]]$obj$Q,  NULL)
  testthat::expect_equivalent(length(qp[[1]]$obj$L[1:(n0*n1)]), n0 * n1)
  
  
})
