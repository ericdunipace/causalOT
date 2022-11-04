testthat::test_that("qp_wass_const correct", {
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
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, 
                      cost = cost)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                         c(1, rep(1/n1,n1), mean(cost^power)))
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 4:5
  
  #### get simulation functions ####
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  
  testthat::expect_error(qp <- qp_wass_const(x,z,K = constraint, p = power,cost = cost))
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 4
  
  #### get simulation functions ####
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
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
  
  qp <- qp_wass_const(x,z,K = constraint, p = power,  cost = cost,
                      sample_weight = sample_weights)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1), mean(cost^power)))
  
  
})

testthat::test_that("qp_wass_const correct, margins", {
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
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  marg.cost <- lapply(1:p, function(i) causalOT:::cost_calc_lp(x[z==0,i,drop = FALSE], x[z==1,i,drop = FALSE]))
  cost.list <- c(marg.cost, list(cost))
  
  add.margins <- TRUE
  add.joint <- TRUE
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, cost = cost.list,
                      add.margins = add.margins, add.joint = add.joint)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1), mean(cost^power),
                                sapply(marg.cost, function(m) mean(m^power))))
  
  
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
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  constraint = 1
  power <- 1:2
  
  #### get simulation functions ####
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  
  testthat::expect_error(qp <- qp_wass_const(x,z,K = constraint, p = power, cost = cost,
                                             add.margins = TRUE, add.joint = TRUE))
  
})

testthat::test_that("qp_wass_const correct, margins with sw", {
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
  data <- causalOT:::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  x <- data$get_x()
  z <- data$get_z()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  cost <- causalOT:::cost_calc_lp(x[z==0,], x[z==1,])
  marg.cost <- lapply(1:p, function(i) causalOT:::cost_calc_lp(x[z==0,i,drop = FALSE], x[z==1,i,drop = FALSE]))
  cost.list <- c(marg.cost, list(cost))
  
  add.margins <- TRUE
  add.joint <- TRUE
  
  w0 <- renormalize(runif(n0))
  w1 <- renormalize(runif(n1))
  w1[seq(1,n1,2)] <- 0
  w1[seq(1,n1,2)] <- 0
  sample_weights <- list(w0 = renormalize(w0),
                         w1 = renormalize(w1))
  
  qp <- qp_wass_const(x,z,K = constraint, p = power, cost = cost.list,
                      add.margins = add.margins, add.joint = add.joint,
                      sample_weight = sample_weights)
  testthat::expect_equivalent((qp$LC$A %*% rep(1/(n0*n1), n0*n1))@x,
                              c(1, rep(1/n1,n1),  mean(cost^power),
                                sapply(marg.cost, function(m) mean(m^power))
                               ))
  
  
})

