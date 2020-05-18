testthat::test_that("qp_wass gives error if power has length > 1 ", {
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
  
  testthat::expect_error(qp <- qp_wass(x,z,K = constraint, p = power, estimand = "ATC", cost = cost))
  
})

