testthat::test_that("grid search function works, dataSim", {
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
  solver <- "gurobi"
  estimates <- c("ATT", "ATC","feasible")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  data$gen_data()
  for(e in estimates) {
    weight <- sbw_grid_search(data, grid = seq(0,0.5, length.out = 10), 
                  estimate = e,
                              n.boot = 100,
                  solver = solver)
    testthat::expect_equal(weight[1:3], do.call("calc_weight_bal", 
                                           list(data = data, 
                                                constraint = weight$standardized.mean.difference,  
                                                estimate = e, 
                                                                   method = "SBW", solve = solver)))
  }
  
})
