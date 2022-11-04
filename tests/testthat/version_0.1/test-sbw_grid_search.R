testthat::test_that("grid search function works, dataSim", {
  set.seed(9870)
  testthat::skip_on_cran()
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  library(causalOT)
  
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  ground_power <- 2
  solver <- "mosek"
  estimates <- c("ATT", "ATC","feasible","ATE","cATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  data$gen_data()
  
  test.fun <- function(g1,g2) {
    testthat::expect_equal(g1,g2)}
  for (e in estimates) {
    weight <- sbw_grid_search(data, grid = seq(0,0.5, length.out = 10), 
                  estimand = e,
                              n.boot = 100,
                  solver = solver)
    weight$args$standardized.mean.difference <- NULL
    mapply(test.fun, g1 = weight, g2 = do.call("calc_weight_bal", 
                                           list(data = data, 
                                                constraint = weight$args$constraint,  
                                                estimand = e, 
                                               method = "SBW", solve = solver)))
  }
  
})

testthat::test_that("grid search deletes extra args", {
  testthat::skip_on_cran()
  set.seed(9870)
  library(causalOT)
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(1,2)
  ground_power <- 2
  solver <- "mosek"
  estimates <- c("ATT", "ATC","feasible","ATE","cATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  
  data$gen_data()
  
  test.fun <- function(g1,g2) {testthat::expect_equal(g1,g2)}
  for(e in estimates) {
    testthat::expect_silent(weight <- sbw_grid_search(data, grid = seq(0,0.5, length.out = 10), 
                              constraint = 50,
                              estimand = e,
                              n.boot = 100,
                              solver = solver))
    weight$args$standardized.mean.difference <- NULL
    mapply(test.fun, g1 = weight, g2 = do.call("calc_weight_bal", 
                                                    list(data = data, 
                                                         constraint = weight$args$constraint,  
                                                         estimand = e, 
                                                         method = "SBW", solve = solver)))
  }
  
})
