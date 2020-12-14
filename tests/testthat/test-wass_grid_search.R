testthat::test_that("grid search actually works", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "Lp"
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  estimand <- "ATT"
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  NNM <- calc_weight(data, estimand = estimand, method = "NNM", tranport.matrix = TRUE, metric = metric)
  minwass <- wasserstein_p(a = NNM, cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), p = power)
  wass_full <- wasserstein_p(a= rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
                             b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
                             cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
                             p = power
  )
  # debugonce(wass_grid_search)
  wsel <- wass_grid_search(data, grid = seq(minwass, wass_full, length.out = 10),
                   estimand = estimand, n.boot = 100, method = "Constrained Wasserstein",
                   metric = metric, p = power, solver = "mosek",
                   wass.method = "networkflow", wass.iter = 0)
})
