testthat::test_that("grid search actually works, cwass", {
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
  wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
                             b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
                             cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
                             p = power
  )
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = seq(minwass, wass_full, length.out = 10),
                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                   metric = metric, p = power, solver = "mosek",
                   wass.method = "networkflow", wass.iter = 0))
  
  testthat::expect_lte(wsel$args$constraint - 1.917939, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel2 <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0))
  testthat::expect_lte(wsel2$args$constraint - 2.327817, 1e-3)
})

testthat::test_that("grid search actually works, cwass mahal", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "mahalanobis"
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
  wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
                             b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
                             cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
                             p = power
  )
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = seq(minwass, wass_full, length.out = 10),
                                                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                   metric = metric, p = power, solver = "mosek",
                                                   wass.method = "networkflow", wass.iter = 0))
  
  testthat::expect_lte(wsel$args$constraint - 2.054565, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel2 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  testthat::expect_lte(wsel2$args$constraint - 2.185474, 1e-3)
  # testthat::expect_equal(wsel, wsel2)
})

testthat::test_that("grid search actually works, marg wass sdlp", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "sdLp"
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  # NNM <- calc_weight(data, estimand = estimand, method = "NNM", tranport.matrix = TRUE, metric = metric)
  # minwass <- wasserstein_p(a = NNM, cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), p = power)
  # wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
  #                            b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
  #                            cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
  #                            p = power
  # )
  # debugonce(wass_grid_search)
  estimand <- "ATT"
  testthat::expect_warning(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = FALSE)
    
    )
  
  testthat::expect_lte(wsel$args$constraint[1] - c(0.8286069), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = FALSE)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = FALSE)
  
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "mahalanobis"
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  # NNM <- calc_weight(data, estimand = estimand, method = "NNM", tranport.matrix = TRUE, metric = metric)
  # minwass <- wasserstein_p(a = NNM, cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), p = power)
  # wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
  #                            b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
  #                            cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
  #                            p = power
  # )
  # debugonce(wass_grid_search)
  estimand <- "ATT"
  testthat::expect_warning(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = FALSE)
    
  )
  
  testthat::expect_lte(wsel$args$constraint[1] - c(0.845143), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = FALSE)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = FALSE)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass sdlp joint", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "sdLp"
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  # NNM <- calc_weight(data, estimand = estimand, method = "NNM", tranport.matrix = TRUE, metric = metric)
  # minwass <- wasserstein_p(a = NNM, cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), p = power)
  # wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
  #                            b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
  #                            cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
  #                            p = power
  # )
  # debugonce(wass_grid_search)
  estimand <- "ATT"
  testthat::expect_warning(
    # debugonce(wass_grid_search)
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE)
    
  )
  
  testthat::expect_lte(wsel$args$constraint[1] - c(1.442966 ), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = TRUE)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal joint", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "mahalanobis"
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  # NNM <- calc_weight(data, estimand = estimand, method = "NNM", tranport.matrix = TRUE, metric = metric)
  # minwass <- wasserstein_p(a = NNM, cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), p = power)
  # wass_full <- wasserstein_p(a = rep(1/nrow(data$get_x0()),nrow(data$get_x0())), 
  #                            b = rep(1/nrow(data$get_x1()), nrow(data$get_x1())),
  #                            cost = cost_calc_lp(data$get_x0(), data$get_x1(), ground_p = power), 
  #                            p = power
  # )
  # debugonce(wass_grid_search)
  estimand <- "ATT"
  testthat::expect_warning(
  # debugonce(wass_grid_search)
  wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = TRUE)

  )
  
  testthat::expect_lte(wsel$args$constraint[1] - c(.952 ), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = TRUE)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = "Wasserstein",
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE)
    
  )
  
})