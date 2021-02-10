testthat::test_that("grid search actually works, cwass", {
  testthat::skip_on_cran()
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
  
  testthat::expect_lte(wsel$args$constraint - 2.191191, 1e-3)
  
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = seq(minwass, wass_full, length.out = 10),
                                                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                   metric = metric, p = power, solver = "mosek",
                                                   wass.method = "greenkhorn", wass.iter = 1000))
  
  testthat::expect_lte(wsel$args$constraint - 2.327817, 1e-3)
  
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = seq(minwass, wass_full, length.out = 10),
                                                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                   metric = metric, p = power, solver = "mosek",
                                                   wass.method = "greenkhorn", wass.iter = 1000,
                                                   unbiased = TRUE))
  
  testthat::expect_lte(wsel$args$constraint - 2.464442, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0)
    )
  testthat::expect_lte(wsel2$args$constraint$joint - 2.31591, 1e-3)
  
  estimand <- "ATC"
  testthat::expect_silent(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  
})

testthat::test_that("grid search actually works, cwass mahal", {
  testthat::skip_on_cran()
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
  
  testthat::expect_lte(wsel$args$constraint - 2.07741, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel2 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  testthat::expect_lte(wsel2$args$constraint$joint - 2.238248, 1e-3)
  # testthat::expect_equal(wsel, wsel2)
})

testthat::test_that("grid search actually works, cwass sdlp", {
  testthat::skip_on_cran(); set.seed(9867)
  
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
  method <- "Constrained Wasserstein"
  add.joint <- TRUE
  add.marginal <- FALSE
  
  
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
  testthat::expect_silent(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
  testthat::expect_lte(wsel$args$constraint$joint[1] - c(1.972782), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = add.joint, add.margins = add.marginal)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
})

testthat::test_that("grid search actually works, cwass lp marg", {
  testthat::skip_on_cran(); 
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
  add.margins <- TRUE
  add.joint <- TRUE
  
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
                                                   wass.method = "networkflow", wass.iter = 0,
                                                   add.margins = FALSE, add.joint = add.joint))
  
  testthat::expect_lte(wsel$args$constraint - 2.327817, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              add.margins = add.margins, add.joint = add.joint)
  )
  testthat::expect_equivalent(wsel2$args$constraint, list(joint = c(1.7188448, 1.2600821, 1.1800639, 2.0879501, 1.9524347, 0.5670994 ),
                                                          margins = c(2.31591)), 
                       tol = 1e-3)
  
  # estimand <- "ATC"
  # testthat::expect_warning(wsel3 <- wass_grid_search(data, grid = NULL,
  #                                                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
  #                                                   metric = metric, p = power, solver = "mosek",
  #                                                   wass.method = "networkflow", wass.iter = 0,
  #                                                   add.margins = add.margins, add.joint = add.joint))
  # 
  # estimand <- "ATE"
  # testthat::expect_warning(wsel4 <- wass_grid_search(data, grid = NULL,
  #                                                   estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
  #                                                   metric = metric, p = power, solver = "mosek",
  #                                                   wass.method = "networkflow", wass.iter = 0,
  #                                                   add.margins = add.margins, add.joint = add.joint))
  # 
  
})

testthat::test_that("grid search actually works, cwass mahal marg", {
  testthat::skip_on_cran(); 
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
  add.margins <- TRUE
  add.joint <- TRUE
  
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
                                                   wass.method = "networkflow", wass.iter = 0,
                                                   add.margins = add.margins, add.joint = add.joint))
  
  testthat::expect_lte(wsel$args$constraint - 2.181136, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel2 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.margins = add.margins, add.joint = add.joint))
  testthat::expect_equivalent(wsel2$args$constraint,
                              list(margins = c(0.8926885, 0.8926885, 0.8926885, 0.8926885, 
                                0.8926885, 0.8926885),
                                joint = 2.370255),
                              1e-2)
  # testthat::expect_equal(wsel, wsel2)
})

testthat::test_that("grid search actually works, cwass sdlp marg", {
  testthat::skip_on_cran(); 
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
  method <- "Constrained Wasserstein"
  add.joint <- TRUE
  add.marginal <- TRUE
  
  
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
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, 
                              list( margins = rep(1.442966, 6), 
                                    joint = 1.972782), tol = 1e-2)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                   estimand = estimand, n.boot = 10, method = method,
                                                   metric = metric, p = power, solver = "mosek",
                                                   wass.method = "networkflow", wass.iter = 0,
                                                   add.joint = add.joint, add.margins = add.marginal)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
})

testthat::test_that("grid search actually works,  wass sdlp", {
  testthat::skip_on_cran(); 
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
  method <- "Wasserstein"
  add.margins <- FALSE
  add.joint <- TRUE
  
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
  testthat::expect_silent(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins
                             )
    
    )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(112.4464), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = method,
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = add.joint, add.margins = add.margins)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = method,
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = add.joint, add.margins = add.margins)
  
  )
  
})

testthat::test_that("grid search actually works,  wass mahal", {
  testthat::skip_on_cran(); 
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
  method <- "Wasserstein"
  
  
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
  testthat::expect_silent(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = FALSE)
    
  )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(111.1607), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = FALSE)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = FALSE)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass sdlp", {
  testthat::skip_on_cran(); 
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
  method <- "Wasserstein"
  add.margins <- TRUE
  add.joint <- TRUE
  
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
  testthat::expect_silent(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins
    )
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint,
                              list(margins = rep(2.804821, 6), 
                                   penalty = 112.4464), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = add.joint, add.margins = add.margins)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal", {
  testthat::skip_on_cran(); 
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
  method <- "Wasserstein"
  add.margins <- TRUE
  add.joint <- TRUE
  
  
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
  testthat::expect_silent(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.margins)
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, list(margins = rep(2.376389, 6), penalty = 111.1607), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.margins)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.margins)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass sdlp bal", {
  testthat::skip_on_cran(); 
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
  method <- "Wasserstein"
  add.joint = TRUE
  add.marginal <- TRUE
  formula <- "~ . + 0"
  balance.constraints <- 0.2
  
  
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
  testthat::expect_silent(
    # debugonce(wass_grid_search)
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal,
                             formula = formula, balance.constraints = balance.constraints)
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, 
                              list(margins = rep(1.654497, 6),
                                   112.4464), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.marginal,
                                                   formula = formula, balance.constraints = balance.constraints)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal bal", {
  testthat::skip_on_cran(); 
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
  formula <- "~. + 0"
  balance.constraints <- 0.1
  add.joint <- TRUE
  add.marginal <- FALSE
  method = "Wasserstein"
  
  
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
  testthat::expect_silent(
  # debugonce(wass_grid_search)
  wsel <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = method,
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 100,
                           add.joint = add.joint, add.margins = add.marginal,
                           formula = formula, balance.constraints = balance.constraints)

  )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(111.1607  ), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.marginal,
                                                    formula = formula, balance.constraints = balance.constraints)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal,
                             formula = formula, balance.constraints = balance.constraints)
    
  )
  
})

testthat::test_that("grid search joint.map, cwass", {
  testthat::skip_on_cran()
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
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              joint.mapping = TRUE, penalty = "L2")
  )
  testthat::expect_lte(wsel2$args$constraint$joint - 2.31591, 1e-3)
  
  estimand <- "ATC"
  testthat::expect_silent(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  
})

testthat::test_that("grid search joint.map, wass", {
  testthat::skip_on_cran()
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
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "Wasserstein",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              joint.mapping = TRUE, penalty = "L2")
  )
  testthat::expect_equivalent(wsel2$args$constraint, list(penalty = 11.65043,
                                                          joint = 0.5000005), 1e-3)
  
  estimand <- "ATC"
  testthat::expect_silent(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    joint.mapping = TRUE,
                                                    wass.method = "networkflow", wass.iter = 0))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    joint.mapping = TRUE,
                                                    wass.method = "networkflow", wass.iter = 0))
  
  
})

testthat::test_that("grid search scm", {
  testthat::skip_on_cran()
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
  method <- "SCM"
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "SCM",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              joint.mapping = TRUE, penalty = "L2")
  )
  testthat::expect_lte(wsel2$args$constraint$penalty - 135.4599, 1e-3)
  
  estimand <- "ATC"
  testthat::expect_silent(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "SCM",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "SCM",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  
  
})

