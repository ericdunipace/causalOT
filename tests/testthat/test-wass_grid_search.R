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
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                           metric = metric, p = power, solver = "mosek",
                           wass.method = "networkflow", wass.iter = 0)
    )
  testthat::expect_lte(wsel2$args$constraint - 2.054565, 1e-3)
  
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
  
  testthat::expect_lte(wsel$args$constraint - 2.057709, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel2 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0))
  testthat::expect_lte(wsel2$args$constraint - 2.185474, 1e-3)
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
  
  testthat::expect_lte(wsel$args$constraint[1] - c(1.910), 1e-3)
  
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
  testthat::skip_on_cran(); set.seed(9867)
  
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
                                                   add.margins = add.margins, add.joint = add.joint))
  
  testthat::expect_lte(wsel$args$constraint - 2.054565, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              add.margins = add.margins, add.joint = add.joint)
  )
  testthat::expect_equivalent(wsel2$args$constraint, c(.954913793028802, 0.700045646077501, 0.655591044438986, 1.15997231098634, 1.08468596963478, 0.315055236968164, 2.45367921448841), 
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
  testthat::skip_on_cran(); set.seed(9867)
  
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
  
  testthat::expect_lte(wsel$args$constraint - 2.057709, 1e-3)
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel2 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Constrained Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.margins = add.margins, add.joint = add.joint))
  testthat::expect_equivalent(wsel2$args$constraint,
                              c(0.9520001, 0.9520001, 0.9520001, 0.9520001, 
                                0.9520001, 0.9520001, 1.97),
                              1e-2)
  # testthat::expect_equal(wsel, wsel2)
})

testthat::test_that("grid search actually works, cwass sdlp marg", {
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
                              c(rep(0.8016476, 6), 1.86), tol = 1e-2)
  
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
  
  testthat::expect_lte(wsel$args$constraint[1] - c(421.8279), 1e-3)
  
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
  testthat::skip_on_cran(); set.seed(9867)
  
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
  
  testthat::expect_lte(wsel$args$constraint[1] - c(421.4782), 1e-3)
  
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
  testthat::expect_warning(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins
    )
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint,
                              c(rep(0.802, 6), 421.697+06), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = add.joint, add.margins = add.margins)
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal", {
  testthat::skip_on_cran(); set.seed(9867)
  
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
  testthat::expect_warning(
    
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.margins)
    
  )
  
  testthat::expect_lte(wsel$args$constraint[1] - c( 0.9752433), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
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
  testthat::expect_warning(
    # debugonce(wass_grid_search)
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal,
                             formula = formula, balance.constraints = balance.constraints)
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, c(rep(0.8129247, 6),
                                                   421.6965034), tol = 1e-3)
  
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
                             add.joint = add.joint, add.margins = add.marginal)
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal bal", {
  testthat::skip_on_cran(); set.seed(9867)
  
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
  
  testthat::expect_lte(wsel$args$constraint[1] - c(421.697  ), 1e-3)
  
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
