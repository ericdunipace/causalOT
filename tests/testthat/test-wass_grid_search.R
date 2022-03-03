
testthat::test_that("grid search actually works,  wass sdlp", {
  testthat::skip_on_cran(); 
  set.seed(9867)
  # testthat::skip_if_not_installed("gurobi")
  #  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  # causalOT:::skip_if_no_geomloss()
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  metric <- "sdLp"
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "osqp"
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
    
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = solver,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins,
                             eval.method = "bootstrap"
                             )
    
    )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(20238.58), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = method,
                           metric = metric, p = power, solver = solver,
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = add.joint, add.margins = add.margins,
                           eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                           estimand = estimand, n.boot = 10, method = method,
                           metric = metric, p = power, solver = solver,
                           wass.method = "networkflow", wass.iter = 0,
                           add.joint = add.joint, add.margins = add.margins,
                           eval.method = "bootstrap")
  
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
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
    
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "osqp",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = FALSE,
                             eval.method = "bootstrap")
    
  )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(2605.169), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "osqp",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = FALSE,
                                                   eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "osqp",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = FALSE,
                             eval.method = "bootstrap")
    
  )
  
})

testthat::test_that("grid search actually works, marg wass sdlp", {
  testthat::skip_on_cran(); 
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
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
    
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = solver,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins,
                             eval.method = "bootstrap"
    )
    
  )
  testthat::expect_equivalent(wsel$args$constraint,
                              list(margins = rep(1.074658, 6), 
                                   penalty = 640000), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_silent(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = solver,
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    add.joint = add.joint, add.margins = add.margins,
                                                   eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = solver,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = add.joint, add.margins = add.margins,
                             eval.method = "bootstrap")
    
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "osqp"
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
    
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = solver,
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.margins,
                             eval.method = "bootstrap")
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, 
                              list(margins = rep(1.074658, 6),
                                   penalty = 20238.58), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = solver,
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.margins,
                                                   eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  testthat::expect_warning(
    wsel <- causalOT:::wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.margins,
                             eval.method = "bootstrap")
    
  )
  
})

testthat::test_that("grid search actually works, marg wass sdlp bal", {
  testthat::skip_on_cran(); 
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
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
                             formula = formula, balance.constraints = balance.constraints,
                             eval.method = "bootstrap")
    
  )
  
  testthat::expect_equivalent(wsel$args$constraint, 
                              list(margins = rep(1.435403, 6),
                                   penalty = 20238.58), tol = 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.marginal,
                                                   formula = formula, balance.constraints = balance.constraints,
                                                   eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal,
                             eval.method = "bootstrap")
    
  )
  
})

testthat::test_that("grid search actually works, marg wass mahal bal", {
  testthat::skip_on_cran(); 
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
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
                           formula = formula, balance.constraints = balance.constraints,
                           eval.method = "bootstrap")

  )
  
  testthat::expect_lte(wsel$args$constraint$penalty - c(2605.169  ), 1e-3)
  
  estimand <- "ATC"
  # debugonce(wass_grid_search)
  testthat::expect_warning(wsel <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = method,
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 100,
                                                    add.joint = add.joint, add.margins = add.marginal,
                                                    formula = formula, balance.constraints = balance.constraints,
                                                    eval.method = "bootstrap")
  )
  
  estimand <- "ATE"
  # debugonce(wass_grid_search)
  testthat::expect_warning(
    wsel <- wass_grid_search(data, grid = NULL,
                             estimand = estimand, n.boot = 10, method = method,
                             metric = metric, p = power, solver = "mosek",
                             wass.method = "networkflow", wass.iter = 100,
                             add.joint = add.joint, add.margins = add.marginal,
                             formula = formula, balance.constraints = balance.constraints,
                             eval.method = "bootstrap")
    
  )
  
})

testthat::test_that("grid search joint.map, wass", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  estimand <- "ATT"
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- causalOT:::wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, method = "Wasserstein",
                              metric = metric, p = power, solver = solver,
                              wass.method = "networkflow", wass.iter = 0,
                              joint.mapping = TRUE, penalty = "L2",
                              eval.method = "bootstrap")
  )
  testthat::expect_equivalent(wsel2$args$constraint, list(penalty = 640,
                                                          joint = 0.02154435), 1e-3)
  
  estimand <- "ATC"
  testthat::expect_warning(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = solver,
                                                    joint.mapping = TRUE,
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    eval.method = "bootstrap"))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = solver,
                                                    joint.mapping = TRUE,
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    eval.method = "bootstrap"))
  
  
})
 
testthat::test_that("grid search scm", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
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
                              joint.mapping = TRUE, penalty = "L2",
                              eval.method = "bootstrap")
  )
  testthat::expect_lte(wsel2$args$constraint$penalty - 640, 1e-3)
  
  estimand <- "ATC"
  testthat::expect_warning(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "SCM",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    eval.method = "bootstrap"))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10, method = "SCM",
                                                    metric = metric, p = power, solver = "mosek",
                                                    wass.method = "networkflow", wass.iter = 0,
                                                    eval.method = "bootstrap"))
  
  
})

testthat::test_that("grid search joint.map crossvalidation, wass", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
   testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
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
  power <- c(4)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  estimand <- "ATT"
  
  #### get simulation functions ####
  data <- Hainmueller$new(n = n, p = p, 
                          design = design, overlap = overlap)
  data$gen_data()
  
  #don't specify grid
  # debugonce(wass_grid_search)
  testthat::expect_silent(
    wsel2 <- wass_grid_search(data, grid = NULL,
                              estimand = estimand, n.boot = 10, 
                              K = 10, R = 1, method = "Wasserstein",
                              metric = metric, p = power, solver = "mosek",
                              wass.method = "networkflow", wass.iter = 0,
                              joint.mapping = TRUE, penalty = "L2",
                              eval.method = "cross.validation")
  )
  testthat::expect_equivalent(wsel2$args$constraint, list(penalty = 20238.58,
                                                          joint = 0.4641589), 1e-3)
  
  estimand <- "ATC"
  testthat::expect_silent(wsel3 <- wass_grid_search(data, grid = NULL,
                                                    estimand = estimand, n.boot = 10,
                                                    K = 3, R = 1,
                                                    method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    joint.mapping = FALSE,
                                                    wass.method = "sinkhorn", wass.iter = 10,
                                                    eval.method = "cross.validation"))
  
  estimand <- "ATE"
  testthat::expect_silent(wsel4 <- wass_grid_search(data, grid = NULL,
                                                    K = 3, R = 1,
                                                    estimand = estimand, n.boot = 10, method = "Wasserstein",
                                                    metric = metric, p = power, solver = "mosek",
                                                    joint.mapping = FALSE,
                                                    wass.method = "sinkhorn", wass.iter = 10,
                                                    eval.method = "cross.validation"))
  
  
})
