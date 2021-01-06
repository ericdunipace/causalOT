testthat::test_that("mapping works, ATT", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  estimand <- "ATT"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM", transport.matrix = TRUE)
  
  pd <- causalOT:::prep_data(original)
  data <- pd$df
  z    <- pd$z
  model.fun <- lm
  formula <- list(treated = "y ~.",
                  control = "y ~.")
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z==1
  c_ind <- z==0
 
  fit_1 <- model.fun(formula$treated, data[t_ind,,drop=FALSE])
  fit_0 <- model.fun(formula$control, data[c_ind,,drop=FALSE])
  f_1   <- predict(fit_1, data)
  f_0   <- predict(fit_0, data)
  
  sw <- get_sample_weight(NULL, z)
  
  # debugonce(causalOT:::mapping)
  testthat::expect_warning(ys <- causalOT:::mapping(data = data, z = z, weights = weights, estimand = estimand, 
                     f1 = f_1, f0 = f_0, sw = sw))
  
  non_map <- causalOT:::.outcome_calc_deprecated(data, z, weights, formula, model.fun, match = TRUE,
                                                 estimand = estimand)
  
  testthat::expect_equal(mean(ys$y1 - ys$y0), non_map)
})

testthat::test_that("mapping works, ATC", {
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
  estimand <- "ATC"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, metric = metric,
                         method = "NNM", transport.matrix = TRUE,
                         p = power)
  
  pd <- causalOT:::prep_data(original)
  data <- pd$df
  z    <- pd$z
  model.fun <- lm
  formula <- list(treated = "y ~.",
                  control = "y ~.")
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z==1
  c_ind <- z==0
  
  fit_1 <- model.fun(formula$treated, data[t_ind,,drop=FALSE])
  fit_0 <- model.fun(formula$control, data[c_ind,,drop=FALSE])
  f_1   <- predict(fit_1, data)
  f_0   <- predict(fit_0, data)
  
  sw <- get_sample_weight(NULL, z)
  
  # debugonce(causalOT:::mapping)
  ys <- causalOT:::mapping(data = data, z = z, weights = weights, estimand = estimand, 
                           f1 = f_1, f0 = f_0, sw = sw)
  #map y0 0.7493255 y1 2.460151
  #f0 1.197229  f1 2.42858
  #resid y0 -0.4479032 y1 0.05772957
  
  # debugonce(causalOT:::.outcome_calc_deprecated)
  non_map <- causalOT:::.outcome_calc_deprecated(data, z, weights, formula, model.fun, match = TRUE,
                                                 estimand = estimand)
  #map y0 0.7493255 y1 2.498858
  #f0 1.197229 f1 2.429185
  #resid y0 -0.4479032 y1 0.06967362
  testthat::expect_equal(mean(ys$y1 - ys$y0), non_map)
  
})

testthat::test_that("mapping works, ATE", {
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
  estimand <- "ATE"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, metric = metric,
                         method = "NNM", transport.matrix = TRUE,
                         p = power)
  
  pd <- causalOT:::prep_data(original)
  data <- pd$df
  z    <- pd$z
  model.fun <- lm
  formula <- list(treated = "y ~.",
                  control = "y ~.")
  
  w0 <- weights$w0
  w1 <- weights$w1
  t_ind <- z==1
  c_ind <- z==0
  
  fit_1 <- model.fun(formula$treated, data[t_ind,,drop=FALSE])
  fit_0 <- model.fun(formula$control, data[c_ind,,drop=FALSE])
  f_1   <- predict(fit_1, data)
  f_0   <- predict(fit_0, data)
  
  sw <- get_sample_weight(NULL, z)
  
  # debugonce(causalOT:::mapping)
  ys <- causalOT:::mapping(data = data, z = z, weights = weights, estimand = estimand, 
                           f1 = f_1, f0 = f_0, sw = sw)
  
  ysATT <- causalOT:::mapping(data = data, z = z, weights = weights, estimand = "ATT", 
                           f1 = f_1, f0 = f_0, sw = sw)
  ysATC <- causalOT:::mapping(data = data, z = z, weights = weights, estimand = "ATC", 
                              f1 = f_1, f0 = f_0, sw = sw)
  yscate <- (sum(ysATT$y1 - ysATT$y0) + sum(ysATC$y1 - ysATC$y0))/length(z)
  #orig y0 0.7493255 y1 2.93044
  #map y0 1.451089 y1 2.460151
  #f0 y0 1.591257  y1 2.210824
  #f1 y0 1.883628  y1 2.42858
  #resid y0 -0.3851154 y1 0.03190556
  
  # debugonce(causalOT:::.outcome_calc_deprecated)
  non_map <- causalOT:::.outcome_calc_deprecated(data, z, weights, formula, model.fun, match = TRUE,
                                                 estimand = estimand)
  #map y0 1.451089 y1 2.460151
  #f0 y0 1.591257  y1 2.210824
  #f1 y0 1.197229  y1 2.42858
  #resid y0 -0.3216016 y1 0.01760001
  testthat::expect_equal(mean(ys$y1-ys$y0), yscate)
  testthat::expect_equal(yscate, mean(f_1) - mean(f_0) + sum((data$y - f_1)[z==1]*weights$w1)- sum((data$y - f_0)[z==0]*weights$w0))
  
  
})
