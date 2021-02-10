testthat::test_that("ot_imputer works", {
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
  total <- rep(NA_real_, nrow(data))
  t_ind <- z==1
  c_ind <- z==0
  total[z==1] <- w1
  total[z==0] <- w0
  
  sample_weight <- list(a = w0,
                        b = w1,
                        total = renormalize(total))
  
  # debugonce(ot_imputer)
  fit_1 <- ot_imputer(formula$treated, data[t_ind,,drop = FALSE], 
                      weights = sample_weight$b)
  fit_0 <- ot_imputer(formula$control, data[c_ind,,drop = FALSE], 
                      weights = sample_weight$a)
  
  
  
  # debugonce(predict.ot_imputer)
  testthat::expect_silent(
  f_1   <- predict.ot_imputer(fit_1, data, weights = sample_weight$total,
                              verbose = FALSE, niter = 10))
  testthat::expect_silent(f_0   <- predict.ot_imputer(fit_0, data, weights = sample_weight$total,
                              verbose = FALSE, niter = 10))
  testthat::expect_true(inherits(fit_1, "ot_imputer"))
  testthat::expect_true(inherits(fit_0, "ot_imputer"))
  
})

testthat::test_that("ot_imputer coef works", {
  testthat::skip_on_cran()
  testthat::skip("Interactive only")
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
  total <- rep(NA_real_, nrow(data))
  t_ind <- z==1
  c_ind <- z==0
  total[z==1] <- w1
  total[z==0] <- w0
  
  sample_weight <- list(a = w0,
                        b = w1,
                        total = renormalize(total))
  
  # debugonce(ot_imputer)
  fit <- ot_imputer(formula[[1]], 
                   data =  cbind(data, 
                                 z = z), 
                   weights = sample_weight$total)
  
  
  # debugonce(coef.ot_imputer)
  tx_effect <- coef.ot_imputer(fit, tx.name = "z", estimand = estimand)["z"]
  
  testthat::expect_equal(c(z = .0517), tx_effect, tol = 1e-3)
  testthat::expect_true(inherits(fit, "ot_imputer"))
  
})

