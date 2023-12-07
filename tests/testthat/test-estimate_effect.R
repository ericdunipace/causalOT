testthat::test_that("bp works, ATE", {
  causalOT:::torch_check()
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
  estimand <- "ATE"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(x = original, estimand = estimand, method = "NNM")
  
  df <- data.frame(y = weights@data@y, 
                   z = weights@data@z, 
                   weights@data@x)
  bp <- barycentric_projection(formula = "y ~ .", 
                               data = df, 
                               separate.samples.on = "z", 
                               weights = weights)
  
  df0 <- df
  df1 <- df
  df0$z <- 0L
  df1$z <- 1L
  predictions0 <- predict(bp, newdata = df0, source.sample = df$z)
  predictions1 <- predict(bp, newdata = df1, source.sample = df$z)
  w     <- rep(NA_real_, nrow(df0))
  w[df$z == 1] <- weights@w1
  w[df$z == 0] <- weights@w0
  delta_bp <- predictions1 - predictions0
  
  tau_bp <- sum(delta_bp * weights@data@weights)
  
  tau_est_sep <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = TRUE)
  tau_est_tgthr <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = FALSE)
  
  testthat::expect_equal(tau_est_tgthr@estimate, tau_est_sep@estimate)
  testthat::expect_equal(tau_bp, tau_est_sep@estimate)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_0,
                         predictions0)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1,
                         predictions1)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1 - tau_est_sep@augmentedData$y_hat_0, delta_bp)
  testthat::expect_equal(causalOT:::renormalize(w),
                         causalOT:::renormalize(tau_est_sep@augmentedData$weights))
  
  
  # augmented 
  tau_aug <- estimate_effect(weights, model.function = barycentric_projection, 
                                 augmented.model = TRUE)

  y <- weights@data@y
  z <- weights@data@z
  y0<- y[z==0]
  y1<- y[z==1]
  
  testthat::expect_equal(tau_aug@estimate,
                         sum(weights@w1 * (y1 - tau_est_sep@augmentedData$y_hat_1[z==1])) -
                         sum(weights@w0 * (y0 - tau_est_sep@augmentedData$y_hat_0[z==0])) +
    mean(tau_est_sep@augmentedData$y_hat_1 - tau_est_sep@augmentedData$y_hat_0), tol = 1e-5)
  
    
})

testthat::test_that("bp works, ATT", {
  causalOT:::torch_check()
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
  estimand <- "ATT"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(x = original, estimand = estimand, method = "NNM")
  
  df <- data.frame(y = weights@data@y, 
                   z = weights@data@z, 
                   weights@data@x)
  bp <- barycentric_projection(formula = "y ~ .", 
                               data = df, 
                               separate.samples.on = "z", 
                               weights = weights)
  
  df0 <- df
  df1 <- df
  df0$z <- 0L
  df1$z <- 1L
  predictions0 <- predict(bp, newdata = df0, source.sample = df$z)
  predictions1 <- predict(bp, newdata = df1, source.sample = df$z)
  w     <- rep(NA_real_, nrow(df0))
  w[df$z == 1] <- weights@w1
  w[df$z == 0] <- weights@w0
  delta_bp <- predictions1 - predictions0
  
  tau_bp <- sum(delta_bp[df$z==1] * weights@w1)
  
  tau_est_sep <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = TRUE)
  tau_est_tgthr <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = FALSE)
  
  testthat::expect_equal(tau_est_tgthr@estimate, tau_est_sep@estimate)
  testthat::expect_equal(tau_bp, tau_est_sep@estimate)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_0,
                         predictions0)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1,
                         predictions1)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1 - tau_est_sep@augmentedData$y_hat_0, delta_bp)
  testthat::expect_equal(causalOT:::renormalize(w),
                         causalOT:::renormalize(tau_est_sep@augmentedData$weights))
  
})

testthat::test_that("bp works, ATC", {
  causalOT:::torch_check()
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
  estimand <- "ATC"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(x = original, estimand = estimand, method = "NNM")
  
  df <- data.frame(y = weights@data@y, 
                   z = weights@data@z, 
                   weights@data@x)
  bp <- barycentric_projection(formula = "y ~ .", 
                               data = df, 
                               separate.samples.on = "z", 
                               weights = weights)
  
  df0 <- df
  df1 <- df
  df0$z <- 0L
  df1$z <- 1L
  predictions0 <- predict(bp, newdata = df0, source.sample = df$z)
  predictions1 <- predict(bp, newdata = df1, source.sample = df$z)
  w     <- rep(NA_real_, nrow(df0))
  w[df$z == 1] <- weights@w1
  w[df$z == 0] <- weights@w0
  delta_bp <- predictions1 - predictions0
  
  tau_bp <- sum(delta_bp[df$z==0] * weights@w0)
  
  tau_est_sep <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = TRUE)
  tau_est_tgthr <- estimate_effect(weights, model.function = barycentric_projection, estimate.separately = FALSE)
  
  testthat::expect_equal(tau_est_tgthr@estimate, tau_est_sep@estimate)
  testthat::expect_equal(tau_bp, tau_est_sep@estimate)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_0,
                         predictions0)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1,
                         predictions1)
  testthat::expect_equal(tau_est_sep@augmentedData$y_hat_1 - tau_est_sep@augmentedData$y_hat_0, delta_bp)
  testthat::expect_equal(causalOT:::renormalize(w),
                         causalOT:::renormalize(tau_est_sep@augmentedData$weights))
  
})

testthat::test_that("estimate effect works lm, ATT", {
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATT"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights, model.function = lm)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w <- rep(NA_real_, length(z))
  w[z==0] <- w0
  w[z==1] <- w1
  w       <- w/sum(w)
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                     w1=w1, w = w), 
                                parent=environment(form))
  
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum(delta[z==1] * w1), 
                         ee@estimate)
  testthat::expect_equal(ee@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         estimate.separately = TRUE)@estimate)
  
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                      w1=w1), 
                                 parent=environment(form))
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred0[z==1]) * w1) -
                           sum((y0 - pred0[z==0]) * w0), 
                         ee2@estimate)
  testthat::expect_equal(ee2@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         augment.estimate = TRUE,
                                         estimate.separately = TRUE)@estimate)
  
  # non-augmented, joint
  ee <- estimate_effect(weights, model.function = lm,
                        estimate.separately = FALSE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==0] <- w0
  w[z==1] <- w1
  w       <- w/sum(w)
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                     w1=w1, w = w), 
                                parent=environment(form))
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum(delta[z==1] * w1), 
                         ee@estimate)
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         estimate.separately = FALSE,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                     w1=w1, w = w), 
                                parent=environment(form))
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred0[z==1]) * w1) -
                           sum((y0 - pred0[z==0]) * w0), 
                         ee2@estimate)
  
})

testthat::test_that("estimate effect works lm, ATC", {
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATC"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights, model.function = lm)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, nrow(x))
  w[z==1] <- w1
  w[z==0] <- w0
  w       <- w/sum(w)
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                     w1=w1, w = w), 
                                parent=environment(form))
  
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum(delta[z==0] * w0), 
                         ee@estimate)
  testthat::expect_equal(ee@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         estimate.separately = TRUE)@estimate)
  
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred1[z==1]) * w1) -
                           sum((y0 - pred1[z==0]) * w0), 
                         ee2@estimate)
  testthat::expect_equal(ee2@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         augment.estimate = TRUE,
                                         estimate.separately = TRUE)@estimate)
  
  # non-augmented, joint
  ee <- estimate_effect(weights, model.function = lm,
                        estimate.separately = FALSE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==0] <- w0
  w[z==1] <- w1
  
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum(delta[z==0] * w0), 
                         ee@estimate)
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         estimate.separately = FALSE,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  
  
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred0[z==1]) * w1) -
                           sum((y0 - pred0[z==0]) * w0), 
                         ee2@estimate,
                         tol = 1e-5)
  
})

testthat::test_that("estimate effect works lm, ATE", {
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATE"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights, model.function = lm)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==1] <- w1
  w[z==0] <- w0
  w  <- w/sum(w)
  
  form <- formula("y ~ .")
  environment(form) <- list2env(list(w0=w0,
                                     w1=w1, w = w), 
                                parent=environment(form))
  
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(weighted.mean(delta, 
                                       w = weights@data@weights), 
                         ee@estimate)
  testthat::expect_equal(ee@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         estimate.separately = TRUE)@estimate)
  
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==1] <- w1
  w[z==0] <- w0
  w  <- w/sum(w)
  
  fit0 <- lm(form, data = data.frame(x0, y = y0), weights = w0)
  fit1 <- lm(form, data = data.frame(x1, y = y1), weights = w1)
  
  pred1 <- predict(fit1, newdata = data.frame(x))
  pred0 <- predict(fit0, newdata = data.frame(x))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred1[z==1]) * w1) -
                           sum((y0 - pred0[z==0]) * w0) + 
                           weighted.mean(delta, 
                                         w = weights@data@weights), 
                         ee2@estimate)
  testthat::expect_equal(ee2@estimate,
                         estimate_effect(weights, 
                                         model.function = lm, 
                                         augment.estimate = TRUE,
                                         estimate.separately = TRUE)@estimate)
  
  # non-augmented, joint
  ee <- estimate_effect(weights, model.function = lm,
                        estimate.separately = FALSE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==0] <- w0
  w[z==1] <- w1
  
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum(delta[z==0] * w0), 
                         ee@estimate)
  
  # augmented, separate
  ee2 <- estimate_effect(weights, model.function = lm,
                         estimate.separately = FALSE,
                         augment.estimate = TRUE)
  
  x  <- weights@data@x
  z  <- weights@data@z
  y  <- weights@data@y
  
  x0 <- x[z==0,]
  x1 <- x[z==1,]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights@w0
  w1 <- weights@w1
  w  <- rep(NA_real_, length(z))
  w[z==1] <- w1
  w[z==0] <- w0
  w  <- w/sum(w)
  
  fit <- lm(form, data = data.frame(x, z = z, y = y), weights = w)
  
  pred1 <- predict(fit, newdata = data.frame(x, z = 1L))
  pred0 <- predict(fit, newdata = data.frame(x, z = 0L))
  
  delta <- pred1 - pred0
  
  testthat::expect_equal(sum((y1 - pred1[z==1]) * w1) -
                           sum((y0 - pred0[z==0]) * w0) + 
                           weighted.mean(delta, 
                                         w = weights@data@weights), 
                         ee2@estimate)
  
})


testthat::test_that("ATT give proper var",{
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATT"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights)
  v <- vcov(ee)
  
  w0 <- weights@w0
  w1 <- weights@w1
  n1 <- length(w1)
  n0 <- length(w0)
  
  y <- weights@data@y
  z <- weights@data@z
  
  mu1 <- sum(y[z==1] * w1)
  mu0 <- sum(y[z==0] * w0)
  testthat::expect_equal(coef(ee), c(estimate = mu1-mu0))
  
  v1 <- sum((w1 * n1 * (y[z==1] - mu1))^2)/(n1-1)
  v0 <- sum((w0 * n1 * (y[z==0] - mu0))^2)/(n1-1)
  testthat::expect_equal(v1, 
  var(y[z==1]) )
  testthat::expect_equal(v0/n1 + v1/n1, as.numeric(v))
})

testthat::test_that("ATT give proper var lm",{
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATT"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights, model.function = lm, augment = TRUE)
  v <- vcov(ee)
  
  w0 <- weights@w0
  w1 <- weights@w1
  n1 <- length(w1)
  n0 <- length(w0)
  
  y <- weights@data@y
  z <- weights@data@z
  
  mu1 <- ee@augmentedData$y_hat_1
  mu0 <- ee@augmentedData$y_hat_0
  tau <- sum(w1 * (y - mu1)[z==1]) - sum(w0 * (y - mu0)[z==0]) +
    sum(w1 * (mu1 - mu0)[z==1])
  testthat::expect_equal(coef(ee), c(estimate = tau))
  
  v1 <- sum((w1 * n1 * (y - mu1)[z==1])^2)/(n1-1)
  v0 <- sum((w0 * n1 * (y - mu0)[z==0])^2)/(n1-1)
  vm <- sum(((mu1 - mu0 - coef(ee)) * z)^2)/(n1-1)
  testthat::expect_equal(v0/n1 + v1/n1 + vm/n1, as.numeric(v))
})

testthat::test_that("ATC give proper var",{
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATC"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights)
  v <- vcov(ee)
  
  w0 <- weights@w0
  w1 <- weights@w1
  n1 <- length(w1)
  n0 <- length(w0)
  
  y <- weights@data@y
  z <- weights@data@z
  
  mu1 <- sum(y[z==1] * w1)
  mu0 <- sum(y[z==0] * w0)
  testthat::expect_equal(coef(ee), c(estimate = mu1-mu0))
  
  v1 <- sum((w1 * n0 * (y[z==1] - mu1))^2)/(n0-1)
  v0 <- sum((w0 * n0 * (y[z==0] - mu0))^2)/(n0-1)
  testthat::expect_equal(v0, 
                         var(y[z==0]) )
  testthat::expect_equal(v0/n0 + v1/n0, as.numeric(v))
})

testthat::test_that("ATC give proper var lm",{
  causalOT:::torch_check()
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  estimand <- "ATC"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  original$gen_data()
  weights <- calc_weight(original, estimand = estimand, method = "NNM")
  
  # non-augmented, separate
  ee <- estimate_effect(weights, model.function = lm, augment = TRUE)
  v <- vcov(ee)
  
  w0 <- weights@w0
  w1 <- weights@w1
  n1 <- length(w1)
  n0 <- length(w0)
  
  y <- weights@data@y
  z <- weights@data@z
  
  mu1 <- ee@augmentedData$y_hat_1
  mu0 <- ee@augmentedData$y_hat_0
  tau <- sum(w1 * (y - mu1)[z==1]) - sum(w0 * (y - mu0)[z==0]) +
    sum(w0 * (mu1 - mu0)[z==0])
  testthat::expect_equal(coef(ee), c(estimate = tau))
  
  v1 <- sum((w1 * n0 * (y - mu1)[z==1])^2)/(n0-1)
  v0 <- sum((w0 * n0 * (y - mu0)[z==0])^2)/(n0-1)
  vm <- sum(((mu1 - mu0 - coef(ee)) * (1-z))^2)/(n0-1)
  testthat::expect_equal(v0/n0 + v1/n0 + vm/n0, as.numeric(v))
})
