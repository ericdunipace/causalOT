testthat::test_that("barycentric_projection works, p = 2 tensor", {
  causalOT:::torch_check()
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^5
  pp <- 6
  overlap <- "low"
  design <- "A"
  estimate <- "ATE"
  power <- 2
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = pp,
                                    design = design, overlap = overlap)
  data$gen_data()
  weights <- causalOT::calc_weight(x = data,
                         z = NULL, y = NULL,
                         estimand = estimate,
                         method = "NNM")
  df <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_silent(fit <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights))
  testthat::expect_equal(fit$a, weights@w0)
  testthat::expect_equal(fit$b, weights@w1)
  testthat::expect_equal(fit$p, power)
  
  
  # check fit works if provide  weights as a vector
  w <- c(weights@w0, weights@w1)[order(order(data$get_z()))]
  fit2 <- causalOT:::barycentric_projection(y ~ ., data = df, weight = w)
  fit3 <- fit
  fit3$data@weights <- w
  testthat::expect_equal(fit3, fit2)
  
  # for same data
  preds <- causalOT:::predict.bp(fit)
  testthat::expect_equal(preds, predict(fit)) # make sure S3 registered
  
  # compare to online formulation
  causalOT:::rkeops_check()
  if(packageVersion("rkeops") >= "2.0" ) {
    rkeops::rkeops_use_float64()
  } else {
    rkeops::compile4float64()
  }
  mess <- testthat::capture_output(fito <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, cost.online = "online"))
  mess <- testthat::capture_output(predo <- predict(fito))
  testthat::expect_equal(preds, predo, tol = 1e-5)
  
  # compare to manual est
  y_hat <- rep(NA_real_, causalOT:::get_n(fit$data))
  y0 <- causalOT:::get_y0(fit$data)
  y1 <- causalOT:::get_y1(fit$data)
  x0 <- causalOT:::get_x0(fit$data)
  x1 <- causalOT:::get_x1(fit$data)
  z  <- causalOT:::get_z(fit$data)
  c_pot <- fit$potentials$f_aa
  t_pot <- fit$potentials$g_bb
  C_00  <- causalOT:::cost(x0,x0, power, cost_function = fit$cost_function)
  C_11  <- causalOT:::cost(x1,x1, power, cost_function = fit$cost_function)
  lambda <- fit$penalty
  w0_log  <- causalOT:::log_weights(weights@w0)
  w1_log  <- causalOT:::log_weights(weights@w1)
  
  eta_00 <- c_pot/lambda + w0_log - C_00$data/lambda
  y_hat[z == 0] <- as.matrix(eta_00$log_softmax(2)$exp()) %*% y0
  
  eta_11 <- t_pot/lambda + w1_log - C_11$data/lambda
  y_hat[z == 1] <- as.matrix(eta_11$log_softmax(2)$exp()) %*% y1
  
  testthat::expect_equal(y_hat, preds)
  
  # for new 
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  preds2 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z())
  testthat::expect_equal(preds, preds2, tol = 1e-5)
  preds2o <- causalOT:::predict.bp(fito, newdata=df.new, source.sample = data$get_z())
  testthat::expect_equal(preds, preds2o, tol = 1e-5)
  
  
  
  data$gen_data()
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_silent(preds3 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z()))
  
  fit_debias <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, debias = TRUE)
  predict_debias <- predict(fit_debias)
  testthat::expect_equivalent(df$y, predict_debias)
  
})

testthat::test_that("barycentric_projection works, p = 1", {
  causalOT:::torch_check()
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^5
  p <- 6
  power <- 1
  overlap <- "low"
  design <- "A"
  estimate <- "ATE"
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
                                    design = design, overlap = overlap)
  data$gen_data()
  weights <- causalOT::calc_weight(x = data,
                                   z = NULL, y = NULL,
                                   estimand = estimate,
                                   method = "NNM")
  df <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_silent(fit <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power))
  testthat::expect_equal(fit$a, weights@w0)
  testthat::expect_equal(fit$b, weights@w1)
  testthat::expect_equal(fit$p, power)
  
  # check fit works if provide  weights as a vector
  w <- c(weights@w0, weights@w1)[order(order(data$get_z()))]
  fit2 <- causalOT:::barycentric_projection(y ~ ., data = df, weight = w, p = power)
  fit3 <- fit
  fit3$data@weights <- w
  testthat::expect_equal(fit3, fit2)
  
  # for same data
  preds <- causalOT:::predict.bp(fit)
  testthat::expect_equal(preds, predict(fit)) # make sure S3 registered
  
  causalOT:::rkeops_check()
  # give error for the p = 1 online
  if(packageVersion("rkeops") >= "2.0" ) {
    rkeops::rkeops_use_float64()
  } else {
    rkeops::compile4float64()
  }
  testthat::expect_warning(fito <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, cost.online = "online"))
  testthat::expect_error(predo <- predict(fito))
  fun <- function(x1, x2, p) {
              ((1/p) * torch::torch_cdist(x1 = torch::torch_tensor(x1, dtype = torch::torch_double()), 
                                          x2 = torch::torch_tensor(x2, dtype = torch::torch_double()), 
                                          p = p)^p)$contiguous()
  }
  predo <- predict(fito, cost_function = fun)
  testthat::expect_equal(preds, predo)
  
  # compare to manual est
  y_hat <- rep(NA_real_, causalOT:::get_n(fit$data))
  y0 <- causalOT:::get_y0(fit$data)
  y1 <- causalOT:::get_y1(fit$data)
  x0 <- causalOT:::get_x0(fit$data)
  x1 <- causalOT:::get_x1(fit$data)
  z  <- causalOT:::get_z(fit$data)
  c_pot <- fit$potentials$f_aa
  t_pot <- fit$potentials$g_bb
  C_00  <- causalOT:::cost(x0,x0, power, cost_function = fit$cost_function)
  C_11  <- causalOT:::cost(x1,x1, power, cost_function = fit$cost_function)
  lambda <- fit$penalty
  w0_log  <- causalOT:::log_weights(weights@w0)
  w1_log  <- causalOT:::log_weights(weights@w1)
  
  eta_00 <- c_pot/lambda + w0_log - C_00$data/lambda
  y_hat[z == 0] <- apply(as.matrix(eta_00$log_softmax(2)$exp()), 1, function(w) matrixStats::weightedMedian(x = y0, w = w))
  
  eta_11 <- t_pot/lambda + w1_log - C_11$data/lambda
  y_hat[z == 1] <- apply(as.matrix(eta_11$log_softmax(2)$exp()), 1, function(w) matrixStats::weightedMedian(x = y1, w = w))
  
  testthat::expect_equal(y_hat, preds)
  
  # for new 
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  preds2 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z())
  testthat::expect_equal(preds, preds2)
  
  data$gen_data()
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_silent(preds3 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z()))
  
  fit_debias <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, debias = TRUE)
  predict_debias <- predict(fit_debias)
  testthat::expect_equivalent(df$y, predict_debias)
  
})

testthat::test_that("barycentric_projection works, p = 3", {
  causalOT:::torch_check()
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^5
  p <- 6
  power <- 3L
  overlap <- "low"
  design <- "A"
  estimate <- "ATE"
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
                                    design = design, overlap = overlap)
  data$gen_data()
  weights <- causalOT::calc_weight(x = data,
                                   z = NULL, y = NULL,
                                   estimand = estimate,
                                   method = "NNM")
  df <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_silent(fit <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power))
  testthat::expect_equal(fit$a, weights@w0)
  testthat::expect_equal(fit$b, weights@w1)
  testthat::expect_equal(fit$p, power)
  
  # check fit works if provide  weights as a vector
  w <- c(weights@w0, weights@w1)[order(order(data$get_z()))]
  fit2 <- testthat::expect_silent(causalOT:::barycentric_projection(y ~ ., data = df, weight = w, p = power))
  fit3 <- fit
  fit3$data@weights <- w
  testthat::expect_equal(fit3, fit2)
  
  # for same data
  testthat::expect_warning(preds <- causalOT:::predict.bp(fit))
  testthat::expect_warning(preds2 <- predict(fit))
  testthat::expect_equal(preds, preds2) # make sure S3 registered
  
  causalOT:::rkeops_check()
  if(packageVersion("rkeops") >= "2.0" ) {
    rkeops::rkeops_use_float64()
  } else {
    rkeops::compile4float64()
  }
  mess <- testthat::capture_output((fito <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, cost.online = "online")))
  testthat::expect_error(causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = as.numeric(power), cost.online = "online"))
  mess <- testthat::capture_output(testthat::expect_warning(predo <- predict(fito)))
  testthat::expect_equal(preds, predo, tol = 1e-5)
  
  # compare to manual est
  closure_test <- function() {
    opt$zero_grad()
    n <- length(y_source)
    m <- length(y_targ)
    dist_mat <- fit$cost_function(y_source$view(c(n,1)), y_targ$view(c(m,1)), power)
    col_loss <- (dist_mat * wt)$sum(2)
    # tot_loss <- col_loss$dot(exp(g/eps + b_log))
    tot_loss <- col_loss$sum()
    tot_loss$backward()
    tot_loss
  }
  y_hat <- rep(NA_real_, causalOT:::get_n(fit$data))
  y0 <- causalOT:::get_y0(fit$data)
  y1 <- causalOT:::get_y1(fit$data)
  x0 <- causalOT:::get_x0(fit$data)
  x1 <- causalOT:::get_x1(fit$data)
  z  <- causalOT:::get_z(fit$data)
  c_pot <- fit$potentials$f_aa
  t_pot <- fit$potentials$g_bb
  C_00  <- causalOT:::cost(x0,x0, power, cost_function = fit$cost_function)
  C_11  <- causalOT:::cost(x1,x1, power, cost_function = fit$cost_function)
  eps <- fit$penalty
  w0_log  <- causalOT:::log_weights(weights@w0)
  w1_log  <- causalOT:::log_weights(weights@w1)
  
  eta_00 <- (c_pot/eps + w0_log) - C_00$data/eps
  wt <- eta_00$log_softmax(2)$exp()
  b_log  <- torch::torch_tensor(w0_log, dtype = torch::torch_double())
  g  <-  torch::torch_tensor(c_pot, dtype = torch::torch_double())
  y_targ <- torch::torch_tensor(y0, dtype = torch::torch_double())
  y_source <- torch::torch_zeros_like(y_targ, dtype = torch::torch_double(), requires_grad = TRUE)
  opt <- torch::optim_lbfgs(y_source)
  
  for(i in 1:5) {
    loss <- opt$step(closure_test)
    # print(loss)
  }
  
  y_hat[z == 0] <- as.numeric(y_source$detach())
  
  eta_11 <- t_pot/eps + w1_log - C_11$data/eps
  wt <- eta_11$log_softmax(2)$exp()
  b_log  <- torch::torch_tensor(w1_log, dtype = torch::torch_double())
  g  <-  torch::torch_tensor(t_pot, dtype = torch::torch_double())
  y_targ <- torch::torch_tensor(y1, dtype = torch::torch_double())
  y_source <- torch::torch_zeros_like(y_targ, dtype = torch::torch_double(), requires_grad = TRUE)
  opt <- torch::optim_lbfgs(y_source)
  
  for(i in 1:5) {
    loss <- opt$step(closure_test)
    # print(loss)
  }
  y_hat[z == 1] <- as.numeric(y_source$detach())
  
  testthat::expect_equal(y_hat, preds)
  
  # for new 
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_warning(preds2 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z()))
  testthat::expect_equal(preds, preds2)
  
  data$gen_data()
  df.new <- data.frame(y = data$get_y(), z = data$get_z(), data$get_x())
  testthat::expect_warning(preds3 <- causalOT:::predict.bp(fit, newdata=df.new, source.sample = data$get_z()))
  
  fit_debias <- causalOT:::barycentric_projection(y ~ ., data = df, weight = weights, p = power, debias = TRUE)
  predict_debias <- predict(fit_debias)
  testthat::expect_equivalent(df$y, predict_debias)
})

