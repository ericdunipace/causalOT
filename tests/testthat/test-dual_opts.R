testthat::test_that("test forward functions", {
  causalOT:::torch_check()
  
  n <- 256
  
  z  <- matrix(rnorm(n/2*2), n/2, 2) + matrix(c(0,.5), n/2,2, byrow = TRUE)
  x  <- matrix(rnorm(n *2), n, 2)
  
  
  m1 <- Measure(x, target.values = colMeans(z), adapt = "weights")
  mt <- Measure(z)
  
  gamma <- torch::torch_tensor(stats::rnorm(n), 
                               device = m1$device,
                               dtype = m1$dtype)
  
  ot_tens <- causalOT:::OT$new(x = x, y = z, debias = TRUE, tensorized = "tensorized", penalty = 10)
  
  C_xy <- ot_tens$C_xy$data
  C_xx <- ot_tens$C_xx$data
  
  a_log<- causalOT:::log_weights(ot_tens$a)
  b_log<- causalOT:::log_weights(ot_tens$b)
  lambda <- ot_tens$penalty
  delta <- 0.01
  
  dual_forwards <- torch::jit_compile(causalOT:::dual_forward_code_tensorized)
  
  a1_script <- dual_forwards$calc_w1(gamma$detach(), C_xy, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  a2_script <- dual_forwards$calc_w2(gamma$detach(), C_xx, a_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  g    <- b_log -  ((gamma$detach() + a_log)$detach()$view(c(n,1))-C_xy/lambda)$logsumexp(1)
  K    <- (gamma$detach() + a_log)$view(c(n,1)) + g - C_xy/lambda
  a1   <- (K )$logsumexp(2)$exp()$detach()
  a1   <- as.numeric((a1/a1$sum())$to(device = "cpu"))
  
  testthat::expect_equal(a1, as.numeric(a1_script$to(device = "cpu")), label = "calc_w1")
  
  f_star <- gamma$detach() + a_log
  K2   <- (f_star$view(c(n,1)) + f_star -C_xx/lambda)
  norm  <-  K2$view(c(n*n,1))$logsumexp(1)
  a2     <- as.numeric((K2 - norm)$logsumexp(1)$exp()$detach()$to(device = "cpu"))
  
  testthat::expect_equal(a2, as.numeric(a2_script$to(device = "cpu")), label = "calc_w2")
  
  testthat::expect_equal(gamma$dot(a1_script-a2_script)$item() * - 1,
                         dual_forwards$cot_dual(gamma$detach(), C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))$loss$item(), label = "loss calc")
  
  beta1 <- torch::torch_tensor(stats::rnorm(2), 
                               device = gamma$device,
                               dtype = gamma$dtype)
  
  f_prime <- gamma$detach() + a_log #- m1$balance_functions$matmul(beta1_det)
  beta1_det <- beta1$detach()
  g    <- b_log -  (f_prime$view(c(n,1))-C_xy/lambda)$logsumexp(1)
  K    <- (f_prime$view(c(n,1)) + g - C_xy/lambda)
  a1   <- (K )$logsumexp(2)$exp()$detach()
  a1   <- as.numeric((a1)$to(device = "cpu"))

  f_star <- gamma$detach() + a_log#- m1$balance_functions$matmul(beta2_det)
  K2   <- (f_star$view(c(n,1)) + f_star -C_xx/lambda)
  norm  <-  K2$view(c(n*n,1))$logsumexp(1)
  a2     <- as.numeric((K2 - norm)$logsumexp(1)$log_softmax(1)$exp()$detach()$to(device = "cpu"))
  
  testthat::expect_equal(a1, as.numeric(a1_script$to(device = "cpu")), label = "calc_w1")
  testthat::expect_equal(a2, as.numeric(a2_script$to(device = "cpu")), label = "calc_w2")
  
  res <- dual_forwards$cot_dual(gamma$detach(), C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))

  loss_gamma <- gamma$dot(a1_script-a2_script)$item()
  
  testthat::expect_equal(loss_gamma * - 1,
                         res$loss$item(), label = "loss calc",
                         tol = 1e-5)
  
  
  testthat::expect_equal( (a1_script-a2_script)$norm()$item(),
                          res$avg_diff$item(),
                          tol = 1e-5)
  
  testthat::expect_equal(res$bf_diff$item(), 0.0,
                         tol = 1e-5)

  
  
  diff1 <- (m1$balance_functions$transpose(2,1)$matmul(a1_script) - m1$balance_target)
  beta_check1 <- diff1 * beta1$detach() - delta * beta1$detach()$abs()

  loss_beta = diff1$dot(beta1) - delta * beta1$abs()$sum()
  
  loss <- loss_gamma + loss_beta

  loss$multiply_(-1.0) #to make min
  
  
  res2 <- dual_forwards$cot_bf_dual(gamma, C_xy, C_xx, a_log,
                                    b_log,
                                   torch::jit_scalar(lambda),
                                   torch::jit_scalar(as.integer(n)),
                                   beta1, m1$balance_functions,
                                   m1$balance_target,
                                   torch::jit_scalar(delta))
  
  testthat::expect_equal(loss$item(), res2$loss$item(), tol = 1e-4)
  testthat::expect_equal(res2$avg_diff$item(), res$avg_diff$item())
  testthat::expect_equal(res2$bf_diff$item(), diff1$abs()$max()$item())
  
  
  # test keops versions
  causalOT:::rkeops_check()
  
  ot_keops <- causalOT:::OT$new(x = x, y = z, debias = TRUE, tensorized = "online", penalty = 10)
  
  C_xy <- ot_keops$C_xy
  C_xx <- ot_keops$C_xx
  
  keops_fun <- causalOT:::dual_forwards_keops
  
  a1_script <- keops_fun$calc_w1(gamma$detach(), C_xy, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  a2_script <- keops_fun$calc_w2(gamma$detach(), C_xx, a_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  res_keops <- keops_fun$cot_dual(
    gamma, C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n))
  )
  res_keops_2 <- keops_fun$cot_bf_dual(
    gamma, C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)),
    beta1, m1$balance_functions,
    m1$balance_target,
    torch::jit_scalar(delta)
  )
  
  testthat::expect_equal(as.numeric(a1_script$to(device = "cpu")), a1,
                         tol = 1e-3)
  testthat::expect_equal(as.numeric(a2_script$to(device = "cpu")), a2,
                         tol = 1e-3)  
  
  testthat::expect_equal(loss_gamma * -1,  res_keops$loss$item(),
                         tol = 1e-5)
  testthat::expect_equal(loss$item(),  res_keops_2$loss$item(), tol = 1e-5)
  testthat::expect_equal(diff1$abs()$max()$item(), res_keops_2$bf_diff$item(), tol = 1e-5)
  testthat::expect_equal(as.numeric(res2$beta_check$to(device = "cpu")),
                         as.numeric(res_keops_2$beta_check$to(device = "cpu")), tol = 1e-5 )

  
})

testthat::test_that("dual nn modules work as expected",{
  causalOT:::torch_check()
  set.seed(1231)
  
  n <- 256
  
  z  <- matrix(rnorm(n/2*2), n/2, 2) + matrix(c(0,.5), n/2,2, byrow = TRUE)
  x  <- matrix(rnorm(n *2), n, 2)
  
  
  m1 <- Measure(x, target.values = colMeans(z), adapt = "weights")
  mt <- Measure(z)
  
  opt <- causalOT:::cotDualOpt$new(n, 2)
  
  gamma <- torch::torch_tensor(stats::rnorm(n), 
                               device = m1$device,
                               dtype = m1$dtype)
  
  torch::with_no_grad(opt$gamma$copy_(gamma))
  
  ot_tens <- causalOT:::OT$new(x = x, y = z, debias = TRUE, tensorized = "tensorized", penalty = 10)
  
  C_xy <- ot_tens$C_xy
  C_xx <- ot_tens$C_xx
  
  a_log<- causalOT:::log_weights(ot_tens$a)
  b_log<- causalOT:::log_weights(ot_tens$b)
  lambda <- ot_tens$penalty
  delta <- 0.01
  
  dual_forwards <- torch::jit_compile(causalOT:::dual_forward_code_tensorized)
  
  res <- dual_forwards$cot_dual(gamma$detach(), C_xy$data, C_xx$data, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  res_mod <- opt$forward(C_xy, C_xx, a_log, b_log, lambda)
  
  tests <- function(res, res_mod, opt, gamma) {
    testthat::expect_equal(res, res_mod)
    testthat::expect_equal(res$loss$item(), res_mod$loss$item(),
                           tol = 1e-5) 
    testthat::expect_equal(res$avg_diff$item(), res_mod$avg_diff$item(),
                           tol = 1e-5)
    testthat::expect_equal(res$bf_diff$item(), res_mod$bf_diff$item(),
                           tol = 1e-5)
    
    
    param <- opt$clone_param()
    testthat::expect_equal(as.numeric(param$gamma$to(device = "cpu")), as.numeric(gamma$to(device = "cpu")), tol = 1e-5)
    testthat::expect_true(param$gamma$requires_grad == FALSE)
    testthat::expect_true(opt$gamma$requires_grad == TRUE)
    
    # test convergence function
    param$gamma<- param$gamma * 0.0
    testthat::expect_true(isFALSE(opt$converged(res_mod,
                                                1e-5, 1e-6, param,
                                                tol = 1e-8, lambda, delta)))
                          
    testthat::expect_true(opt$converged(res_mod,
                                        1e-5, 1e-6, param,
                                        tol = 300, lambda, delta)
                                                
                          )
                          
  }
  
  tests(res, res_mod, opt, gamma)
  
  # bf
  optbf <- causalOT:::cotDualBfOpt$new(n,2)
  
  torch::with_no_grad({
    optbf$beta$copy_(c(1,2))
    optbf$gamma$copy_(gamma)
    }
    )
  beta1 <- optbf$beta$detach()$clone()
  res <- dual_forwards$cot_bf_dual(
    gamma$detach() - m1$balance_functions$matmul(beta1), 
                                   C_xy$data, C_xx$data, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)),
                                   beta1, m1$balance_functions,
                                   m1$balance_target,
                                   torch::jit_scalar(delta)
                                   )
  
  res_mod <- optbf$forward(C_xy, C_xx, a_log, b_log, lambda,
                           m1$balance_functions,
                           m1$balance_target,
                           torch::jit_scalar(delta))
  
  tests(res, res_mod, optbf,  gamma -  m1$balance_functions$matmul(beta1))
  
  #### check keops opt ####
  
  causalOT:::rkeops_check()
  
  ot_keops <- causalOT:::OT$new(x = x, y = z, debias = TRUE, tensorized = "online", penalty = 10)
  
  C_xy <- ot_keops$C_xy
  C_xx <- ot_keops$C_xx
  
  opt <- causalOT:::cotDualOpt_keops$new(n, 2)
  torch::with_no_grad(opt$gamma$copy_(gamma))
  
  keops_fun <- causalOT:::dual_forwards_keops
  
  res <- keops_fun$cot_dual(gamma$detach(), C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)))
  
  res_mod <- opt$forward(C_xy, C_xx, a_log, b_log, lambda)
  tests(res, res_mod, opt, gamma)
  
  
  optbf <- causalOT:::cotDualBfOpt_keops$new(n, 2)
  torch::with_no_grad({
    optbf$gamma$copy_(gamma)
    optbf$beta$copy_(beta1)
    })
  
  res <- keops_fun$cot_bf_dual(
    gamma$detach() - m1$balance_functions$matmul(beta1), 
    C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), torch::jit_scalar(as.integer(n)),
    beta1, m1$balance_functions,
    m1$balance_target,
    torch::jit_scalar(delta)
  )
  
  res_mod <- optbf$forward(C_xy, C_xx, a_log, b_log, lambda,
                           m1$balance_functions,
                           m1$balance_target,
                           torch::jit_scalar(delta))
  tests(res, res_mod, optbf, gamma - m1$balance_functions$matmul(beta1))
  testthat::expect_true(all(as.logical((optbf$beta == c(1,2))$to(device = "cpu"))))
  
})

testthat::test_that("training function works for dual optimizer",{
  causalOT:::torch_check()
  set.seed(1231)
  
  n <- 256
  
  z  <- matrix(rnorm(n/2*2), n/2, 2) + matrix(c(0,.5), n/2,2, byrow = TRUE)
  x  <- matrix(rnorm(n *2), n, 2)
  
  
  m1 <- Measure(x, target.values = colMeans(z), adapt = "weights")
  mt <- Measure(z)
  
  cot <- causalOT:::cotDualTrain$new(m1,mt)
  otp <- OTProblem(m1,mt)
  
  cot_names <- names(formals(cot$setup_arguments))
  otp_names <- names(formals(otp$setup_arguments))
  
  testthat::expect_equal(cot_names, otp_names)
  
  # test that setup arg makes correct nn_holder
  testthat::expect_silent(cot$setup_arguments())
  testthat::expect_silent(otp$setup_arguments())
  
  testthat::expect_true(inherits(cot$.__enclos_env__$private$nn_holder, "cotDualBfOpt"))
  
  testthat::expect_true(length(cot$.__enclos_env__$private$nn_holder$beta) == ncol(z))
  testthat::expect_true(length(cot$.__enclos_env__$private$nn_holder$beta) == ncol(x))
  
  # no bf, tensor
  testthat::expect_true(inherits(causalOT:::cotDualTrain$new(Measure(x, adapt = "weights"), Measure(z))$setup_arguments()$.__enclos_env__$private$nn_holder, "cotDualOpt"))
  
  #no bf, keops
  causalOT:::rkeops_check()
  testthat::expect_true(inherits(causalOT:::cotDualTrain$new(Measure(x, adapt = "weights"), Measure(z))$setup_arguments(cost.online = "online")$.__enclos_env__$private$nn_holder, "cotDualOpt_keops"))
  
  #no bf, keops
  testthat::expect_true(inherits(causalOT:::cotDualTrain$new(Measure(x, adapt = "weights", target.values = colMeans(z)), Measure(z))$setup_arguments(cost.online = "online")$.__enclos_env__$private$nn_holder, "cotDualBfOpt_keops"))
  
  
  #### test weights function ###
  nnh <- cot$.__enclos_env__$private$nn_holder
  priv <- cot$.__enclos_env__$private
  a1 <- nnh$calc_w1(nnh$gamma, priv$C_xy$data, priv$a_log,
              priv$b_log, torch::jit_scalar(priv$lambda),
              torch::jit_scalar(as.integer(n)))
  a2 <- nnh$calc_w2(nnh$gamma, priv$C_xx$data, priv$a_log,
                     torch::jit_scalar(priv$lambda),
                    torch::jit_scalar(as.integer(n)))
  # debugonce(cot$.__enclos_env__$.__active__$weights)
  w <- cot$weights
  # testthat::expect_true(length(w) == 3)
  # testthat::expect_equal(as.numeric(w[[2]]), as.numeric(a1))
  # testthat::expect_equal(as.numeric(w[[3]]), as.numeric(a2))
  # testthat::expect_equal(as.numeric(w[[1]]), as.numeric(a2 + a1)*0.5)
  testthat::expect_equal(as.numeric(w$to(device = "cpu")), as.numeric(((a2 + a1)*0.5)$to(device = "cpu")))
  
  testthat::expect_equal(names(cot$.__enclos_env__$private$parameters),
                         c("gamma", "beta"))
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$gamma$to(device = "cpu")),
                         as.numeric(cot$.__enclos_env__$private$parameters$gamma$to(device = "cpu")))
  testthat::expect_equal(rlang::obj_address(cot$.__enclos_env__$private$nn_holder$gamma),
                         rlang::obj_address(cot$.__enclos_env__$private$parameters$gamma))
  torch::with_no_grad(cot$.__enclos_env__$private$nn_holder$beta$copy_(c(1,2)))
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$beta$to(device = "cpu")),
                         c(1,2))
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$parameters$beta$to(device = "cpu")),
                         c(1,2))
  
  
  # test that set_lambda works
  testthat::expect_true(length(cot$penalty$lambda) > 1)
  priv <- cot$.__enclos_env__$private
  testthat::expect_true(priv$lambda == cot$penalty$lambda[1L])
  priv$set_lambda(4)
  testthat::expect_equal(priv$lambda , torch::jit_scalar(4))
  testthat::expect_error(priv$set_lambda(-1))
  
  # test that set_delta works
  testthat::expect_true(length(cot$penalty$delta) > 1)
  priv <- cot$.__enclos_env__$private
  testthat::expect_true(priv$delta == "numeric")
  priv$set_delta(.4)
  testthat::expect_equal(priv$delta , torch::jit_scalar(.4))
  testthat::expect_error(priv$set_lambda(-1))
  
  # test that set_penalties works
  priv <- cot$.__enclos_env__$private
  priv$set_penalties(c(lambda = Inf, delta = .4))
  testthat::expect_equal(priv$delta , torch::jit_scalar(.4))
  testthat::expect_equal(priv$lambda, torch::jit_scalar(359871.9312),
                         tol = 1e-5)
  
  testthat::expect_warning(priv$set_penalties(c(5,5)))
  testthat::expect_silent(priv$set_penalties(5))
  testthat::expect_error(priv$set_penalties(c(steve = 5,5)))
  
  priv$set_penalties(list(lambda = 50, delta = 5))
  testthat::expect_equal(priv$delta , torch::jit_scalar(5),
                         tol = 1e-5)
  testthat::expect_equal(priv$lambda, torch::jit_scalar(50),
                         tol = 1e-5)
  
  # make sure optimization setup works
  # debugonce(priv$torch_optim_setup)
  priv$torch_optim_setup(torch_optim = torch::optim_rmsprop,
                         torch_scheduler = torch::lr_multiplicative,
                         torch_args = NULL)
  testthat::expect_true(
    inherits(priv$opt, "optim_rmsprop")
  )
  testthat::expect_true(
    inherits(priv$sched, "lr_multiplicative")
  )
  # testthat::expect_equal(
  #   capture.output(print(priv$sched$lr_lambdas[[1]]))[1], 
  #   "function(epoch) {0.99}"
  # )
  
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$gamma$to(device = "cpu")),
                         as.numeric(cot$.__enclos_env__$private$parameters$gamma$params$to(device = "cpu")),
                         tol = 1e-5)
  
  testthat::expect_equal(1e-2, #priv$lambda/100,
                         cot$.__enclos_env__$private$parameters$gamma$lr)
  
  
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$beta$to(device = "cpu")),
                         as.numeric(cot$.__enclos_env__$private$parameters$beta$params$to(device = "cpu")))
  
  testthat::expect_equal(0.01,
                         cot$.__enclos_env__$private$parameters$beta$lr)
  
  # torch_optim_reset
  # debugonce(priv$torch_optim_reset)
  priv <- cot$.__enclos_env__$private
  old_add <- rlang::obj_address(priv$opt)
  priv$torch_optim_reset(0.44)
  testthat::expect_equal(0.44, #priv$lambda/100,
                         cot$.__enclos_env__$private$parameters$gamma$lr)
  testthat::expect_equal(0.44,
                         cot$.__enclos_env__$private$parameters$beta$lr)
  
  
  testthat::expect_true(rlang::obj_address(priv$opt) != old_add)
  testthat::expect_equal(rlang::obj_address(priv$nn_holder$gamma),
                        rlang::obj_address(priv$parameters$gamma$params))
  
  
  # optimization_loop
  # debugonce(priv$optimization_loop)
  out <- priv$optimization_loop(2, 1e-4)
  testthat::expect_true(out$iter == 2)
  
  testthat::expect_equal(rlang::obj_address(priv$nn_holder$gamma),
                         rlang::obj_address(priv$parameters$gamma$params))
  testthat::expect_equal(as.numeric(cot$.__enclos_env__$private$nn_holder$gamma$to(device = "cpu")),
                         as.numeric(cot$.__enclos_env__$private$parameters$gamma$params$to(device = "cpu")))
  testthat::expect_true(all(as.numeric(cot$.__enclos_env__$private$nn_holder$gamma$to(device = "cpu")) != 0) )
  
  # test parameters get set
  pars <- priv$parameters
  testthat::expect_true(pars$gamma$params$requires_grad == TRUE)
  testthat::expect_equal(as.numeric(pars$gamma$params$to(device = "cpu")), as.numeric(cot$.__enclos_env__$private$nn_holder$gamma$to(device = "cpu")))
  
  
  pars <- priv$parameters_get_set()
  ws    <- pars[[ls(pars)]]
  w2    <- cot$weights
  # testthat::expect_equal(as.numeric(ws[[1]]), as.numeric(w2[[1]]))
  testthat::expect_equal(as.numeric(ws$to(device = "cpu")), as.numeric(w2$to(device = "cpu")))
  
  ms <- cot$.__enclos_env__$private$measures
  m  <- NULL
  for (i in ls(ms)) {
    if(ms[[i]]$adapt == "weights") {
      m <- ms[[i]]
      break
    }
  }
  
  testthat::expect_error(priv$parameters_get_set(ws ))
  testthat::expect_error(priv$parameters_get_set(list(ws,ws) ))
  testthat::expect_silent(priv$parameters_get_set(list(ws) ))
  testthat::expect_equal(as.numeric(m$weights$to(device = "cpu")),
                         as.numeric(ws$to(device = "cpu")), tol = 1e-5)
  
  # testthat::expect_true(rlang::obj_address(cot$.__enclos_env__$private$nn_holder$gamma) == rlang::obj_address(pars$gamma))
  # 
  # pars <- priv$parameters_get_set(clone = TRUE)
  # testthat::expect_true(pars$gamma$requires_grad == FALSE)
  # testthat::expect_equal(as.numeric(pars$gamma), as.numeric(cot$.__enclos_env__$private$nn_holder$gamma))
  # testthat::expect_true(rlang::obj_address(cot$.__enclos_env__$private$nn_holder$gamma) != rlang::obj_address(pars$gamma))
  # 
  # pars$gamma <- pars$gamma * 0 + 1
  # priv$parameters_get_set(pars)
  # testthat::expect_equal(as.numeric(pars$gamma), as.numeric(priv$nn_holder$gamma))
  # testthat::expect_true(rlang::obj_address(cot$.__enclos_env__$private$nn_holder$gamma) != rlang::obj_address(pars$gamma))
  # testthat::expect_true(inherits(pars, "weightEnv"))
  
  
  # hyperparam
  cot <- causalOT:::cotDualTrain$new(m1,mt)
  cot$setup_arguments()
  # debugonce(cot$solve)
  cot$solve(niter = 1L, torch_optim = torch::optim_rmsprop, torch_scheduler = torch::lr_multiplicative)
  # debugonce(private$parameters_get_set)
  # debugonce(private$iterate_over_delta)
  # f
  # debugonce(cot$choose_hyperparameters)
  # cot$choose_hyperparameters()
  # debugonce(private$setup_choose_hyperparameters)
  
  
  
  testthat::expect_silent( cot$choose_hyperparameters(n_boot_lambda = 10, n_boot_delta = 10) )
  testthat::expect_true(is.numeric(cot$selected_delta[[1]]))
  testthat::expect_true(cot$selected_lambda < 359871.93116805560749)
})
