testthat::test_that("measure forms", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  x <- matrix(rnorm(100,10), 100, 10)
  m <- Measure(x = x)
  
  testthat::expect_error(m$weights <- 1)
  testthat::expect_error(m$weights <- NA)
  testthat::expect_error(m$weights <- NULL)
  testthat::expect_true(m$probability_measure)
  testthat::expect_true(m$adapt == "none")
  testthat::expect_true( is.na(m$balance_functions) )
  testthat::expect_true( is.na(m$balance_target) )
  testthat::expect_true(isFALSE(m$weights$requires_grad))
  
  mc <- m$clone(deep = TRUE)
  testthat::expect_true(rlang::obj_address(mc) != rlang::obj_address(m))
  testthat::expect_true(rlang::obj_address(m$x) != rlang::obj_address(mc$x))
  testthat::expect_equal(as_matrix(m$x), as_matrix(mc$x))
  testthat::expect_equal(as_matrix(m$weights), as_matrix(mc$weights))
  
  m_w <- Measure(x = x, adapt = "weights", dtype = torch::torch_double())
  testthat::expect_true(isTRUE(as_logical(m_w$weights$requires_grad)))
  testthat::expect_true( is.na(m_w$balance_functions) )
  testthat::expect_true( is.na(m_w$balance_target) )
  testthat::expect_silent(m_w$weights <- rep(1.0/100,100))
  first_weights <- m_w$weights$clone()
  error_weights <- rep(1.0,100)
  m_w$weights <- error_weights
  testthat::expect_equal(first_weights, m_w$weights)
  error_weights[1] <- -1
  testthat::expect_error(m_w$weights <- error_weights)
  loss <- sum(m_w$weights*5)
  loss$backward()
  testthat::expect_true(all(as_logical(m_w$.__enclos_env__$private$mass_$grad$to(device = "cpu") >0)))
  mcw <- m_w$clone(deep = TRUE)
  testthat::expect_true(rlang::obj_address(mcw) != rlang::obj_address(m_w))
  testthat::expect_true(rlang::obj_address(m_w$x) != rlang::obj_address(mcw$x))
  testthat::expect_equal(as_matrix(m_w$x), as_matrix(mcw$x))
  testthat::expect_equal(as_matrix(m_w$weights), as_matrix(mcw$weights))
  testthat::expect_true(m_w$requires_grad)
  testthat::expect_true(mcw$requires_grad)
  
  mcdw <- m_w$detach()
  testthat::expect_true(rlang::obj_address(mcdw) != rlang::obj_address(m_w))
  testthat::expect_true(rlang::obj_address(m_w$x) != rlang::obj_address(mcdw$x))
  testthat::expect_equal(as_matrix(m_w$x), as_matrix(mcdw$x))
  testthat::expect_equal(as_matrix(m_w$weights), as_matrix(mcdw$weights))
  testthat::expect_true(m_w$requires_grad)
  testthat::expect_true(!mcdw$requires_grad)
  
  m_x <- Measure(x = x, adapt = "x")
  testthat::expect_true(isFALSE(as_logical(m_x$weights$requires_grad)))
  testthat::expect_true(isTRUE(as_logical(m_x$x$requires_grad)))
  testthat::expect_true( is.na(m_x$balance_functions) )
  testthat::expect_true( is.na(m_x$balance_target) )
  testthat::expect_equal(as_numeric(m_x$weights), rep(1.0/100,100))
  testthat::expect_silent(m_x$weights <- rep(1.0/100,100))
  first_weights <- m_x$weights$clone()
  error_weights <- rep(1.0,100)
  m_x$weights <- error_weights
  testthat::expect_equal(first_weights, m_x$weights)
  error_weights[1] <- -1
  testthat::expect_error(m_x$weights <- error_weights)
  loss <- sum(m_x$x*5)
  loss$backward()
  testthat::expect_true(all(as_logical(m_x$x$grad == 5)))
  
  
  # checking the balance targets stuff
  y <- matrix(rnorm(1024*10), 1024,10) + 2
  m_t <- Measure(x = x, balance.functions = x, target.values = colMeans(y))
  testthat::expect_true(isFALSE(as_logical(m_t$weights$requires_grad)))
  testthat::expect_true(isFALSE(as_logical(m_t$x$requires_grad)))
  testthat::expect_true( is.na(m_t$balance_functions) )
  testthat::expect_true( is.na(m_t$balance_target) )
  
  testthat::expect_warning(m_wt <- Measure(x = x, adapt = "weights", balance.functions = matrix(0, 100,10), target.values = colMeans(y)))
  testthat::expect_true(is.na(m_wt$balance_functions))
  m_wt <- Measure(x = x, adapt = "weights", balance.functions = x, target.values = colMeans(y))
  testthat::expect_true(isTRUE(as_logical(m_wt$weights$requires_grad)))
  testthat::expect_true(isFALSE(as_logical(m_wt$x$requires_grad)))
  testthat::expect_equal(ncol(m_wt$balance_functions), ncol(x))
  testthat::expect_true( !all(as_logical(m_wt$balance_functions$isnan()) ))
  testthat::expect_true( !all(is.na(m_wt$balance_target) ))
  
  testthat::expect_warning(m_xt <- Measure(x = x, adapt = "x", balance.functions = x, target.values = colMeans(y)))
  testthat::expect_true(isFALSE(as_logical(m_xt$weights$requires_grad)))
  testthat::expect_true(isTRUE(as_logical(m_xt$x$requires_grad)))
  testthat::expect_equal(ncol(m_xt$balance_functions), ncol(x))
  testthat::expect_true( !all(as_logical(m_xt$balance_functions$isnan()) ))
  testthat::expect_true( !all(is.na(m_xt$balance_target) ))
  testthat::expect_true(rlang::obj_address(m_xt$balance_functions)!= rlang::obj_address(m_xt$x))
  
  testthat::expect_silent(m_xt <- Measure(x = x, adapt = "x",  target.values = colMeans(y)))
  testthat::expect_true(isFALSE(as_logical(m_xt$weights$requires_grad)))
  testthat::expect_true(isTRUE(as_logical(m_xt$x$requires_grad)))
  testthat::expect_equal(ncol(m_xt$balance_functions), ncol(x))
  testthat::expect_true( !all(as_logical(m_xt$balance_functions$isnan()) ))
  testthat::expect_true( !all(is.na(m_xt$balance_target) ))
  testthat::expect_true(rlang::obj_address(m_xt$balance_functions)== rlang::obj_address(m_xt$x))
  
  
  m_1 <- Measure(x[,1], target.values = colMeans(y)[1], adapt = "weights")
  testthat::expect_equal(dim(m_1$balance_functions), c(100L, 1L))
  
})

testthat::test_that("OTProblem tests",{
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  x  <- matrix(rnorm(128*2) + 5, 128, 2)
  m1 <- Measure(x = x)
  
  y  <- matrix(rnorm(256*2), 256, 2) + matrix(c(0,2), 256,2, byrow = TRUE)
  m2 <- Measure(x = y)
  
  addresses <- c(rlang::obj_address(m1), rlang::obj_address(m2))
  
  ot <- OTProblem(m1, m2)
  
  m_wrong_dtype <- Measure(x = y, dtype = torch::torch_float16())
  testthat::expect_error(
    OTProblem(m1, m_wrong_dtype), label = "wrong type data storage throws error"
  )
  testthat::expect_error(
    OTProblem(x, m2), label = "Detect non Measure object 1"
  )
  testthat::expect_error(
    OTProblem(m1, y), label = "Detect non Measure object 2"
  )
  
  m_wrong_ncol <- Measure(x = matrix(0, 100, 20))
  testthat::expect_error(
    OTProblem(m1, m_wrong_ncol), label = "Detect wrong number of columns"
  )
  
  testthat::expect_silent(
    OTProblem(m1, m2, y)
  )
  
  testthat::expect_equal(
    sort(ls(ot$.__enclos_env__$private$measures)), 
    sort(addresses),
    label = "Check object addresses same original object"
  )
  
  testthat::expect_equal(
    sort(ls(ot$.__enclos_env__$private$problems)),
    sort(c(paste0(addresses[1], ", ", addresses[2]))),
    label = "Check problem addresses same as original object"
  )
  
  
  # check can add
  ot_x <- OTProblem(m1, m1)
  ot_y <- OTProblem(m2, m2)
  
  mult_check <- ot_x * 0.5
  
  testthat::expect_true(rlang::obj_address(mult_check) != rlang::obj_address(ot_x), label = "Check that multiplication creates new object")
  
  testthat::expect_equal(
    ls(mult_check$.__enclos_env__$private$measures),
    ls(ot_x$.__enclos_env__$private$measures)
  )
  testthat::expect_equal(
    ls(mult_check$.__enclos_env__$private$measures),
    rlang::obj_address(m1)
  )
  
  # debugonce(causalOT:::binaryop.OTProblem)
  ot_final <- ot - mult_check
  orig_final_address <- rlang::obj_address(ot_final)
  ot_final <- ot_final - 0.5 * ot_y
  testthat::expect_true(rlang::obj_address(ot_final) != orig_final_address)
  
  testthat::expect_equal(
    sort(ls(ot$.__enclos_env__$private$measures)), 
    sort(addresses),
    label = "Check object addresses same after final object"
  )
  
  testthat::expect_equal(
    sort(ls(ot$.__enclos_env__$private$problems)),
    sort(c(paste0(addresses[1], ", ", addresses[2]))),
    label = "Check problem addresses same as original object")
  
  testthat::expect_equal(
    sort(ls(ot_final$.__enclos_env__$private$problems)),
    sort(c(paste0(addresses[1], ", ", addresses[2]), 
           paste0(addresses[1], ", ", addresses[1]),
           paste0(addresses[2], ", ", addresses[2]))),
    label = "Check problem addresses appropriate for new object"
  )
  
  testthat::expect_equal(
    sort(ls(ot_final$.__enclos_env__$private$measures)),
    sort(ls(ot$.__enclos_env__$private$measures)),
    label = "Check final measures same as original"
  )

  # barycenter check
  m3 <- Measure(x = matrix(runif(64*2), 64, 2), adapt = "x")
  
  ot_1 <- OTProblem(m1, m3)
  ot_2 <- OTProblem(m2, m3)
  
  ot_bary <- 0.5 * ot_1 + 0.5 * ot_2
  
  testthat::expect_error(ot_bary$solve(niter = 1, tol = 1e-4),
                         label = "make sure the args are set")
  
  # debugonce(ot_bary$setup_arguments)
  ot_bary$setup_arguments()
  # debugonce(causalOT:::inf_sinkhorn_dist)
  
  testthat::expect_silent(
    ot_bary$solve(niter = 1, tol = 1e-4, torch_args = list(line_search_fn = "strong_wolfe"))
  )
  testthat::expect_warning(
    ot_bary$solve(niter = 1, tol = 1e-4)
  )
  
  testthat::expect_warning(ot_bary$setup_arguments(lambda = 10, debias = TRUE))
  
  # debugonce(ot_bary$solve)
  # needs lr about 1e-1
  testthat::expect_silent(
    ot_bary$solve(niter = 1, tol = 1e-8, torch_opt = torch::optim_rmsprop, torch_args = list(lr=1e-1))
  )
})

testthat::test_that("weights adapatiation", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  z  <- matrix(rnorm(64*2), 64, 2) + matrix(c(0,.5), 64,2, byrow = TRUE)
  x  <- matrix(rnorm(128*2), 128, 2)
  y  <- matrix(rnorm(256*2), 256, 2) + matrix(c(.5,1), 256,2, byrow = TRUE)
  
  mt <- Measure(x = z)
  m1 <- Measure(x = x, adapt = "weights")
  m2 <- Measure(x = y, adapt = "weights")
  ot <- OTProblem(m1, m2)
  
  # debugonce(ot$setup_arguments)
  ot$setup_arguments(debias = TRUE)
  
  # debugonce(ot$solve)
  # Rprof(tmp<-tempfile())
  # ot$solve(niter = 100L, tol = 1e-7,
  #          torch_args = list(lr = 1,
  #                            line_search_fn = "strong_wolfe",
  #                            history_size = 5))
  # Rprof(NULL)
  # summaryRprof(tmp)
  # unlink(tmp)
  
  testthat::expect_silent(ot$solve(niter = 1L, tol = 1e-7,
                                    torch_opt = torch::optim_rmsprop,
                                   torch_args = list(lr = 1e-3)))
  
  # debugonce(ot$choose_hyperparameters)
  ot$choose_hyperparameters()
  info <- ot$info()
  testthat::expect_named(info)
  testthat::expect_true(all(names(info) %in% c("loss", "hyperparam.metrics", 
                                               "iterations", "balance.function.differences")))
  testthat::expect_true(all(info$iterations==1))
  
  #### with targets ####
  # torch optim
  z  <- matrix(rnorm(64*2), 64, 2) + matrix(c(0,.5), 64,2, byrow = TRUE)
  x  <- matrix(rnorm(128*2), 128, 2)
  y  <- matrix(rnorm(256*2), 256, 2) + matrix(c(.5,1), 256,2, byrow = TRUE)
  
  mt <- Measure(x = z)
  m1 <- Measure(x = x, target.values = colMeans(z), adapt = "weights")
  m2 <- Measure(x = y,  target.values = colMeans(z), adapt = "weights")
  ot <- OTProblem(m1, m2)
  
  # debugonce(ot$setup_arguments)
  ot$setup_arguments(debias = TRUE)
  
  old_delta <- ot$penalty$delta
  
  ot$.__enclos_env__$private$set_delta(5)
  adds <- ls(ot$.__enclos_env__$private$target_objects)
  testthat::expect_equal(ot$.__enclos_env__$private$target_objects[[adds[1]]]$delta, 5)
  testthat::expect_equal(ot$.__enclos_env__$private$target_objects[[adds[2]]]$delta, 5)
  
  # debugonce(ot$.__enclos_env__$private$balance_function_check)
  ot$.__enclos_env__$private$delta_values_setup(run.quick=FALSE, osqp_args = NULL)
  testthat::expect_equal(ot$penalty$delta, old_delta)
  
  # debugonce(ot$.__enclos_env__$private$balance_function_check)
  ot$.__enclos_env__$private$delta_values_setup(run.quick=TRUE, osqp_args = NULL)
  to_names <- ls(ot$.__enclos_env__$private$target_objects)
  testthat::expect_equal(ot$.__enclos_env__$private$target_objects[[to_names[[1L]]]]$delta, 1e-04)
  testthat::expect_equal(ot$.__enclos_env__$private$target_objects[[to_names[[2L]]]]$delta, 1e-04)
  testthat::expect_true(is.na(ot$.__enclos_env__$private$penalty_list$delta))
 
  # debugonce(ot$solve)
  # Rprof(tmp<-tempfile())
  # ot$solve(niter = 100L, tol = 1e-7,
  #                                   quick.balance.function = TRUE,
  #                                   # torch_optim = torch::optim_rmsprop,
  #                                   torch_args = list(lr = 1,
  #                                                     line_search_fn = "strong_wolfe",
  #                                                     history_size = 5))
  # Rprof(NULL)
  # summaryRprof(tmp)
  # unlink(tmp)
  testthat::expect_silent(ot$solve(niter = 1L, tol = 1e-7,
           quick.balance.function = TRUE,
           # torch_optim = torch::optim_rmsprop,
           torch_args = list(lr = 1,
                             line_search_fn = "strong_wolfe",
                             history_size = 5)))
  
  # debugonce(ot$choose_hyperparameters)
  ot$choose_hyperparameters()
  testthat::expect_true(is.numeric(ot$selected_lambda))
  info <- ot$info()
  testthat::expect_named(info)
  testthat::expect_true(all(names(info) %in% c("loss", "hyperparam.metrics", 
                                               "iterations", "balance.function.differences")))
  testthat::expect_true(all(info$iterations==1))
  
  testthat::expect_warning(ot$choose_hyperparameters(n_boot_lambda = 1L))
  
  mt <- Measure(x = z)
  m1 <- Measure(x = x, target.values = colMeans(z), adapt = "weights")
  m2 <- Measure(x = y,  target.values = colMeans(z), adapt = "weights")
  ot <- OTProblem(m1, m2)
  ot$setup_arguments()
  
  # debugonce(ot$solve)
  testthat::expect_warning(ot$solve(niter = 1L, tol = 1e-7,
           quick.balance.function = FALSE,
           torch_args = list(lr = 1e-5,
                             max_iter = 1L,
                             max_eval = 1L,
                             history_size = 5)))
  # debugonce(ot$choose_hyperparameters)
  ot$choose_hyperparameters(n_boot_lambda = 10L, n_boot_delta  = 10L)
  
  testthat::expect_true(is.numeric(ot$selected_lambda))
  info <- ot$info()
  testthat::expect_named(info)
  testthat::expect_true(all(names(info) %in% c("loss", "hyperparam.metrics", 
                                               "iterations", "balance.function.differences")))
  testthat::expect_true(all(info$iterations==1))
})

testthat::test_that("bary with muilt groups", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  #bary center + two groups
  z  <- matrix(runif(64*2), 64, 2) + matrix(c(0,.5), 64,2, byrow = TRUE)
  x1  <- matrix(rnorm(128*2)-0.01, 128, 2)
  y1  <- matrix(rnorm(256*2), 256, 2) + matrix(c(.25,0), 256,2, byrow = TRUE)
  x2  <- matrix(rnorm(128*2)+0.25, 128, 2)
  y2  <- matrix(rnorm(256*2), 256, 2) + matrix(c(.5,.25), 256,2, byrow = TRUE)
  
  mt <- Measure(x = z, adapt = "x")
  m1 <- Measure(x = x1, adapt = "weights")
  m2 <- Measure(x = y1, adapt = "weights")
  m3 <- Measure(x = x2, adapt = "weights")
  m4 <- Measure(x = y2, adapt = "weights")
  ot <- OTProblem(m1, mt) * 0.5 + OTProblem(m2, mt) * 0.5+ 
    OTProblem(m3,mt) * 0.5 + OTProblem(m4, mt) * 0.5
  ot$setup_arguments(debias = TRUE)
  # debugonce(ot$solve)
  # ot$solve(niter = 100L, tol = 1e-5, torch_args = list(line_search_fn = "strong_wolfe"))
  testthat::expect_silent(ot$solve(niter = 1L, tol = 1e-3, torch_optim = torch::optim_rmsprop, torch_args = list(lr = 1e-5)))
  # debugonce(ot$choose_hyperparameters)
  ot$choose_hyperparameters(n_boot_lambda = 2L)
  testthat::expect_true(is.numeric(ot$selected_lambda))
  
  # bary center also targeting a mean
  mt <- Measure(x = z, target.values = rep(0.25,2), adapt = "x")
  m1 <- Measure(x = x1, target.values = rep(0.25,2), adapt = "weights")
  m2 <- Measure(x = y1, target.values = rep(0.25,2), adapt = "weights")
  m3 <- Measure(x = x2, target.values = rep(0.25,2), adapt = "weights")
  m4 <- Measure(x = y2, target.values = rep(0.25,2), adapt = "weights")
  ot <- OTProblem(m1, mt) * 0.5 + OTProblem(m2, mt) * 0.5+ 
    OTProblem(m3,mt) * 0.5 + OTProblem(m4, mt) * 0.5
  
  ot$setup_arguments(debias = TRUE)
  ot$solve(niter = 1L, torch_optim = torch::optim_rmsprop, torch_args = list(lr = 1e-5))
  ot$choose_hyperparameters(n_boot_lambda = 10)
  testthat::expect_true(is.numeric(ot$selected_lambda))
  
  # frank-wolfe
  # without balance functions
  z  <- matrix(rnorm(64*2), 64, 2) + matrix(c(0,.5), 64,2, byrow = TRUE)
  x  <- matrix(rnorm(128*2), 128, 2)
  y  <- matrix(rnorm(256*2), 256, 2) + matrix(c(.5,1), 256,2, byrow = TRUE)
  
  mt <- Measure(x = z)
  m1 <- Measure(x = x, adapt = "weights")
  m2 <- Measure(x = y, adapt = "weights")
  
  ot <- OTProblem(m1, mt)+ OTProblem(m2,mt)
  
  # debugonce(ot$setup_arguments)
  ot$setup_arguments(debias = TRUE)
  
  # debugonce(ot$solve)
  # debugonce(ot$.__enclos_env__$private$frankwolfe_step)
  ot$solve(niter = 1L, optimizer = "frank-wolfe", osqp_args = list(verbose = FALSE))
  # debugonce(ot$choose_hyperparameters)
  ot$choose_hyperparameters(n_boot_lambda = 10L)
  testthat::expect_true(is.numeric(ot$selected_lambda))
  
  # with balance functions!
  mt <- Measure(x = z)
  m1 <- Measure(x = x, target.values = colMeans(z), adapt = "weights")
  m2 <- Measure(x = y,  target.values = colMeans(z), adapt = "weights")
  
  ot <- OTProblem(m1, mt)+ OTProblem(m2,mt)
  
  # debugonce(ot$setup_arguments)
  ot$setup_arguments(debias = TRUE)
  
  # debugonce(ot$solve)
  # debugonce(ot$.__enclos_env__$private$frankwolfe_step)
  ot$solve(niter = 1L, optimizer = "frank-wolfe", osqp_args = list(verbose = FALSE))
  ot$choose_hyperparameters(n_boot_lambda = 10L)
  testthat::expect_true(is.numeric(ot$selected_lambda))
  
  # with barycenter
  #dble opt expect error
  mt <- Measure(x = z * 3 + 1, target.values = colMeans(z), adapt = "x")
  m1 <- Measure(x = x, target.values = colMeans(z), adapt = "weights")
  m2 <- Measure(x = y,  target.values = colMeans(z), adapt = "weights")
  
  ot <- OTProblem(m1, mt)+ OTProblem(m2,mt)
  ot$setup_arguments(debias = TRUE)
  # debugonce(ot$solve)
  testthat::expect_error(ot$solve(niter = 1L, optimizer = "frank-wolfe", osqp_args = list(verbose = FALSE)))
  #opt one at at time should work
  
  m10 <- Measure(x = x, target.values = colMeans(z))
  m20 <- Measure(x = y,  target.values = colMeans(z))
  ot_z <- OTProblem(m1$detach(), mt)+ OTProblem(m2$detach(),mt)
  ot_z$setup_arguments(debias = TRUE)
  # debugonce(ot_z$.__enclos_env__$private$torch_optim_step)
  ot_z$solve(niter = 1L, torch_optim = torch::optim_rmsprop)
  
  mt0 <- mt$detach()
  testthat::expect_true(!mt0$requires_grad)
  testthat::expect_true(mt$requires_grad)
  
  ot <- OTProblem(m1, mt0)+ OTProblem(m2,mt0)
  ot$setup_arguments(debias = TRUE)
  # debugonce(ot$solve)
  ot$solve(niter = 1L, optimizer = "frank-wolfe", osqp_args = list(verbose = FALSE))
  ot$choose_hyperparameters(n_boot_lambda = 10L)
  testthat::expect_true(is.numeric(ot$selected_lambda))
  
})
