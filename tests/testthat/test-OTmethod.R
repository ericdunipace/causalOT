testthat::test_that("OT object forms tensor", {
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  # giving masses
  testthat::expect_silent(ot1 <- causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
  cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
  diameter=NULL))
  
  # no masses
  testthat::expect_silent(ot2 <- causalOT:::OT$new(x = x, y = y, penalty = penalty, 
                                       cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                                       diameter=NULL))
  
  testthat::expect_equal(ot1, ot2)
  
})

testthat::test_that("OT object forms online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  # giving masses
  causalOT:::rkeops_check()
  ot1 <- causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                                        cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                                        diameter=NULL)
  
  # no masses
  testthat::expect_silent(ot2 <- causalOT:::OT$new(x = x, y = y, penalty = penalty, 
                                        cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                                        diameter=NULL))
  
  testthat::expect_equal(ot1, ot2)
  testthat::expect_equal(ot1$diameter,  21.12901, tolerance = 1e-4)
  
})

testthat::test_that("sinkhorn_loop runs, tensor", {
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  output <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  loss <- sum(output$f_xy * at) + 
    sum(output$g_yx * bt)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(), 2.786224, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 0, tolerance = 1e-3)
  # compare
  # pot <- causalOT::sinkhorn(x = x*sqrt(0.5), y = y*sqrt(0.5), a = a, b = b, power = 2, blur = 10, scaling = 0.99, debias = FALSE )
  # sum(pot$f * a) + sum(pot$g * b)
  #2.786224 total, f = 1.390419, g = 1.395806
})

testthat::test_that("sinkhorn_loop runs, online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  niter <- 1000
  tol <- 1e-8
  
  # giving masses
  causalOT:::rkeops_check()
  ot <- causalOT:::OT$new(x = x, y = y, 
                          a = a, b = b, 
                          penalty = penalty,
                          cost = NULL, p = 2, 
                          debias = TRUE, tensorized = "online",
                          diameter=NULL)
  output <- ot$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  output1 <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  
  loss <- sum(output$f_xy * torch::torch_tensor(a,
                                        dtype = output$f_xy$dtype,
                                        device = output$f_xy$device)) + 
    sum(output$g_yx * torch::torch_tensor(b,
                                          dtype = output$g_yx$dtype,
                                          device = output$g_yx$device))
  loss1 <- sum(output1$f_xy * 
                 torch::torch_tensor(a,
                                     dtype = output1$f_xy$dtype,
                                     device = output1$f_xy$device) ) + 
    sum(output1$g_yx * 
       torch::torch_tensor(b,
                           dtype = output1$g_xy$dtype,
                           device = output1$f_xy$device))
  testthat::expect_equal(loss$item(), loss1$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * torch::torch_tensor(a,
                                                               dtype = output$f_xy$dtype,
                                                               device = output$f_xy$device) )$item(),
                         sum(output1$f_xy * torch::torch_tensor(a,
                                                                dtype = output$f_xy$dtype,
                                                                device = output$f_xy$device) )$item(), tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * torch::torch_tensor(b,
                                                               dtype = output$g_yx$dtype,
                                                               device = output$g_yx$device))$item(), 
                         sum(output1$g_yx * torch::torch_tensor(b,
                                                                dtype = output1$g_xy$dtype,
                                                                device = output1$f_xy$device))$item(), tolerance = 1e-3)
  
})

testthat::test_that("sinkhorn_self runs, tensor", {
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output_x <- ot1$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y <- ot1$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  loss_x <- sum(output_x * torch::torch_tensor(a,
                                               dtype = output_x$dtype,
                                               device = output_x$device) ) * 2
  testthat::expect_equal(loss_x$item(), 2.17266, tolerance = 1e-4)
  
  loss_y <- sum(output_y * torch::torch_tensor(b,
                                               dtype = output_y$dtype,
                                               device = output_y$device) ) * 2
  testthat::expect_equal(loss_y$item(), 2.99741, tolerance = 1e-4)
  
})

testthat::test_that("sinkhorn_self runs, online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  niter <- 1000
  tol <- 1e-8
  
  # giving masses
  causalOT:::rkeops_check()
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = TRUE, tensorized = "online",
               diameter=NULL)
  output_x <-  ot$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y <-  ot$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output_x1 <-  ot1$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y1 <-  ot1$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  loss_x <- sum(output_x * torch::torch_tensor(a,
                                                dtype = output_x$dtype,
                                                device = output_x$device) ) * 2
  loss_x1 <- sum(output_x1 * torch::torch_tensor(a,
                                                  dtype = output_x1$dtype,
                                                  device = output_x1$device) ) * 2
  testthat::expect_equal(loss_x$item(), loss_x1$item(), tolerance = 1e-4)
  
  loss_y <- sum(output_y * torch::torch_tensor(b,
                                               dtype = output_y$dtype,
                                               device = output_y$device) ) * 2
  loss_y1 <- sum(output_y1 * torch::torch_tensor(b,
                                                 dtype = output_y1$dtype,
                                                 device = output_y1$device) ) * 2
  testthat::expect_equal(loss_y$item(), loss_y1$item(), tolerance = 1e-4)
  
})

testthat::test_that("sinkhorn_cot runs, tensor", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  ot$sinkhorn_cot(niter, tol)
  output <- ot$potentials
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  
  loss <- sum(output$f_xy *  at) + sum(output$g_yx * bt)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(),2.786224, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 0, tolerance = 1e-3)
  testthat::expect_true(all(output$f_xx==0))
  testthat::expect_true(all(as.logical(output$g_yy==0)))
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                    cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                    diameter=NULL)
  
  output1 <- ot1$sinkhorn_cot(niter, tol)$potentials
  loss1 <- sum(output$f_xy * at ) + sum(output$g_yx * bt)
  
  testthat::expect_equal(loss1$item(), loss$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output1$f_xx * at )$item()*2, 2.17266, tolerance = 1e-4)
  testthat::expect_true(all(as.logical((output1$g_yy==0)$to(device = "cpu"))))
})

testthat::test_that("sinkhorn_cot runs, online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  causalOT:::rkeops_check()
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = FALSE, tensorized = "online",
               diameter=NULL)
  
  output <- ot$sinkhorn_cot(niter, tol)$potentials
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  
  loss <- sum(output$f_xy *  at) + sum(output$g_yx * bt)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(),2.786224, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 0, tolerance = 1e-3)
  testthat::expect_true(all(as.logical(output$f_xx==0)))
  testthat::expect_true(all(as.logical(output$g_yy==0)))
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                diameter=NULL)
  
  output1 <- ot1$sinkhorn_cot( niter, tol)$potentials
  loss1 <- sum(output$f_xy * at ) + sum(output$g_yx * bt)
  
  testthat::expect_equal(loss1, loss, tolerance = 1e-4)
  testthat::expect_equal(sum(output1$f_xx * at )$item()*2, 2.17266, tolerance = 1e-4)
  testthat::expect_true(all(as.logical((output1$g_yy==0)$to(device = "cpu"))))
  
})

testthat::test_that("sinkhorn_dist runs, tensor", {
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  ot1 <- causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  testthat::expect_error( causalOT:::sinkhorn_dist(ot1))
  ot1$sinkhorn_opt(niter, tol)
  loss <- causalOT:::sinkhorn_dist(ot1)
  output <- ot1$potentials
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(),2.786224, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 0, tolerance = 1e-3)
  testthat::expect_equal(loss$item(), sum(output$f_xy * at )$item() +
                           sum(output$g_yx * bt)$item(), tol = 1e-4)
  # check primal
  primal_loss <- (exp((output$f_xy$view(c(n,1)) - ot1$C_xy$data + output$g_yx$view(c(1,ot1$m)))/ot1$penalty + ot1$.__enclos_env__$private$a_log$view(c(n,1)) + ot1$.__enclos_env__$private$b_log$view(c(1,m))) * ot1$C_xy$data)$sum()
  testthat::expect_true((primal_loss < loss)$item())
  
  # compare
  # pot <- causalOT::sinkhorn(x = x*sqrt(0.5), y = y*sqrt(0.5), a = a, b = b, power = 2, blur = 10, scaling = 0.99, debias = FALSE )
  # sum(pot$f * a) + sum(pot$g * b)
  #2.786224 total, f = 1.390419, g = 1.395806
})

testthat::test_that("sinkhorn_dist runs, online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  niter <- 1000
  tol <- 1e-8
  
  # giving masses
  causalOT:::rkeops_check()
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = TRUE, tensorized = "online",
               diameter=NULL)
  output <- ot$sinkhorn_opt(niter, tol)$potentials
  loss <- causalOT:::sinkhorn_dist(ot)
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output1 <- ot1$sinkhorn_opt(niter, tol)$potentials
  loss1 <- causalOT:::sinkhorn_dist(ot1)
  
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  
  testthat::expect_equal(loss$item(), loss1$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(),
                         sum(output1$f_xy * at )$item(), tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 
                         sum(output1$g_yx * bt)$item(), tolerance = 1e-3)
  
})

testthat::test_that("sinkhorn_loop runs, tensor", {
  causalOT:::torch_check()
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  niter <- 1000
  tol <- 1e-10
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                          cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                          diameter=NULL)
  
  output <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  at <- torch::torch_tensor(a,
                            dtype = output$f_xy$dtype,
                            device = output$f_xy$device)
  bt <- torch::torch_tensor(b,
                            dtype = output$g_yx$dtype,
                            device = output$g_yx$device)
  loss <- sum(output$f_xy * at) + 
    sum(output$g_yx * bt)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * at )$item(), 2.786224, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * bt)$item(), 0, tolerance = 1e-3)
  # compare
  # pot <- causalOT::sinkhorn(x = x*sqrt(0.5), y = y*sqrt(0.5), a = a, b = b, power = 2, blur = 10, scaling = 0.99, debias = FALSE )
  # sum(pot$f * a) + sum(pot$g * b)
  #2.786224 total, f = 1.390419, g = 1.395806
})

testthat::test_that("sinkhorn_loop gradient", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- torch::torch_tensor(matrix(stats::rnorm(n*d), n, d), 
                           dtype = torch::torch_double(),
                           requires_grad = TRUE)
  y <- torch::torch_tensor(matrix(stats::rnorm(m*d), m, d), 
                           dtype = torch::torch_double(),
                           requires_grad = TRUE)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  niter <- 1000
  tol <- 1e-8
  
  # giving masses
  causalOT:::rkeops_check()
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                         cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                         diameter=NULL)
  output <- ot$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  testthat::expect_true(output$f_xy$requires_grad)
  testthat::expect_true(output$g_yx$requires_grad)
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                          cost = NULL, p = 2, 
                          debias = FALSE, tensorized = "tensorized",
                          diameter=NULL)
  
  output1 <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  testthat::expect_true(output1$f_xy$requires_grad)
  testthat::expect_true(output1$g_yx$requires_grad)
  
})

testthat::test_that("sinkhorn_loop gradient", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- 10
  x <- torch::torch_tensor(matrix(stats::rnorm(n*d), n, d), 
                           dtype = torch::torch_double(),
                           requires_grad = TRUE)
  y <- torch::torch_tensor(matrix(stats::rnorm(m*d), m, d), 
                           dtype = torch::torch_double(),
                           requires_grad = TRUE)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  niter <- 1000
  tol <- 1e-8
  
  ot1 <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                          cost = NULL, p = 2, 
                          debias = TRUE, tensorized = "tensorized",
                          diameter=NULL)
  
  output1x <- ot1$.__enclos_env__$private$sinkhorn_self_loop(which.margin = "x", niter, tol)
  output1y <- ot1$.__enclos_env__$private$sinkhorn_self_loop(which.margin = "y", niter, tol)
  testthat::expect_true(output1x$requires_grad)
  testthat::expect_true(output1y$requires_grad)
  
  # giving masses
  causalOT:::rkeops_check()
  ot <-causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                         cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                         diameter=NULL)
  outputx <- ot$.__enclos_env__$private$sinkhorn_self_loop(which.margin = "x", niter, tol)
  outputy <- ot$.__enclos_env__$private$sinkhorn_self_loop(which.margin = "y", niter, tol)
  testthat::expect_true(outputx$requires_grad)
  testthat::expect_true(outputy$requires_grad)
  
  
  
  
})

testthat::test_that("OT infinite penalty distances online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  causalOT:::torch_check()
  
  set.seed(1231)
  n <- 15
  m <- 13
  d <- 3
  penalty <- Inf
  x <- matrix(stats::rnorm(n*d), n, d)
  y <- matrix(stats::rnorm(m*d), m, d)
  a <- rep(1/n, n)
  b <- rep(1/m, m)
  
  # giving masses
  causalOT:::rkeops_check()
  ot1 <- causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                           cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                           diameter=NULL)
  
  testthat::expect_silent(loss1 <- causalOT:::energy_dist(ot1))
  testthat::expect_silent(loss2 <- causalOT:::inf_sinkhorn_dist(ot1))
  testthat::expect_equal(loss1,loss2)
  
  causalOT:::rkeops_check()
  ot2 <- causalOT:::OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                           cost = NULL, p = 2, debias = FALSE, tensorized = "online",
                           diameter=NULL)
  
  testthat::expect_silent(loss3 <- causalOT:::inf_sinkhorn_dist(ot2))
  testthat::expect_true(ot2$a$dtype == loss3$dtype)
  
  # debugonce(causalOT:::inf_sinkhorn_dist)
  ot2$a <- torch::torch_tensor(ot2$a, requires_grad = TRUE)
  loss <- causalOT:::inf_sinkhorn_dist(ot2)
  testthat::expect_silent(loss$backward())
})