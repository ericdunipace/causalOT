testthat::test_that("OT object forms tensor", {
  
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
  testthat::expect_silent(ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
  cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
  diameter=NULL))
  
  # no masses
  testthat::expect_silent(ot2 <- OT$new(x = x, y = y, penalty = penalty, 
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
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  output <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  
  loss <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),1.390419, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 1.395806, tolerance = 1e-3)
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
  ot <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                                         cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                                         diameter=NULL)
  output <- ot$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  output1 <- ot1$.__enclos_env__$private$sinkhorn_loop(niter, tol)
  
  loss <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  loss1 <- sum(output1$f_xy * a ) + sum(output1$g_yx * b)
  testthat::expect_equal(loss$item(), loss1$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),
                         sum(output1$f_xy * a )$item(), tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 
                         sum(output1$g_yx * b)$item(), tolerance = 1e-3)
  
})

testthat::test_that("sinkhorn_self runs, tensor", {
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
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output_x <- ot1$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y <- ot1$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  loss_x <- sum(output_x * a ) * 2
  testthat::expect_equal(loss_x$item(), 2.17266, tolerance = 1e-4)
  
  loss_y <- sum(output_y * b ) * 2
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
  ot <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = TRUE, tensorized = "online",
               diameter=NULL)
  output_x <-  ot$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y <-  ot$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output_x1 <-  ot1$.__enclos_env__$private$sinkhorn_self_loop("x", niter, tol)
  output_y1 <-  ot1$.__enclos_env__$private$sinkhorn_self_loop("y", niter, tol)
  
  loss_x <- sum(output_x * a ) * 2
  loss_x1 <- sum(output_x1 * a ) * 2
  testthat::expect_equal(loss_x$item(), loss_x1$item(), tolerance = 1e-4)
  
  loss_y <- sum(output_y * b ) * 2
  loss_y1 <- sum(output_y1 * b ) * 2
  testthat::expect_equal(loss_y$item(), loss_y1$item(), tolerance = 1e-4)
  
})

testthat::test_that("sinkhorn_cot runs, tensor", {
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
  
  ot <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = FALSE, tensorized = "auto",
                diameter=NULL)
  
  ot$sinkhorn_cot(niter, tol)
  output <- ot$potentials
  
  loss <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),1.390419, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 1.395806, tolerance = 1e-3)
  testthat::expect_true(all(output$f_xx==0))
  testthat::expect_true(all(as.logical(output$g_yy==0)))
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                    cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                    diameter=NULL)
  
  output1 <- ot1$sinkhorn_cot(niter, tol)$potentials
  loss1 <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  
  testthat::expect_equal(loss1$item(), loss$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output1$f_xx * a )$item()*2, 2.17266, tolerance = 1e-4)
  testthat::expect_true(all(as.logical(output1$g_yy==0)))
})

testthat::test_that("sinkhorn_cot runs, online", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed(pkg="rkeops")
  testthat::skip_on_ci()
  
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
  
  ot <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = FALSE, tensorized = "online",
               diameter=NULL)
  
  output <- ot$sinkhorn_cot(niter, tol)$potentials
  
  loss <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),1.390419, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 1.395806, tolerance = 1e-3)
  testthat::expect_true(all(as.logical(output$f_xx==0)))
  testthat::expect_true(all(as.logical(output$g_yy==0)))
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "online",
                diameter=NULL)
  
  output1 <- ot1$sinkhorn_cot( niter, tol)$potentials
  loss1 <- sum(output$f_xy * a ) + sum(output$g_yx * b)
  
  testthat::expect_equal(loss1, loss, tolerance = 1e-4)
  testthat::expect_equal(sum(output1$f_xx * a )$item()*2, 2.17266, tolerance = 1e-4)
  testthat::expect_true(all(as.logical(output1$g_yy==0)))
  
})

testthat::test_that("sinkhorn_dist runs, tensor", {
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
  
  testthat::expect_error( sinkhorn_dist(ot1))
  ot1$sinkhorn_opt(niter, tol)
  loss <- causalOT:::sinkhorn_dist(ot1)
  output <- ot1$potentials
  
  testthat::expect_equal(loss$item(), 2.786224, tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),1.390419, tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 1.395806, tolerance = 1e-3)
  testthat::expect_equal(loss$item(), sum(output$f_xy * a )$item() +
                           sum(output$g_yx * b)$item(), tol = 1e-4)
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
  ot <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
               cost = NULL, p = 2, debias = TRUE, tensorized = "online",
               diameter=NULL)
  output <- ot$sinkhorn_opt(niter, tol)$potentials
  loss <- sinkhorn_dist(ot)
  
  ot1 <- OT$new(x = x, y = y, a = a, b = b, penalty = penalty, 
                cost = NULL, p = 2, debias = TRUE, tensorized = "auto",
                diameter=NULL)
  
  output1 <- ot1$sinkhorn_opt(niter, tol)$potentials
  loss1 <- sinkhorn_dist(ot1)
  
  testthat::expect_equal(loss$item(), loss1$item(), tolerance = 1e-4)
  testthat::expect_equal(sum(output$f_xy * a )$item(),
                         sum(output1$f_xy * a )$item(), tolerance = 1e-3)
  testthat::expect_equal(sum(output$g_yx * b)$item(), 
                         sum(output1$g_yx * b)$item(), tolerance = 1e-3)
  
})