vcov_test_ot <- function(x, y, a = NULL, b = NULL, debias = TRUE,
                         dtype = torch::torch_double(), tol = 1e-7) {
  ot <- causalOT:::OT$new(
    x = x, y = y, a = a, b = b, penalty = 1, cost = NULL, p = 2,
    debias = debias, tensorized = "tensorized", diameter = NULL,
    dtype = dtype
  )
  ot$sinkhorn_opt(niter = 1e4, tol = tol)
  ot
}

vcov_test_data <- function() {
  set.seed(1)
  list(
    x = matrix(stats::rnorm(20), 10, 2),
    y = matrix(stats::rnorm(16) + 0.5, 8, 2)
  )
}

expect_valid_vcov <- function(V, dim) {
  testthat::expect_equal(dim(V), c(dim, dim))
  testthat::expect_equal(V, t(V), tolerance = 1e-12)
  testthat::expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
}

# Plug-in of the limiting covariance of the potentials in Theorem 2.2 of
# Gonzalez-Sanz, Loubes & Niles-Weed (2022), arXiv:2207.07427:
#   [I    A_Q] [df]        [Z_y]   Z_y(x) = (1/m) sum_j xi(x, Y_j) - E xi(x, Y)
#   [A_P  I  ] [dg] = -eps [Z_x]   Z_x(y) = (1/n) sum_i xi(X_i, y) - E xi(X, y)
# with g normalised so that Q(g) = 0 and Z_x, Z_y independent.
gln_vcov <- function(ot) {
  eps <- ot$penalty
  a <- as.numeric(ot$a)
  b <- as.numeric(ot$b)
  n <- length(a)
  m <- length(b)
  P <- as.matrix(ot$primal("xy")$xy)

  xi <- P / outer(a, b) # density of the plan w.r.t. a x b
  A_Q <- xi %*% diag(b) # (A_Q g)(x_i) = sum_j b_j xi_ij g_j
  A_P <- t(xi) %*% diag(a) # (A_P f)(y_j) = sum_i a_i xi_ij f_i
  Abar_P <- A_P - matrix(1, m, 1) %*% (t(b) %*% A_P) # centred in L2_0(Q)

  Sig_y <- (xi - 1) %*% diag(b^2) %*% t(xi - 1) # cov of Z_y at the x_i
  Sig_x <- t(xi - 1) %*% diag(a^2) %*% (xi - 1) # cov of Z_x at the y_j

  M <- solve(diag(n) - A_Q %*% Abar_P)
  L_fx <- eps * M %*% A_Q
  L_fy <- -eps * M
  L_x <- rbind(L_fx, -eps * diag(m) - Abar_P %*% L_fx)
  L_y <- rbind(L_fy, -Abar_P %*% L_fy)

  L_x %*% Sig_x %*% t(L_x) + L_y %*% Sig_y %*% t(L_y)
}

testthat::test_that("vcov matches the GLN (2022) plug-in covariance, not debiased", {
  causalOT:::torch_check()
  dat <- vcov_test_data()

  # GLN assumes the marginals hold exactly, so solve Sinkhorn tightly
  ot <- vcov_test_ot(dat$x, dat$y, debias = FALSE, tol = 1e-12)
  V <- as.matrix(ot$vcov())
  expect_valid_vcov(V, 18)
  testthat::expect_equal(V, gln_vcov(ot), tolerance = 1e-10)

  # non-uniform weights on y
  set.seed(4)
  b <- stats::runif(8)
  ot <- vcov_test_ot(dat$x, dat$y, b = b / sum(b), debias = FALSE, tol = 1e-12)
  testthat::expect_equal(as.matrix(ot$vcov()), gln_vcov(ot), tolerance = 1e-10)
})

testthat::test_that("vcov respects the gauge sum_j b_j g_j = 0", {
  causalOT:::torch_check()
  dat <- vcov_test_data()
  ot <- vcov_test_ot(dat$x, dat$y, debias = FALSE)

  V <- as.matrix(ot$vcov())
  b <- as.numeric(ot$b)
  g <- 10 + 1:8

  # sum_j b_j g_j is always 0, so it has no variance
  testthat::expect_equal(drop(t(b) %*% V[g, g] %*% b), 0, tolerance = 1e-12)
  # but the last g still varies with the others
  testthat::expect_gt(V[18, 18], 1e-3)
})

testthat::test_that("debiased vcov has the xy block and nonzero cross terms", {
  causalOT:::torch_check()
  dat <- vcov_test_data()

  V0 <- as.matrix(vcov_test_ot(dat$x, dat$y, debias = FALSE)$vcov())
  V1 <- as.matrix(vcov_test_ot(dat$x, dat$y, debias = TRUE)$vcov())
  expect_valid_vcov(V1, 36)

  # hessian is block diagonal across problems, so the xy block is unchanged
  testthat::expect_equal(V1[1:18, 1:18], V0, tolerance = 1e-12)

  # x drives both g_yx and f_xx, y drives both f_xy and g_yy
  testthat::expect_gt(max(abs(V1[1:18, 19:36])), 1e-3)
})

testthat::test_that("vcov zeroes out points with no mass", {
  causalOT:::torch_check()
  dat <- vcov_test_data()

  b <- rep(1 / 7, 8)
  b[3] <- 0
  V <- as.matrix(vcov_test_ot(dat$x, dat$y, b = b)$vcov())
  expect_valid_vcov(V, 36)
  testthat::expect_true(all(V[10 + 3, ] == 0)) # g_yx
  testthat::expect_true(all(V[28 + 3, ] == 0)) # g_yy

  a <- rep(1 / 9, 10)
  a[5] <- 0
  V <- as.matrix(vcov_test_ot(dat$x, dat$y, a = a)$vcov())
  expect_valid_vcov(V, 36)
  testthat::expect_true(all(V[5, ] == 0)) # f_xy
  testthat::expect_true(all(V[18 + 5, ] == 0)) # f_xx
})

testthat::test_that("vcov handles a single point", {
  causalOT:::torch_check()
  dat <- vcov_test_data()
  y1 <- dat$y[1, , drop = FALSE]
  x1 <- dat$x[1, , drop = FALSE]

  # one y: g_yx fixed by the gauge and f_xy has no variance
  V <- as.matrix(vcov_test_ot(dat$x, y1, debias = FALSE)$vcov())
  expect_valid_vcov(V, 11)
  testthat::expect_equal(max(abs(V)), 0, tolerance = 1e-12)

  V <- as.matrix(vcov_test_ot(dat$x, y1, debias = TRUE)$vcov())
  expect_valid_vcov(V, 22)

  # one point on each side: nothing varies
  V <- as.matrix(vcov_test_ot(x1, y1, debias = FALSE)$vcov())
  testthat::expect_equal(V, matrix(0, 2, 2))

  V <- as.matrix(vcov_test_ot(x1, y1, debias = TRUE)$vcov())
  testthat::expect_equal(V, matrix(0, 4, 4))
})

testthat::test_that("vcov returns the OT object's dtype", {
  causalOT:::torch_check()
  dat <- vcov_test_data()
  ot <- vcov_test_ot(dat$x, dat$y, dtype = torch::torch_float())

  V <- ot$vcov()
  testthat::expect_true(V$dtype == torch::torch_float())
  testthat::expect_equal(dim(V), c(36, 36))
})
