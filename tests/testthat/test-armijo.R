testthat::test_that("armijo takes the full step when it is good enough", {
  causalOT:::torch_check()
  phi <- function(x, dx, alpha) {
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor((z - 1)^2, dtype = torch::torch_double())
  }

  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 1, derphi0 = -2, x = 0, dx = 1)
  testthat::expect_equal(as.numeric(res$alpha), 1)
  testthat::expect_equal(as.numeric(res$phi1), 0)
})

testthat::test_that("armijo quadratic step is exact on a quadratic", {
  causalOT:::torch_check()
  evals <- 0
  phi <- function(x, dx, alpha) {
    evals <<- evals + 1
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor((z - 0.1)^2, dtype = torch::torch_double())
  }

  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 0.01, derphi0 = -0.2, x = 0, dx = 1)
  testthat::expect_equal(as.numeric(res$alpha), 0.1, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(res$phi1), 0, tolerance = 1e-6)
  testthat::expect_equal(evals, 2)
})

testthat::test_that("armijo cubic step is exact on a cubic", {
  causalOT:::torch_check()
  evals <- 0
  phi <- function(x, dx, alpha) {
    evals <<- evals + 1
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor(1 - z + 10 * z^2 - 8 * z^3, dtype = torch::torch_double())
  }

  # the quadratic step at 0.25 raises the loss, the cubic finds the local minimum
  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 1, derphi0 = -1, x = 0, dx = 1)
  alpha <- as.numeric(res$alpha)
  testthat::expect_equal(alpha, (10 - sqrt(76)) / 24, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(res$phi1), 1 - alpha + 10 * alpha^2 - 8 * alpha^3, tolerance = 1e-6)
  testthat::expect_equal(evals, 3)
})

testthat::test_that("armijo keeps backtracking until the condition holds", {
  causalOT:::torch_check()
  evals <- 0
  phi <- function(x, dx, alpha) {
    evals <<- evals + 1
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor(ifelse(z <= 1e-3, 1 - z, 5), dtype = torch::torch_double())
  }

  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 1, derphi0 = -1, x = 0, dx = 1)
  alpha <- as.numeric(res$alpha)
  testthat::expect_gt(alpha, 0)
  testthat::expect_lte(alpha, 1e-3)
  testthat::expect_lte(as.numeric(res$phi1), 1 - 1e-4 * alpha)
  testthat::expect_equal(as.numeric(res$phi1), 1 - alpha, tolerance = 1e-6)
  testthat::expect_equal(evals, 4)
})

testthat::test_that("armijo returns a NULL step when nothing decreases", {
  causalOT:::torch_check()
  phi <- function(x, dx, alpha) {
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor(z, dtype = torch::torch_double())
  }

  # slope says descent but the function only goes up
  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 0, derphi0 = -1, x = 0, dx = 1)
  testthat::expect_null(res$alpha)
})

testthat::test_that("armijo stops when the interpolated step is NaN", {
  causalOT:::torch_check()
  evals <- 0
  phi <- function(x, dx, alpha) {
    evals <<- evals + 1
    if (evals > 50) stop("armijo did not stop")
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor(ifelse(abs(z - 1) < 1e-9, 10, NaN), dtype = torch::torch_double())
  }

  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = 1, derphi0 = -1, x = 0, dx = 1)
  testthat::expect_null(res$alpha)
  testthat::expect_equal(evals, 3)
})

testthat::test_that("armijo does not round phi0 to float32", {
  causalOT:::torch_check()
  phi0 <- 1 + 1e-7
  phi <- function(x, dx, alpha) {
    z <- x + dx * as.numeric(alpha)
    torch::torch_tensor(ifelse(abs(z - 1) < 1e-12, phi0 + 5, phi0 + 1e-8), dtype = torch::torch_double())
  }

  # every step raises the loss, but only by less than float32 can see
  res <- causalOT:::scalar_search_armijo(phi = phi, phi0 = phi0, derphi0 = -1e-6, x = 0, dx = 1)
  testthat::expect_null(res$alpha)
})
