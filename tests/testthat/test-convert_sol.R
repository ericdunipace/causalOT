testthat::test_that("convert sol doesn't give error", {
  sol <- list(rep(1,10))
  testthat::expect_silent(convert_sol(sol, estimand = "ATT", method = "SBW", 5, 5))
  testthat::expect_silent(convert_sol(sol, estimand = "ATC", method = "SBW", 5, 5))
  testthat::expect_silent(convert_sol(list(sol[[1]],sol[[1]]), estimand = "cATE", method = "SBW", 5, 5))
  
})

testthat::test_that("convert sol does give error", {
  sol <- rep(1,10)
  testthat::expect_error(convert_sol(sol, estimand = "ATE", method = "SBW", 5, 5))
  testthat::expect_error(convert_sol(sol, estimand = "ATE", method = "SBW", 5, 5))
  testthat::expect_error(convert_sol(sol, estimand = "ATT", method = "SBW", 5, 5))
  testthat::expect_error(convert_sol(sol, estimand = "ATC", method = "SBW", 5, 5))
  testthat::expect_error(convert_sol(sol, estimand = "cATE", method = "SBW", 5, 5))
})
