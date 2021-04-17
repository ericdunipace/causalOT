testthat::test_that("convert sol doesn't give error", {
  sol <- list(list(sol = rep(1,10), dual = NULL))
  sample_weight <- get_sample_weight(NULL, c(rep(0,5), rep(1,5)))
  testthat::expect_silent(convert_sol(sol, estimand = "ATT", method = "SBW", 5, 5,sample_weight))
  testthat::expect_silent(convert_sol(sol, estimand = "ATC", method = "SBW", 5, 5, sample_weight))
  testthat::expect_silent(convert_sol(list(sol[[1]],sol[[1]]), estimand = "cATE", method = "SBW", 5, 5,sample_weight))
  
})

testthat::test_that("convert sol does give error", {
  sol <- list(sol = rep(1,10), dual = NULL)
  sample_weight <- get_sample_weight(NULL, c(rep(0,5), rep(1,5)))
  
  testthat::expect_error(convert_sol(sol, estimand = "ATE", method = "SBW", 5, 5,sample_weight))
  testthat::expect_error(convert_sol(sol, estimand = "ATE", method = "SBW", 5, 5,sample_weight))
  testthat::expect_error(convert_sol(sol, estimand = "ATT", method = "SBW", 5, 5,sample_weight))
  testthat::expect_error(convert_sol(sol, estimand = "ATC", method = "SBW", 5, 5,sample_weight))
  testthat::expect_error(convert_sol(sol, estimand = "cATE", method = "SBW", 5, 5,sample_weight))
})
