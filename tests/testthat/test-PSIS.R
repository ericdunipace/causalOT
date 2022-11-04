testthat::test_that("PSIS diagnostics work", {
  testthat::skip_on_cran()
  causalOT:::torch_check()
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
                                    design = design, overlap = overlap)
  data$gen_data()
  test1 <- causalOT::calc_weight(x = data, estimand = "ATE", method = "NNM")
  test2 <- causalOT::calc_weight(x = data, estimand = "ATE", method = "Logistic")
  
  weights <- list(NNM = test1,
                  IPW = test2
  )
  testthat::expect_warning(ps <- lapply(weights, PSIS))
  testthat::expect_warning(ps.check <- PSIS(weights))
  testthat::expect_equivalent(ps, ps.check)
  
  diag.check <- lapply(ps, PSIS_diag)
  diag <- PSIS_diag(ps)
  testthat::expect_warning(diag.check2 <- PSIS_diag(weights))
  testthat::expect_equal(diag, diag.check)
  testthat::expect_equal(diag, diag.check2)
})

