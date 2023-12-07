testthat::test_that("calc_weight works", {
  causalOT:::torch_check()
  
  set.seed(23483)
  n <- 2^5
  p <- 6
  overlap <- "low"
  design <- "A"
  estimate <- "ATE"
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p,
        design = design, overlap = overlap)
  data$gen_data()
  weights <- calc_weight(x = data,
        z = NULL,
        estimand = estimate,
        method = "NNM")
  
  testthat::expect_equal(weights@penalty,
                         list(w0 = c(0), w1 = c(0)))
  testthat::expect_equal(weights@call, 
                         as.call(str2lang("calc_weight(x = data, z = NULL, estimand = estimate, method = 'NNM')")))
  testthat::expect_s4_class(weights, "causalWeights")
  
  mess <- testthat::capture_output(weights <- calc_weight(x = data,
                         estimand = estimate,
                         method = "SCM"))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "SCM")
  testthat::expect_s4_class(weights, "causalWeights")
  
  
  weights <- calc_weight(x = data,
                        estimand = estimate,
                        method = "Logistic")
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "Logistic")
  testthat::expect_s4_class(weights, "causalWeights")
  
  weights <- calc_weight(x = data,
                         estimand = estimate,
                         method = "Probit")
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "Probit")
  testthat::expect_s4_class(weights, "causalWeights")
  
  testthat::expect_warning(weights <- calc_weight(x = data,
                                                  estimand = estimate,
                                                  method = "EntropyBW"))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "EntropyBW")
  testthat::expect_s4_class(weights, "causalWeights")
  
  testthat::skip_on_cran()
  
  mess <- testthat::capture_output(weights <- calc_weight(x = data,
                                                          estimand = estimate,
                                                          method = "EnergyBW",
                                                          options =list(niter = 2L)))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "EnergyBW")
  testthat::expect_s4_class(weights, "causalWeights")
  
  
  mess <- testthat::capture_output(weights <- calc_weight(x = data,
                                                          estimand = estimate,
                                                          method = "COT",
                                                          options =list(lambda = 100, niter = 2L)))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "COT")
  testthat::expect_s4_class(weights, "causalWeights")
  
  testthat::skip_if_not_installed("CBPS")
  mess <- testthat::capture_output(weights <- calc_weight(x = data,
                        estimand = estimate,
                        method = "CBPS"))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "CBPS")
  testthat::expect_s4_class(weights, "causalWeights")
  
  testthat::skip_if_not_installed("osqp")
  mess <- testthat::capture_output(testthat::expect_warning(
    weights <- calc_weight(x = data,
                         estimand = estimate,
                         method = "SBW")
    ))
  testthat::expect_equal(weights@estimand, estimate)
  testthat::expect_equal(weights@method, "SBW")
  testthat::expect_s4_class(weights, "causalWeights")
  
  
  
})

