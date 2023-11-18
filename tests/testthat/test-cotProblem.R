testthat::test_that("cotProblem CBPS", {
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "CBPS"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_silent(gs <- causalOT:::cotProblem(data = data,
                              estimand = "ATT",
                              method = method,
                              options = options))
  
  testthat::expect_true(inherits(gs, "likelihoodMethods"))
  testthat::expect_true(inherits(gs@solver, "function"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem logistic", {
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "Logistic"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_silent(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options))
  
  testthat::expect_true(inherits(gs, "likelihoodMethods"))
  testthat::expect_true(inherits(gs@solver, "function"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem probit", {
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "Probit"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_silent(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options))
  
  testthat::expect_true(inherits(gs, "likelihoodMethods"))
  testthat::expect_true(inherits(gs@solver, "function"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem SBW", {
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "SBW"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  mess <- testthat::capture_output(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options))
  
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(inherits(gs@solver, "SBW"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem EBW", {
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "EntropyBW"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_silent(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options))
  
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(inherits(gs@solver, "EntropyBW"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem COT", {
  causalOT:::torch_check()
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "COT"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_silent(testthat::expect_warning(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options)))
  
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(inherits(gs@solver, "COT"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})

testthat::test_that("cotProblem NNM", {
  causalOT:::torch_check()
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "NNM"
  options = list(niter = 2L,
                 nboot = 2,
                 debias = TRUE,
                 torch.optimizer = torch::optim_lbfgs)
  
  testthat::expect_warning(gs <- causalOT:::cotProblem(data = data,
                                                      estimand = "ATT",
                                                      method = method,
                                                      options = options))
  
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(inherits(gs@solver, "COT"))
  testthat::expect_equal(gs@estimand, "ATT")
  testthat::expect_equal(gs@method, method)
})
