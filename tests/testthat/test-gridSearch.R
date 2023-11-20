testthat::test_that("gridSearch class forms", {
  causalOT:::torch_check()
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(128)
  hain$gen_data()
  data <- dataHolder(hain)
  
  # by raw class
  testthat::expect_silent(gs <- new("gridSearch", penalty_list = 0, nboot = 100L,
  solver = hain,
  method = "COT",
  estimand = "ATE"))
  
  testthat::expect_true(inherits(gs, "gridSearch"))
  
})

testthat::test_that("gridSearch SBW", {
  testthat::skip_if_not_installed("osqp")
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "SBW"
  
  # by function SBW
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                              estimand = "ATT",
                              method = method,
                              options = list()))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_equal(gs@nboot, 1000L)
  testthat::expect_true(inherits(gs@penalty_list, "numeric"))
  testthat::expect_true(length(gs@penalty_list) == 20)
  testthat::expect_true(inherits(gs@solver, "SBW"))
  testthat::expect_true(gs@method == "SBW")
  testthat::expect_true(gs@estimand == "ATT")
  mess <- testthat::capture_output(suppressWarnings(cw <- causalOT:::cot_solve(gs)))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                              estimand = "ATC",
                              method = method,
                              options = list()))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(gs@estimand == "ATC")
  mess <- testthat::capture_output(suppressWarnings(cw <- causalOT:::cot_solve(gs)))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                              estimand = "ATE",
                              method = method,
                              options = list()))
  testthat::expect_true(inherits(gs, "ateClass"))
  
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
})

testthat::test_that("gridSearch EBW", {
  testthat::skip_if_not_installed("lbfgsb3c")
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "EntropyBW"
  
  # by function SBW
  gs <- causalOT:::gridSearch(data = data,
         estimand = "ATT",
         method = method,
         options = list())
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_equal(gs@nboot, 1000L)
  testthat::expect_true(inherits(gs@penalty_list, "numeric"))
  testthat::expect_true(length(gs@penalty_list) == 20)
  testthat::expect_true(inherits(gs@solver, "EntropyBW"))
  testthat::expect_true(gs@method == "EntropyBW")
  testthat::expect_true(gs@estimand == "ATT")
  
  mess <- testthat::capture_output(testthat::expect_warning(cw <- causalOT:::cot_solve(gs)))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATC",
                                                               method = method,
                                                               options = list()))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(gs@estimand == "ATC")
  mess <- testthat::capture_output(testthat::expect_warning(cw <- causalOT:::cot_solve(gs)))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATE",
                                                               method = method,
                                                               options = list()))
  testthat::expect_true(inherits(gs, "ateClass"))
  
  mess <- testthat::capture_output(testthat::expect_warning(cw <- causalOT:::cot_solve(gs)))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
})

testthat::test_that("gridSearch COT", {
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
  
  # by function
  gs <- testthat::expect_warning(causalOT:::gridSearch(data = data,
                              estimand = "ATT",
                              method = method,
                              options = options))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_equal(gs@nboot, 2)
  testthat::expect_true(inherits(gs@penalty_list, "numeric"))
  testthat::expect_true(length(gs@penalty_list) == 8)
  testthat::expect_true(inherits(gs@solver, "COT"))
  testthat::expect_true(gs@method == "COT")
  testthat::expect_true(gs@estimand == "ATT")
  
  cw <- causalOT:::cot_solve(gs)
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(testthat::expect_warning(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATC",
                                                               method = method,
                                                               options = options)))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(gs@estimand == "ATC")
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(testthat::expect_warning(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATE",
                                                               method = method,
                                                               options = options)))
  testthat::expect_true(inherits(gs, "ateClass"))
  
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
})

testthat::test_that("gridSearch NNM", {
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
  
  # by function
  gs <- testthat::expect_warning(causalOT:::gridSearch(data = data,
                              estimand = "ATT",
                              method = method,
                              options = options))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_equal(gs@nboot, 2)
  testthat::expect_true(inherits(gs@penalty_list, "numeric"))
  testthat::expect_true(length(gs@penalty_list) == 1)
  testthat::expect_true(inherits(gs@solver, "COT"))
  testthat::expect_true(gs@method == "NNM")
  testthat::expect_true(gs@estimand == "ATT")
  
  cw <- causalOT:::cot_solve(gs)
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(testthat::expect_warning(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATC",
                                                               method = method,
                                                               options = options)))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(gs@estimand == "ATC")
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(testthat::expect_warning(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATE",
                                                               method = method,
                                                               options = options)))
  testthat::expect_true(inherits(gs, "ateClass"))
  
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
})


testthat::test_that("gridSearch SCM", {
  testthat::skip_if_not_installed("osqp")
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(n=64)
  hain$gen_data()
  data <- dataHolder(hain)
  method <- "SCM"
  
  # by function
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                              estimand = "ATT",
                              method = method,
                              options = list(NULL)))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_equal(gs@nboot, 1000)
  testthat::expect_true(inherits(gs@penalty_list, "numeric"))
  testthat::expect_true(is.na(gs@penalty_list))
  testthat::expect_true(length(gs@penalty_list) == 1)
  testthat::expect_true(inherits(gs@solver, "SCM"))
  testthat::expect_true(gs@method == "SCM")
  testthat::expect_true(gs@estimand == "ATT")
  
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATC",
                                                               method = method,
                                                               options = options))
  testthat::expect_true(inherits(gs, "gridSearch"))
  testthat::expect_true(gs@estimand == "ATC")
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
  mess <- testthat::capture_output(gs <- causalOT:::gridSearch(data = data,
                                                               estimand = "ATE",
                                                               method = method,
                                                               options = options))
  testthat::expect_true(inherits(gs, "ateClass"))
  
  mess <- testthat::capture_output(cw <- causalOT:::cot_solve(gs))
  testthat::expect_true(inherits(cw, "causalWeights"))
  
})


