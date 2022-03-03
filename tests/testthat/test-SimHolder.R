testthat::skip_on_cran()
testthat::skip("Interactive only")

warn.fun <- function() {
  warn <- warnings()
  pw <- c("Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.\n",
          "Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.\n",
          "Warning: Not enough tail samples to fit the generalized Pareto distribution.\n")
  if(!is.null(warn) ) {
    if(!all(names(warn) %in% pw)){
      idx <- !(names(warn) %in% pw)
      print(warn[idx])
    }
  } 
}

methods <- c('Logistic', 
             'SBW', 
             'SCM',
             'CBPS', 
             'RKHS',
             'NNM', 
             "Wasserstein", 
             # 'Constrained Wasserstein',
             'None',
             'gp'
)

testthat::test_that("SimHolder generates object", {
  set.seed(9867)
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  #### Load Packages ####
  library(causalOT)

  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  sh <- SimHolder$new(nsim = nsims,
            dataSim = original,
            grid.search = TRUE,
            truncations = std_mean_diff,
            standardized.difference.means = std_mean_diff,
            outcome.model = list("lm"),
            outcome.formula = list(none = NULL,
                                  augmentation = NULL),
            model.augmentation = "both",
            match = "both",
            solver = "mosek",
            Wass = list(wass_powers = 2,
                        ground_powers = 2,
                        metrics = "Lp"))
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
  
  
  psform <-  list(Logistic = list("z~."),
                  SBW = list("~. + 0"),
                  Wasserstein = list(NA, "~. + 0"),
                  "Constrained Wasserstein" = list(NA, "~. + 0"))
  # SimHolder$debug("method.setup")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      propensity.formula = psform,
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = 2,
                                  ground_powers = 2,
                                  metrics = "Lp",
                                  add.joint = TRUE, 
                                  add.margins = c(FALSE, TRUE)))
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})

testthat::test_that("SimHolder generates object, Kallus2018", {
  set.seed(9867)
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 4
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Kallus2018$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = 2,
                      ground_powers = 2,
                      metrics = "Lp"))
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})

testthat::test_that("SimHolder generates object, Sonabend2020", {
  set.seed(9867)
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 4
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Sonabend2020$new(n = n, p = p, 
                             design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = 2,
                      ground_powers = 2,
                      metrics = "Lp"))
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})

testthat::test_that("SimHolder runs", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)

  #### Load Packages ####
  library(causalOT)

  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"

  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = FALSE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                        ground_powers = ground_power,
                        metrics = distance,
                        constrained.wasserstein.target = c("SBW"),
                        add.margins = c(TRUE, FALSE),
                        joint.mapping = c(FALSE, TRUE),
                        penalty = c("none", "L2", "entropy", "variance")
                      ))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
      {
        
        sh$run()
        warn <- warnings()
      }
    )
  if (!is.null(warn)) warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  testthat::expect_silent(
    {out <- sh$get.output()
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out)
  wass <- sh$get.wass(out)}
  )
  testthat::expect_equal(unique(out$method), methods)
  testthat::expect_true("E_Y1" %in% colnames(out))
  testthat::expect_true("E_Y0" %in% colnames(out))
})

testthat::test_that("SimHolder runs, only ATE", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(234028)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = FALSE,
                      truncations = std_mean_diff,
                      estimands = "ATE",
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      # verbose = TRUE,
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  add.margins = c(FALSE)
                      ))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  out <- sh$get.output()
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out)
  wass <- sh$get.wass(out)
  testthat::expect_equal(unique(out$method ), methods)
  testthat::expect_true("E_Y1" %in% colnames(out))
  testthat::expect_true("E_Y0" %in% colnames(out))
  
  
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      estimands = "ATE",
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  add.margins = c(FALSE)
                      ))
  
  
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
})

testthat::test_that("SimHolder runs ot imputer", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = FALSE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      verbose = TRUE,
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  add.margins = c(TRUE, FALSE)
                      ))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  out <- sh$get.output()
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out)
  wass <- sh$get.wass(out)
  testthat::expect_equal(unique(out$method ), methods)
})

testthat::test_that("SimHolder runs with formula options", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  methods <- c("Logistic","SBW","NNM","Wasserstein")
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      methods = methods,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = solver,
                      propensity.formula = list(Logistic = list("z ~ . + .*.",
                                                                "z ~ . "),
                                                SBW = list("~.+0",
                                                           "~ . + .*. + I(.^2) + 0"),
                                                Wasserstein = list(NA, "~.+0")),
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  niter = 10,
                                  method = "greenkhorn",
                                  constrained.wasserstein.target = c("SBW")
                      ))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  out <- sh$get.output()
  testthat::expect_equal(unique(out$formula), 
                         c("z ~ . + .*.", "z ~ . ",  "~.+0",  "~ . + .*. + I(.^2) + 0",  NA ))
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out); wass <- sh$get.wass(out)
  testthat::expect_equal(unique(out$method ), c('Logistic', 'SBW',"NNM",'Wasserstein'))
})

testthat::test_that("SimHolder runs,verbose", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(1)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  methods <- c("Logistic", "SBW", "NNM", "Wasserstein", "Constrained Wasserstein")
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = FALSE,
                      methods = methods,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                          ground_powers = ground_power,
                          metrics = distance,
                          constrained.wasserstein.target = c("SBW")),
                      verbose = TRUE)
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      file.create("temp.txt")
      sink(file = "temp.txt")
      testthat::expect_message(sh$run())
      warn <- warnings()
      sink()
      file.remove("temp.txt")
    }
  )
  # if (!is.null(warn)) print(warn)
  warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  testthat::expect_silent({
  out <- sh$get.output()
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out); wass <- sh$get.wass(out)
  })
})

testthat::test_that("SimHolder runs while targeting RKHS", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("RKHS")
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0.2, 0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = FALSE,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list( wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance,
                      constrained.wasserstein.target = c("RKHS")) )
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  # if(!is.null(warn)) print(warn)
  warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  testthat::expect_silent({
  out <- sh$get.output()
  outcome <- sh$get.outcome(out)
  ess <- sh$get.ESS.frac(out)
  diag <- sh$get.diagnostics(out)
  psis <- sh$get.psis(out); wass <- sh$get.wass(out)
  })
})

testthat::test_that("SimHolder with grid works", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(082374)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c( "mahalanobis")
  power <- c(1)
  ground_power <- 2
  std_mean_diff <- c(0, 0.1, 1)
  solver <- "mosek"
  methods <- c("SBW","Wasserstein",
               "SCM")
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("max.cond.calc")
  # SimHolder$debug("weight.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      methods = methods,
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                        ground_powers = ground_power,
                        metrics = distance,
                        niter = 5,
                        method = "greenkhorn",
                        wasserstein.distance.constraints = c(10,11)),
                      verbose = FALSE)
  # options(error = recover)
  # sh$run()
  testthat::expect_warning(
    {
      sh$run()
      # warn <- warnings()
    })
  # if(!is.null(warn)) print(warn)
  warn.fun()
  sh2 <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      methods = methods,
                      standardized.difference.means = NULL,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                        ground_powers = ground_power,
                        metrics = distance,
                        niter = 5,
                        method = "sinkhorn",
                        wasserstein.distance.constraints = NULL,
                        joint.mapping = TRUE,
                        penalty = c("none","L2","variance","entropy"))
                      , verbose = TRUE
                      )
  testthat::expect_warning(
    {
      sh2$run()
      warn <- warnings()
    })
  # if(!is.null(warn)) print(warn)
  warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  testthat::expect_equal(unique(sh$get.output()$method), c("SBW","Wasserstein",
                                                           # "Constrained Wasserstein",
                                                           "SCM"))
  testthat::expect_true(unique(sh2$get.output()$penalty) %in%
                         c(NA, "none","L2","variance", "entropy"))
  testthat::expect_equal(unique(sh2$get.output()$method), c("SBW","Wasserstein",
                                                           # "Constrained Wasserstein",
                                                           "SCM"))
  
  
  sh3 <- SimHolder$new(nsim = nsims,
                       dataSim = original,
                       grid.search = TRUE,
                       truncations = std_mean_diff,
                       methods = methods,
                       standardized.difference.means = NULL,
                       outcome.model = list("lm"),
                       outcome.formula = list(none = NULL,
                                              augmentation = NULL),
                       model.augmentation = "both",
                       match = "both",
                       solver = "mosek",
                       Wass = list(wass_powers = power,
                                   ground_powers = ground_power,
                                   metrics = distance,
                                   niter = 5,
                                   method = "sinkhorn",
                                   wasserstein.distance.constraints = NULL,
                                   joint.mapping = FALSE,
                                   penalty = c("none"))
                       , verbose = TRUE
  )
  testthat::expect_warning(sh3$run())
  testthat::expect_equal(unique(sh3$get.output()$penalty),
                         c(NA, "none"))
  testthat::expect_equal(unique(sh3$get.output()$method),
                         c("SBW","Wasserstein",
                           # "Constrained Wasserstein",
                           "SCM"))
  
  
  
  sh4 <- SimHolder$new(nsim = 1,
                       dataSim = original,
                       grid.search = TRUE,
                       truncations = std_mean_diff,
                       methods = c("None", "SCM"),
                       standardized.difference.means = NULL,
                       outcome.model = list("lm"),
                       outcome.formula = list(none = NULL,
                                              augmentation = NULL),
                       model.augmentation = "both",
                       match = "both",
                       solver = "mosek",
                       Wass = list(wass_powers = power,
                                   ground_powers = ground_power,
                                   metrics = distance,
                                   niter = 5,
                                   method = "sinkhorn",
                                   wasserstein.distance.constraints = NULL,
                                   joint.mapping = TRUE,
                                   penalty = c("none","L2"))
                       , verbose = TRUE
  )
  # debugonce(wass_grid_search)
  testthat::expect_warning(sh4$run())
  
  # SimHolder$debug("estimate")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("initialize")
  sh5 <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      methods = c("NNM", "Wasserstein"),
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  niter = 5,
                                  method = "sinkhorn",
                                  penalty = "entropy",
                                  wasserstein.distance.constraints = c(10,11),
                                  add.divergence = c(TRUE,FALSE)),
                      verbose = TRUE)
  testthat::expect_warning(sh5$run())
  output <- sh5$get.output()
  testthat::expect_true("add.divergence" %in% colnames(output))
  testthat::expect_true(all(c(NA, FALSE, TRUE) %in% as.logical(unique(output[,"add.divergence"])$add.divergence)))
  testthat::expect_true(all(unique(output$add.divergence) %in% c(NA, TRUE, FALSE)))
  testthat::expect_true("E_Y1" %in% colnames(output))
  testthat::expect_true("E_Y0" %in% colnames(output))
  testthat::expect_true("E_Y1" %in% colnames(output))
  testthat::expect_true("E_Y0" %in% colnames(output))
  
  testthat::expect_silent(
    {out <- sh5$get.output()
    outcome <- sh5$get.outcome(out)
    ess <- sh5$get.ESS.frac(out)
    diag <- sh5$get.diagnostics(out)
    psis <- sh5$get.psis(out)
    wass <- sh5$get.wass(out)}
  )
  
  
})

testthat::test_that("SimHolder with grid works, opt.hyperparam", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(082374)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("RKHS")
  power <- c(2)
  ground_power <- 2
  std_mean_diff <- c(0, 0.1, 1)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("max.cond.calc")
  # SimHolder$debug("weight.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      RKHS = list(opt = TRUE, opt.method = "stan", iter = 10),
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                          ground_powers = ground_power,
                          metrics = distance,
                          method = "greenkhorn",
                          niter = 5,
                          wasserstein.distance.constraints = c(10,11)))
  testthat::expect_warning(
    {
      sh$run()
      warn <- warnings()
    })
  # if(!is.null(warn)) print(warn)
  warn.fun()
  # sh2 <- SimHolder$new(nsim = nsims,
  #                      dataSim = original,
  #                      grid.search = TRUE,
  #                      RKHS = list(opt = TRUE, opt.method = "optim",
  #                                  kernel = "polynomial"),
  #                      truncations = std_mean_diff,
  #                      standardized.difference.means = NULL,
  #                      outcome.model = list("lm"),
  #                      outcome.formula = list(none = NULL,
  #                                             augmentation = NULL),
  #                      model.augmentation = "both",
  #                      match = "both",
  #                      solver = "mosek",
  #                      wass_powers = power,
  #                      ground_powers = ground_power,
  #                      metrics = distance)
  # testthat::expect_warning(
  #   {
  #     sh2$run()
  #     warn <- warnings()
  #   })
  # if(!is.null(warn)) print(warn)
  # testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  # testthat::expect_type(original$get_x0(), "double")
  # testthat::expect_type(original$get_x1(), "double")
  # testthat::expect_type(original$get_z(), "double")
  # testthat::expect_type(original$get_y(), "double")
})

testthat::test_that("SimHolder runs confidence intervals", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(234028)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      methods = c("NNM", "Wasserstein"),
                      estimands = "ATE",
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  add.margins = c(FALSE),
                                  confidence.interval = "asymptotic"
                      ))
  
  
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
  
  output <- sh$get.output()
  
  testthat::expect_true(inherits(output$confidence.interval[1], "list"))
  
  
  
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 1
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  ground_power <- 1
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      methods = c("NNM"),
                      estimands = "ATE",
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "mosek",
                      Wass = list(wass_powers = power,
                                  ground_powers = ground_power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  add.margins = c(FALSE),
                                  confidence.interval = "bootstrap"
                      ))
  
  
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
  
  output <- sh$get.output()
  
  testthat::expect_true(inherits(output$confidence.interval[1], "list"))
  testthat::expect_true("E_Y1" %in% colnames(output))
  testthat::expect_true("E_Y0" %in% colnames(output))
  
  testthat::expect_silent(
    {out <- sh$get.output()
    outcome <- sh$get.outcome(out)
    ess <- sh$get.ESS.frac(out)
    diag <- sh$get.diagnostics(out)
    psis <- sh$get.psis(out)
    wass <- sh$get.wass(out)}
  )
  
})

testthat::test_that("SimHolder wass entropy turns to lbfgs", {
  testthat::skip_on_cran()
  # testthat::skip_if_not_installed("gurobi")
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  # testthat::skip("Interactive only")
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^5
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("sdLp")
  power <- c(2)
  std_mean_diff <- c(0.2,0.3)
  solver <- "mosek"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p,
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  # SimHolder$debug("update")
  # SimHolder$debug("estimate")
  # SimHolder$debug("model_estimate")
  # SimHolder$debug("get_delta")
  # SimHolder$debug("method.setup")
  # SimHolder$debug("cost.setup")
  # SimHolder$debug("get_cost")
  # SimHolder$debug("max.cond.calc")
  sh <- SimHolder$new(nsim = nsims,
                      methods = c("SBW","Wasserstein"),
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      RKHS = list(opt = FALSE),
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      solver = solver,
                      verbose = TRUE,
                      Wass = list(wass_powers = power,
                                  metrics = distance,
                                  constrained.wasserstein.target = c("SBW"),
                                  penalty = c("entropy")
                      ))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if (!is.null(warn)) warn.fun()
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
  
  testthat::expect_silent(
    {out <- sh$get.output()
    outcome <- sh$get.outcome(out)
    ess <- sh$get.ESS.frac(out)
    diag <- sh$get.diagnostics(out)
    psis <- sh$get.psis(out)
    wass <- sh$get.wass(out)}
  )
  testthat::expect_equal(unique(out$method), c("SBW","Wasserstein"))
  testthat::expect_true("E_Y1" %in% colnames(out))
  testthat::expect_true("E_Y0" %in% colnames(out))
  testthat::expect_true("lbfgs" %in% out$solver)
})

