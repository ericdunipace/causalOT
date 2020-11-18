testthat::test_that("SimHolder generates object", {
  set.seed(9867)
  
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
  solver <- "gurobi"
  
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
            solver = "gurobi",
            wass_powers = 2,
            ground_powers = 2,
            metrics = "Lp")
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})

testthat::test_that("SimHolder generates object, Kallus2018", {
  set.seed(9867)
  
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
  solver <- "gurobi"
  
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
                      solver = "gurobi",
                      wass_powers = 2,
                      ground_powers = 2,
                      metrics = "Lp")
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})

testthat::test_that("SimHolder generates object, Sonabend2020", {
  set.seed(9867)
  
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
  solver <- "gurobi"
  
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
                      solver = "gurobi",
                      wass_powers = 2,
                      ground_powers = 2,
                      metrics = "Lp")
  testthat::expect_equivalent(class(sh), c("SimHolder", "R6"))
})


testthat::test_that("SimHolder runs", {
  set.seed(9867)

  #### Load Packages ####
  library(causalOT)

  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis","RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0.2,0.3)
  solver <- "gurobi"

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
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance,
                      constrained.wasserstein.target = c("SBW"))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
      {
        
        sh$run()
        warn <- warnings()
      }
    )
  if(!is.null(warn)) print(warn)
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
})

testthat::test_that("SimHolder runs,verbose", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis","RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0.2,0.3)
  solver <- "gurobi"
  
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
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance,
                      constrained.wasserstein.target = c("SBW"),
                      verbose = TRUE)
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      testthat::expect_message(sh$run())
      warn <- warnings()
    }
  )
  if(!is.null(warn)) print(warn)
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
})

testthat::test_that("SimHolder runs while targeting RKHS", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0.2, 0.3)
  solver <- "gurobi"
  
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
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance,
                      constrained.wasserstein.target = c("RKHS"))
  # the cost of one was all NA and the weights too...
  # sh$run()
  testthat::expect_warning(
    {
      
      sh$run()
      warn <- warnings()
    }
  )
  if(!is.null(warn)) print(warn)
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
})

testthat::test_that("SimHolder with grid works", {
  set.seed(082374)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0, 0.1, 1)
  solver <- "gurobi"
  
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
                      truncations = std_mean_diff,
                      standardized.difference.means = std_mean_diff,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance)
  testthat::expect_warning(
    {
      sh$run()
      warn <- warnings()
    })
  if(!is.null(warn)) print(warn)
  sh2 <- SimHolder$new(nsim = nsims,
                      dataSim = original,
                      grid.search = TRUE,
                      truncations = std_mean_diff,
                      standardized.difference.means = NULL,
                      outcome.model = list("lm"),
                      outcome.formula = list(none = NULL,
                                             augmentation = NULL),
                      model.augmentation = "both",
                      match = "both",
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance)
  testthat::expect_warning(
    {
      sh2$run()
      warn <- warnings()
    })
  if(!is.null(warn)) print(warn)
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
})

testthat::test_that("SimHolder with grid works, opt.hyperparam", {
  set.seed(082374)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0, 0.1, 1)
  solver <- "gurobi"
  
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
                      solver = "gurobi",
                      wass_powers = power,
                      ground_powers = ground_power,
                      metrics = distance)
  testthat::expect_warning(
    {
      sh$run()
      warn <- warnings()
    })
  if(!is.null(warn)) print(warn)
  sh2 <- SimHolder$new(nsim = nsims,
                       dataSim = original,
                       grid.search = TRUE,
                       RKHS = list(opt = TRUE, opt.method = "optim"),
                       truncations = std_mean_diff,
                       standardized.difference.means = NULL,
                       outcome.model = list("lm"),
                       outcome.formula = list(none = NULL,
                                              augmentation = NULL),
                       model.augmentation = "both",
                       match = "both",
                       solver = "gurobi",
                       wass_powers = power,
                       ground_powers = ground_power,
                       metrics = distance)
  testthat::expect_warning(
    {
      sh2$run()
      warn <- warnings()
    })
  if(!is.null(warn)) print(warn)
  testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
  testthat::expect_type(original$get_x0(), "double")
  testthat::expect_type(original$get_x1(), "double")
  testthat::expect_type(original$get_z(), "double")
  testthat::expect_type(original$get_y(), "double")
})


# testthat::test_that("SimHolder with with multiple kernels works, RKHS", {
#   set.seed(082374)
#   
#   #### Load Packages ####
#   library(causalOT)
#   
#   #### Sim param ####
#   n <- 2^6
#   p <- 6
#   nsims <- 2
#   overlap <- "high"
#   design <- "A"
#   distance <- c("mahalanobis")
#   power <- c(2)
#   ground_power <- 2
#   std_mean_diff <- c(0, 0.1, 1)
#   solver <- "gurobi"
#   
#   #### get simulation functions ####
#   original <- Hainmueller$new(n = n, p = p, 
#                               design = design, overlap = overlap)
#   # SimHolder$debug("initialize")
#   # SimHolder$debug("update")
#   # SimHolder$debug("estimate")
#   # SimHolder$debug("get_delta")
#   # SimHolder$debug("method.setup")
#   # SimHolder$debug("cost.setup")
#   # SimHolder$debug("max.cond.calc")
#   # SimHolder$debug("weight.calc")
#   sh <- SimHolder$new(nsim = nsims,
#                       dataSim = original,
#                       grid.search = TRUE,
#                       RKHS = list(opt = TRUE, opt.method = "stan", iter = 10,
#                                   kernel = c("RBF","linear","polynomial")),
#                       truncations = std_mean_diff,
#                       standardized.difference.means = std_mean_diff,
#                       outcome.model = list("lm"),
#                       outcome.formula = list(none = NULL,
#                                              augmentation = NULL),
#                       model.augmentation = "both",
#                       match = "both",
#                       solver = "gurobi",
#                       wass_powers = power,
#                       ground_powers = ground_power,
#                       metrics = distance)
#   testthat::expect_warning(
#     {
#       sh$run()
#       warn <- warnings()
#     })
#   if(!is.null(warn)) print(warn)
#   sh2 <- SimHolder$new(nsim = nsims,
#                        dataSim = original,
#                        grid.search = TRUE,
#                        RKHS = list(opt = TRUE, opt.method = "optim"),
#                        truncations = std_mean_diff,
#                        standardized.difference.means = NULL,
#                        outcome.model = list("lm"),
#                        outcome.formula = list(none = NULL,
#                                               augmentation = NULL),
#                        model.augmentation = "both",
#                        match = "both",
#                        solver = "gurobi",
#                        wass_powers = power,
#                        ground_powers = ground_power,
#                        metrics = distance)
#   testthat::expect_warning(
#     {
#       sh2$run()
#       warn <- warnings()
#     })
#   if(!is.null(warn)) print(warn)
#   testthat::expect_equal(class(sh$get.output()), c("data.table", "data.frame"))
#   testthat::expect_type(original$get_x0(), "double")
#   testthat::expect_type(original$get_x1(), "double")
#   testthat::expect_type(original$get_z(), "double")
#   testthat::expect_type(original$get_y(), "double")
# })
