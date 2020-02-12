test_that("SimHolder generates object", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)

  #### Sim param ####
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp", "mahalanobis")
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


test_that("SimHolder generates object", {
  set.seed(9867)
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp", "mahalanobis")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "gurobi"
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  # SimHolder$debug("initialize")
  SimHolder$debug("update")
  SimHolder$debug("estimate")
  # SimHolder$debug("method.setup")
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
                      metrics = distance)
  sh$run()
  
})
