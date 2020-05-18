testthat::test_that("sim.function works", {
  set.seed(224893390) #from random.org
  
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
  trunc <- std_mean_diff <- c(0.001, 0.01, 0.1)
  agumentation <- match <- "both"
  solver <- "gurobi"
  grid.search <- TRUE
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  #### Simulations ####
  # debugonce(sim.function)
  testthat::expect_s3_class(output <- sim.function(dataGen = original, 
                          nsims = nsims, 
                          ground_p = ground_power, 
                          p = power, 
                          grid.search = grid.search,
                          augmentation = agumentation,
                          match = match,
                          standardized.mean.difference = std_mean_diff,
                          truncations = trunc,
                          distance = distance, 
                          calculate.feasible = FALSE,
                          solver = solver), "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
  
  # testthat::expect_silent
  (testthat::expect_s3_class(output <- sim.function(dataGen = original, 
                         nsims = nsims, 
                         ground_p = ground_power, 
                         p = power, 
                         grid.search = grid.search,
                         augmentation = agumentation,
                         match = match,
                         standardized.mean.difference = std_mean_diff,
                         truncations = trunc,
                         distance = distance, 
                         calculate.feasible = TRUE,
                         solver = solver), "simOutput"))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.03))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
})
