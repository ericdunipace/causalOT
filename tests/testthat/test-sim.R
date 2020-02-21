test_that("sim.function works", {
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
  output <- sim.function(dataGen = original, 
                          nsims = nsims, 
                          ground_p = ground_power, 
                          p = power, 
                          grid.search = grid.search,
                          augmentation = agumentation,
                          match = match,
                          standardized.mean.difference = std_mean_diff,
                          truncations = trunc,
                          distance = distance, 
                          solver = solver)
  
})
