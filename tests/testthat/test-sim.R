warn.fun <- function() {
  warn <- warnings()
  pw <- c("Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.\n",
          "Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.\n")
  if(!is.null(warn) ) {
    if(!all(names(warn) %in% pw)){
      idx <- !(names(warn) %in% pw)
      print(warn[idx])
    }
  } 
}

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
  trunc <- std_mean_diff <- c(0.001, 0.01, 1)
  agumentation <- match <- "both"
  solver <- "gurobi"
  grid.search <- TRUE
  wdc <- c(10:11)
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  #### Simulations ####
  # debugonce(sim.function)
  testthat::expect_warning({
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
                         calculate.feasible = FALSE,
                         solver = solver,
                         wass.method  = "greenkhorn",
                         wass.niter = 5,
                         wasserstein.distance.constraints = wdc)
  })
  warn.fun()
  # if(!is.null(warn)) print(warn)
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.03))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] <= 1.01))
  testthat::expect_true(is.numeric(unlist(output$PSIS[,c("psis.k.control","psis.k.treated")])))
  
  # testthat::expect_silent
  testthat::expect_warning({
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
                         calculate.feasible = TRUE,
                         solver = solver,
                         wass.method  = "greenkhorn",
                         wass.niter = 5,
                         wasserstein.distance.constraints = wdc
                         )
  })
  warn.fun()
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.03))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] <= 1.01))
  testthat::expect_true(is.numeric(unlist(output$PSIS[,c("psis.k.control","psis.k.treated")])))
  
})

testthat::test_that("sim.function works, sonabend2020", {
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
  trunc <- std_mean_diff <- c(0.001, 0.01, 1)
  agumentation <- match <- "both"
  solver <- "gurobi"
  grid.search <- TRUE
  wdc <- 10:11
  
  #### get simulation functions ####
  original <- Sonabend2020$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  #### Simulations ####
  # debugonce(sim.function)
  testthat::expect_warning({
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
                           calculate.feasible = FALSE,
                           solver = solver,
                           wass.method  = "greenkhorn",
                           wass.niter = 5,
                           wasserstein.distance.constraints = wdc)
  })
  warn.fun()
  # if(!is.null(warn)) print(warn)
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.03))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] <= 1.01))
  testthat::expect_true(is.numeric(unlist(output$PSIS[,c("psis.k.control","psis.k.treated")])))
  
  # testthat::expect_silent
  testthat::expect_warning({
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
                           calculate.feasible = TRUE,
                           solver = solver,
                           wass.method  = "greenkhorn",
                           wass.niter = 5,
                           wasserstein.distance.constraints = wdc)
  })
  warn.fun()
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.36))
  testthat::expect_true(all(output$Wasserstein[,"dist"] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] <= 1.01))
  testthat::expect_true(is.numeric(unlist(output$PSIS[,c("psis.k.control","psis.k.treated")])))
  
})

testthat::test_that("sim.function works with RKHS", {
  set.seed(224893390) #from random.org
  
  #### Load Packages ####
  library(causalOT)
  
  #### Sim param ####
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp", "mahalanobis", "RKHS")
  power <- c(1,2)
  ground_power <- 2
  trunc <- std_mean_diff <- c(0.001, 0.01, 1)
  agumentation <- match <- "both"
  solver <- "gurobi"
  grid.search <- TRUE
  wdc <- 40:41
  
  #### get simulation functions ####
  original <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  
  #### Simulations ####
  # debugonce(sim.function)
  testthat::expect_warning({
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
                         calculate.feasible = FALSE,
                         solver = solver,
                         wass.method  = "greenkhorn",
                         wass.niter = 5,
                         wasserstein.distance.constraints = wdc)
  })
  warn.fun()
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.03))
  testthat::expect_true(all(unlist(output$Wasserstein[,"dist"]) >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$PSIS[,c("psis.ESS.frac.control","psis.ESS.frac.treated")] <= 1.01))
  testthat::expect_true(is.numeric(unlist(output$PSIS[,c("psis.k.control","psis.k.treated")])))
  
  
  testthat::expect_warning({
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
                           calculate.feasible = TRUE,
                           solver = solver,
                           wass.method  = "greenkhorn",
                           wass.niter = 5,
                           wasserstein.distance.constraints = wdc)
  })
  warn.fun()
  testthat::expect_s3_class(output, "simOutput")
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] >= 0))
  testthat::expect_true(all(output$`ESS/N`[,c("ESS.frac.control","ESS.frac.treated")] <= 1.14))
  testthat::expect_true(all(unlist(output$Wasserstein[,"dist"]) >= 0))
})

testthat::test_that("bug in gp code",
                    {
                      testthat::skip("Only run interactively")
                      design <- "B"
                      overlap <- "low"
                      data <- "Hainmueller"
                      seed <- 767362796
                      
                      n <- 2^9
                      p <- 6
                      nsims <- 1
                      distance <- c("mahalanobis", "Lp") #c("Lp", "mahalanobis", "RKHS")
                      wass_power <- 1 #c(1,2)
                      ground_power <- 1 #1:2
                      std_mean_diff <- seq(0, p^(-1/2), length.out = 50)
                      trunc <- c(0, 0.01, 0.05, 0.1, 0.2)
                      solver <- "gurobi"
                      augmentation <- match <- "both"
                      grid.search <- TRUE
                      RKHS <- list(opt = TRUE, opt.method = "stan")
                      calc.feasible <- FALSE #FALSE
                      cwass.targ <- c("SBW")
                      
                      if(data == "Hainmueller") {
                        dataGen <- Hainmueller$new(n = n, p = p, 
                                                   design = design, overlap = overlap)
                      } else if (data == "Sonabend"){
                        dataGen <- Sonabend2020$new(n = n, p = p, 
                                                    design = design, overlap = overlap)
                      }
                      times <- proc.time()
                      testthat::expect_warning(output <- sim.function(dataGen = dataGen,
                                             nsims = nsims, 
                                             ground_p = ground_power, 
                                             p = wass_power,
                                             grid.search = grid.search,
                                             RKHS = RKHS,
                                             match = match,
                                             estimands = c("cATE","ATE"),
                                             augmentation = augmentation,
                                             standardized.mean.difference = std_mean_diff,
                                             truncations = trunc,
                                             distance = distance, 
                                             solver = solver,
                                             calculate.feasible = calc.feasible,
                                             constrained.wasserstein.target = cwass.targ,
                                             seed = seed,
                                             verbose = FALSE))
                      
                      
                    })