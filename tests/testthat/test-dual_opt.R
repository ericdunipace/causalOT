testthat::test_that("dual opt works", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  x0 <- data$get_x0()
  x  <- data$get_x()
  
  # debugonce(dual_opt)
  testthat::expect_silent({
  sbwfit <- dual_opt(x = x0, target = x, 
          init = NULL,
          sample_weights = NULL, 
          method = "SBW",
          wasserstein = list(metric = c("sdLp"),
                             power = 2,
                             cost = NULL,
                             delta = 1.0),
          balance = list(balance.functions = x0,
                         formula = NULL,
                         balance.delta = 0.001),
          marginal.wasserstein = list(marginal.costs = NULL,
                                      marginal.constraints = 1),
            control = NULL)
  
  
  # debugonce(dual_opt)
  wfit <- dual_opt(x = x0, target = x, 
                  init = NULL,
                  sample_weights = NULL, 
                  method = c("Wasserstein"),
                  wasserstein = list(metric = c("mahalanobis"),
                                     power = 2,
                                     cost = NULL,
                                     delta = 1.0),
                  balance = NULL,
                  marginal.wasserstein = NULL,
                  control = NULL)
  
  wfit <- dual_opt(x = x0, target = x, 
                   init = NULL,
                   sample_weights = NULL, 
                   method = c("Wasserstein"),
                   wasserstein = list(metric = c("mahalanobis"),
                                      power = 2,
                                      cost = NULL,
                                      delta = 1),
                   balance = NULL,
                   marginal.wasserstein = list(marginal.costs = NULL,
                                               marginal.constraints = 0.5),
                   control = list(maxit = 10000))
  
  wfit <- dual_opt(x = x0, target = x, 
                   init = NULL,
                   sample_weights = NULL, 
                   method = c("Wasserstein"),
                   wasserstein = list(metric = c("mahalanobis"),
                                      power = 2,
                                      cost = NULL,
                                      delta = 1.0),
                   balance = list(balance.functions = x0,
                                  formula = NULL,
                                  balance.delta = 0.1),
                   marginal.wasserstein = NULL,
                   control = NULL)
  
  wfit <- dual_opt(x = x0, target = x, 
                   init = NULL,
                   sample_weights = NULL, 
                   method = c("Wasserstein"),
                   wasserstein = list(metric = c("mahalanobis"),
                                      power = 2,
                                      cost = NULL,
                                      delta = 1.0),
                   balance = list(balance.functions = x0,
                                  formula = NULL,
                                  balance.delta = 0.1),
                   marginal.wasserstein = list(marginal.costs = NULL,
                                               marginal.constraints = 1),
                   control = NULL)
  
  # debugonce(dual_opt)
  cwfit <- dual_opt(x = x0, target = x, 
                  init = NULL,
                  sample_weights = NULL, 
                  method = c("Constrained Wasserstein"),
                  wasserstein = list(metric = c("mahalanobis"),
                                     power = 2,
                                     cost = NULL,
                                     delta = 1.0),
                  balance = NULL,
                  marginal.wasserstein = NULL,
                  control = NULL)
  
  # debugonce(dual_opt)
  cwfit <- dual_opt(x = x0, target = x, 
                    init = NULL,
                    sample_weights = NULL, 
                    method = c("Constrained Wasserstein"),
                    wasserstein = list(metric = c("mahalanobis"),
                                       power = 2,
                                       cost = NULL,
                                       delta = 1.0),
                    balance = list(balance.functions = x0,
                                   formula = NULL,
                                   balance.delta = 0.001),
                    marginal.wasserstein = NULL,
                    control = NULL)
  
  cwfit <- dual_opt(x = x0, target = x, 
                    init = NULL,
                    sample_weights = NULL, 
                    method = c("Constrained Wasserstein"),
                    wasserstein = list(metric = c("mahalanobis"),
                                       power = 2,
                                       cost = NULL,
                                       delta = 1.0),
                    balance = NULL,
                    marginal.wasserstein = list(marginal.costs = NULL,
                                                marginal.constraints = 1),
                    control = NULL)
  
  
  cwfit <- dual_opt(x = x0, target = x, 
                    init = NULL,
                    sample_weights = NULL, 
                    method = c("Constrained Wasserstein"),
                    wasserstein = list(metric = c("mahalanobis"),
                                       power = 2,
                                       cost = NULL,
                                       delta = 1.0),
                    balance = list(balance.functions = x0,
                                   formula = NULL,
                                   balance.delta = 0.001),
                    marginal.wasserstein = list(marginal.costs = NULL,
                                                marginal.constraints = 1),
                    control = NULL)
  
})
  
})
