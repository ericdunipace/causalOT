testthat::test_that("optimal weighting works, augmentation", {
  set.seed(6464546)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment = FALSE
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  # debugonce( data$opt_weight)
  opt_weights_mosek <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "mosek"))
  opt_weights_gurobi<- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "gurobi"))
  opt_weights_cplex <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "cplex"))
  names(opt_weights_mosek) <-
    names(opt_weights_gurobi) <- 
    names(opt_weights_cplex) <- estimates
  
  m0 <- lm.fit(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop = FALSE], data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE])
  m1 <- lm.fit(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop = FALSE], data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE])
  
  mu_1 = switch(paste0(c(data$.__enclos_env__$private$design, data$.__enclos_env__$private$overlap), collapse =", "),
                "A, high" = 2.0944332,
                "A, low" = 2.3886758,
                "B, high" = 9.003046,
                "B, low" = 9.497603 )
  mu_0 = switch(paste0(c(data$.__enclos_env__$private$design, data$.__enclos_env__$private$overlap), collapse =", "),
                "A, high" = 0.9165373,
                "A, low" = 0.6069783,
                "B, high" = 7.011926,
                "B, low" = 6.528442)
  
  mu = switch(as.character(data$.__enclos_env__$private$design),
              "A" = 1.5,
              "1" = 1.5,
              "B" = 8,
              "2" = 8,
              stop("design must be one of 'A' or 'B'"))
  
  
  #check integral
  eval_fun <- function(weight) {
    for(e in estimates){
      estimand <- e
      if(estimand == "cATE") estimand <- "ATE"
      f = switch(as.character(augment),
                 "FALSE" = switch(estimand,
                                  "ATT" = data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE],
                                  "ATC" = data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE],
                                  "ATE" = list(data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE],
                                               data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE])),
                 "TRUE" = switch(estimand,
                                 "ATT" = m0$residuals,
                                 "ATC" = m1$residuals,
                                 "ATE" = list(m0$residuals,
                                              m1$residuals)))
      const <- switch(as.character(augment),
                      "FALSE" = switch(estimand,
                                       "ATT" = 0,
                                       "ATC" = 0,
                                       "ATE" = c(0,
                                                 0)),
                      "TRUE" = switch(estimand,
                                      "ATT" = mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop= FALSE] %*% m0$coefficients),
                                      "ATC" = mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop= FALSE] %*% m1$coefficients),
                                      "ATE" = c(mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop= FALSE] %*% m0$coefficients),
                                                mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop= FALSE] %*% m1$coefficients))))
      mu_var <- switch(estimand,
                       "ATE" = mu,
                       "ATT" = mu_1,
                       "ATC" = mu_0)
      
      if(estimand != "ATE"){
        testthat::expect_equivalent(
          mu_var, c( weight[[estimand]] %*% f))
      } else {
        testthat::expect_equivalent(
          mu_var, c(weight[["ATE"]]$w0 %*% f[[1]]))
        testthat::expect_equivalent(
          mu_var, c(weight[["ATE"]]$w1 %*% f[[2]]))
      }
    }
  }
  out <- lapply(list(opt_weights_mosek,
                                    opt_weights_gurobi,
                                    opt_weights_cplex),
         eval_fun
    )
})

testthat::test_that("optimal weighting works, no augmentation", {
  set.seed(6464546)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment = TRUE
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  
  # debugonce( data$opt_weight)
  opt_weights_mosek <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "mosek"))
  opt_weights_gurobi<- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "gurobi"))
  opt_weights_cplex <- lapply(estimates, function(e) data$opt_weight(estimand = e, augment = augment, solver = "cplex"))
  names(opt_weights_mosek) <-
    names(opt_weights_gurobi) <- 
    names(opt_weights_cplex) <- estimates
  
  m0 <- lm.fit(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop = FALSE], data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE])
  m1 <- lm.fit(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop = FALSE], data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE])
  
  mu_1 = switch(paste0(c(data$.__enclos_env__$private$design, data$.__enclos_env__$private$overlap), collapse =", "),
                "A, high" = 2.0944332,
                "A, low" = 2.3886758,
                "B, high" = 9.003046,
                "B, low" = 9.497603 )
  mu_0 = switch(paste0(c(data$.__enclos_env__$private$design, data$.__enclos_env__$private$overlap), collapse =", "),
                "A, high" = 0.9165373,
                "A, low" = 0.6069783,
                "B, high" = 7.011926,
                "B, low" = 6.528442)
  
  mu = switch(as.character(data$.__enclos_env__$private$design),
              "A" = 1.5,
              "1" = 1.5,
              "B" = 8,
              "2" = 8,
              stop("design must be one of 'A' or 'B'"))
  
  
  #check integral
  eval_fun <- function(weight) {
    for(e in estimates){
      estimand <- e
      if(estimand == "cATE") estimand <- "ATE"
      f = switch(as.character(augment),
                 "FALSE" = switch(estimand,
                                  "ATT" = data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE],
                                  "ATC" = data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE],
                                  "ATE" = list(data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==0,drop=FALSE],
                                               data$.__enclos_env__$private$y[data$.__enclos_env__$private$z==1,drop=FALSE])),
                 "TRUE" = switch(estimand,
                                 "ATT" = m0$residuals,
                                 "ATC" = m1$residuals,
                                 "ATE" = list(m0$residuals,
                                              m1$residuals)))
      const <- switch(as.character(augment),
                      "FALSE" = switch(estimand,
                                       "ATT" = 0,
                                       "ATC" = 0,
                                       "ATE" = c(0,
                                                 0)),
                      "TRUE" = switch(estimand,
                                      "ATT" = mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop= FALSE] %*% m0$coefficients),
                                      "ATC" = mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop= FALSE] %*% m1$coefficients),
                                      "ATE" = c(mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==1,,drop= FALSE] %*% m0$coefficients),
                                                mean(data$.__enclos_env__$private$x[data$.__enclos_env__$private$z==0,,drop= FALSE] %*% m1$coefficients))))
      mu_var <- switch(estimand,
                       "ATE" = mu,
                       "ATT" = mu_1,
                       "ATC" = mu_0)
      
      if(estimand != "ATE"){
        testthat::expect_equivalent(
          mu_var, const + c( weight[[estimand]] %*% f))
      } else {
        testthat::expect_equivalent(
          mu_var, const[1] + c(weight[["ATE"]]$w0 %*% f[[1]]))
        testthat::expect_equivalent(
          mu_var, const[2] + c(weight[["ATE"]]$w1 %*% f[[2]]))
      }
    }
  }
  out <- lapply(list(opt_weights_mosek,
                     opt_weights_gurobi,
                     opt_weights_cplex),
                eval_fun
  )
})

testthat::test_that("optimal weighting comparison works, no augmentation", {
  set.seed(9847)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- FALSE
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = .8, 
                                                       estimand = e, 
                                                       p = power,
                                                       method = "Wasserstein",
                                                       solver = "gurobi"))
  
  
  
  # opt_weights_mosek <- data$opt_weight(estimand = "ATT", augment = augment, solver = "mosek")
  # opt_weights_gurobi<- data$opt_weight(estimand = "ATE", augment = augment, solver = "gurobi")
  # opt_weights_cplex <- data$opt_weight(estimand = "ATE", augment = augment, solver = "cplex")
  
  # debugonce(data$opt_weight_dist)
  compare_mosek <- mapply(data$opt_weight_dist, weight = weights, estimand = estimates,
                                        augment = augment, solver = "mosek")
  for(cc in compare_mosek) testthat::expect_lt(cc, 0.25)
})

testthat::test_that("optimal weighting comparison works. augmentation", {
  set.seed(9847)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "gurobi"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  augment <- TRUE
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = .8, 
                                                       estimand = e, 
                                                       p = power,
                                                       method = "Wasserstein",
                                                       solver = "gurobi"))
  
  
  
  # opt_weights_mosek <- data$opt_weight(estimand = "ATT", augment = augment, solver = "mosek")
  # opt_weights_gurobi<- data$opt_weight(estimand = "ATE", augment = augment, solver = "gurobi")
  # opt_weights_cplex <- data$opt_weight(estimand = "ATE", augment = augment, solver = "cplex")
  
  # debugonce(data$opt_weight_dist)
  compare_mosek <- mapply(data$opt_weight_dist, weight = weights, estimand = estimates,
                          augment = augment, solver = "mosek")
  for(cc in compare_mosek) testthat::expect_lt(cc, 0.55)
})
