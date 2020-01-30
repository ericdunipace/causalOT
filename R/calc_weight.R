setClass("causalWeights", slots = c(w0 = "numeric", w1="numeric", gamma = "NULL",estimate = "character"))

calc_weight <- function(data, constraint,  estimate = c("ATT", "ATC","feasible"), 
                                method = c("SBW","Wasserstein", "Constrained Wasserstein",
                                           "Logistic"),
                        transport.matrix = FALSE,
                                ...) {
  method <- match.arg(method)
  estimate <- match.arg(estimate)
  dots <- list(...)
  
  sol <- if(method != "Logistic") {
    do.call("calc_weight_bal", list(data, constraint,  estimate = estimate, 
                            method = method,
                            ...))
  } else {
    calc_weight_glm (data, constraint, estimate,...)
  }
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  gamma <- NULL
  ns <- data$get_n()
  if(method == "Wasserstein" | method == "Constrained Wasserstein") {
    gamma <- matrix(sol, ns["n0"], ns["n1"])
    if(estimate == "ATT") {
      sol <- rowSums(gamma)
    } else if (estimate == "ATC") {
      sol <- colSums(gamma)
    } else if (estimate == "feasible") {
      sol <- gamma
    }
    output$gamma <- gamma
  }
  if(estimate == "ATT") {
    output$w0 <- sol
    output$w1 <- rep.int(1/ns["n1"],ns["n1"])
  } else if (estimate == "ATC") {
    output$w0 <- rep.int(1/ns["n0"],ns["n0"])
    output$w1 <- sol
  } else if (estimate == "feasible") {
    if(method == "Wasserstein" | method == "Constrained Wasserstein"){
      output$w0 <- rowSums(gamma)
      output$w1 <- colSums(gamma)
    } else if (method == "SBW") {
      output$w0 <- renormalize(res$xopt[1:ns["n0"]])
      output$w1 <- renormalize(res$xopt[ns["n0"] + 1:ns["n1"]])
    }
  }
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, ...)
  }
  output$estimate <- estimate
  class(output) <- "causalWeights"
  return(output)
}

calc_weight_glm<- function(data, constraint,  estimate = c("ATT", "ATC","ATE"),
                            ...) {
  dots <- list(...)
  pd <- prep_data(data,...)
  z <- pd$z
  df <- pd$df
  
  n1 <- sum(z)
  n0 <- sum(1-z)
  n  <- n1 + n0
  if(any(colnames(df) == "y")) {
    df$y <- NULL
  }
  if(is.null(dots$formula)) {
    dots$formula <- formula(z ~ .)
  }
  mod <- glm(dots$formula, data.frame(z = z, df), family = "binomial")
  pred <- predict(mod, type = "response")
  
  if (constraint > 0 & constraint < 1) {
    Ks  <- sort(c(constraint, 1-constraint))
    up  <- Ks[2]
    low <- Ks[1]
    
    pred[pred > up] <- up
    pred[pred < low]<- low
  }
  
  weight <- rep(NA, n1 + n0)
  if (estimate == "ATT") {
    weight[z==1] <- rep(1/n1, n1)
    weight[z==0] <- pred[z==0]/(1 - pred[z==0]) * 1/n1
  } else if (estimate == "ATC") {
    weight[z==1] <- (1 - pred[z==0])/pred[z==0] * 1/n0
    weight[z==0] <- rep(1/n0, n0)
  } else if (estimate == "ATE") {
    weight[z==1] <- 1/pred[z==1] * 1/n
    weight[z==0] <- 1/(1-pred[z==0]) * 1/n
  }
  return(weight)
}

calc_weight_bal <- function(data, constraint,  estimate = c("ATT", "ATC","feasible"), 
                                method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                            solver = c("cplex","gurobi"),
                                ...) {
  method <- match.arg(method)
  estimate <- match.arg(estimate)
  qp <- quadprog(data, constraint,  estimate, 
                 method,
                 ...)
  dots <- list(...)
  
  sol <- QPsolver(qp, solver = solver,...) # normalize to have closer to sum 1
  return(sol)
}

convert_ATE <- function(weight1, weight2, transport.matrix = FALSE,...) {
  list_weight <- list(weight1, weight2)
  check.vals <- sapply(list_weight, function(w) w$estimate)
  ATT.pres <- check.vals %in% "ATT"
  ATC.pres <- check.vals %in% "ATC"
  both <- (sum(c(ATT.pres, ATC.pres)) == 2)
  if(!both) stop("One set of weights must be from an ATT estimate and one must be from an ATC to combine")
  
  ATC.weight <- list_weight[[which(ATC.pres)]]
  ATT.weight <- list_weight[[which(ATT.pres)]]
  
  output <- list(w0 = NULL, w1 = NULL, gamma= NULL, estimate = NULL)
  class(output) <- "causalWeights"
  
  output$w0 <- ATT.weight$w0
  output$w1 <- ATC.weight$w1
  if(transport.matrix) {
    dots <- list(...)
    cost <- dots$cost
    p <- dots$p
    if(is.null(cost) | is.null(p)) stop("To calculate transportation matrix, 'p' and 'cost' must be specified")
    nzero_row <- output$w0>0
    nzero_col <- output$w1>0
    a <- output$w0[nzero_row]
    b <- output$w1[nzero_col]
    cost <- cost[nzero_row, nzero_col]
    transp_plan <- transport::transport(a, b, p = p, costm = cost)
    output$gamma <- matrix(0, nrow = length(output$w0), ncol = length(output$w1))
    gamma <- matrix(0, nrow = length(a), ncol = length(b))
    for(i in 1:nrow(transp_plan)) { 
      gamma[transp_plan$from[i], transp_plan$to[i]] <- transp_plan$mass[i]
    }
    output$gamma[nzero_row, nzero_col] <- gamma
  }
  output$estimate <- "ATE"
  return(output)
}

calc_gamma <- function(weights, ...) {
  dots <- list(...)
  n1 <- length(weights$w1)
  n0 <- length(weights$w0)
  if(!is.null(dots$cost) & ! is.null(dots$p)) {
    tplan <- transport::transport(a = weights$w0,
                                          b = weights$w1,
                                          p = dots$p,
                                          costm = dots$cost)
    gamma <- matrix(0, nrow = n0, ncol = n1)
    gamma[cbind(tplan[[1]], tplan[[2]])] <- tplan[[3]]
    
  } else {
    gamma <- NULL
  }
  return(gamma)
}

setOldClass("DataSim")
# setGeneric("calc_weight", function(data, ...) UseMethod("calc_weight"))
# setMethod("calc_weight", "DataSim", calc_weight.DataSim)