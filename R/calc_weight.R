setClass("causalWeights", slots = c(w0 = "numeric", w1="numeric", gamma = "NULL",estimate = "character"))

calc_weight <- function(data, constraint,  estimate = c("ATT", "ATC","ATE","feasible"), 
                                method = c("SBW","Wasserstein", "Constrained Wasserstein",
                                           "Logistic"),
                        transport.matrix = FALSE, grid.search = FALSE,
                                ...) {
  method <- match.arg(method)
  estimate <- match.arg(estimate)
  grid.search <- ifelse(isTRUE(grid.search), TRUE, FALSE)
  
  output <- if(method == "SBW" & grid.search) {
    do.call("sbw_grid_search", list(data, constraint,  estimate = estimate, ...))
  } else if( method != "Logistic") {
    do.call("calc_weight_bal", list(data, constraint,  estimate = estimate, 
                                    method = method,
                                    ...))
  } else {
    calc_weight_glm (data, constraint, estimate,...)
  }
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, ...)
  }
  output$estimate <- estimate
  class(output) <- "causalWeights"
  return(output)
}

calc_weight_glm <- function(data, constraint,  estimate = c("ATT", "ATC","ATE"),
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
  mod <- glm(dots$formula, data.frame(z = z, df), family = binomial(link = "logit"))
  pred <- predict(mod, type = "response")
  
  if (constraint > 0 & constraint < 1) {
    Ks  <- sort(c(constraint, 1-constraint))
    up  <- Ks[2]
    low <- Ks[1]
    
    pred[pred > up] <- up
    pred[pred < low]<- low
  }
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  if (estimate == "ATT") {
    output$w1 <- rep(1/n1,n1)
    output$w0 <- pred[z==0]/(1 - pred[z==0]) * 1/n1
  } else if (estimate == "ATC") {
    output$w1 <- (1 - pred[z==1])/pred[z==1] * 1/n0
    output$w0 <- rep(1/n0,n0)
  } else if (estimate == "ATE") {
    output$w1 <- 1/pred[z==1] * 1/n
    output$w0 <- 1/(1-pred[z==0]) * 1/n
  }
  return(output)
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
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  gamma <- NULL
  ns <- get_n(data, ...)
  
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
      output$w0 <- renormalize(sol[1:ns["n0"]])
      output$w1 <- renormalize(sol[ns["n0"] + 1:ns["n1"]])
    }
  }
  return(output)
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
    # nzero_row <- output$w0>0
    # nzero_col <- output$w1>0
    # a <- output$w0[nzero_row]
    # b <- output$w1[nzero_col]
    # cost <- cost[nzero_row, nzero_col]
    # transp_plan <- transport::transport(a, b, p = p, costm = cost)
    # output$gamma <- matrix(0, nrow = length(output$w0), ncol = length(output$w1))
    # gamma <- matrix(0, nrow = length(a), ncol = length(b))
    # for(i in 1:nrow(transp_plan)) { 
    #   gamma[transp_plan$from[i], transp_plan$to[i]] <- transp_plan$mass[i]
    # }
    # output$gamma[nzero_row, nzero_col] <- gamma
    output$gamma <- calc_gamma(weights = output, cost = cost, p = p)
  }
  output$estimate <- "ATE"
  return(output)
}

calc_gamma <- function(weights, ...) {
  dots <- list(...)
  n1 <- length(weights$w1)
  n0 <- length(weights$w0)
  if(!is.null(dots$cost) & ! is.null(dots$p)) {
    nzero_row <- weights$w0>0
    nzero_col <- weights$w1>0
    
    a <- weights$w0[nzero_row]
    b <- weights$w1[nzero_col]
    n_a <- length(a)
    n_b <- length(b)
    
    if(!isTRUE(all.equal(sum(a) ,1))) a <- renormalize(a)
    if(!isTRUE(all.equal(sum(b) ,1))) b <- renormalize(b)
    
    temp_gamma <- matrix(0, nrow = n_a, ncol = n_b)
    gamma <- matrix(0, nrow=n0,ncol = n1)
    if(n_a == 1 | n_b == 1) {
      if(n_a == 1) {
        temp_gamma[1,1:n_b] <- b
      } else if ( n_b == 1) {
        temp_gamma[1:n_a,1] <- a
      }
    } else {
      cost <- dots$cost[nzero_row, nzero_col, drop = FALSE]
      transp_plan <- transport::transport(a, b, p = p, costm = cost)
      tplan <- transport::transport(a = a,
                                    b = b,
                                    p = dots$p,
                                    costm = cost)
     
      temp_gamma[cbind(tplan[[1]], tplan[[2]])] <- tplan[[3]]
      
    }
    gamma[nzero_row, nzero_col] <- temp_gamma
    
    
  } else {
    gamma <- NULL
  }
  return(gamma)
}

get_n.DataSim <- function(data,...) {
  return(data$get_n())
}

get_n.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z==1), n1 = sum(df$z==0))
  return(ns)
}

setOldClass("DataSim")
setGeneric("get_n", function(data, ...) UseMethod("get_n"))
setMethod("get_n", "DataSim", get_n.DataSim)
setMethod("get_n", "data.frame", get_n.data.frame)