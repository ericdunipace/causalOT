setClass("causalWeights", slots = c(w0 = "numeric", w1="numeric", gamma = "NULL",estimand = "character"))

calc_weight <- function(data, constraint=NULL,  estimand = c("ATT", "ATC","ATE", "cATE", "feasible"), 
                        method = c("SBW", "RKHS", "RKHS.dose", "Wasserstein", "Constrained Wasserstein",
                                   "NNM",
                                   "Logistic"),
                        transport.matrix = FALSE, grid.search = FALSE,
                        ...) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  grid.search <- isTRUE(grid.search)
  args <- list(data = data, constraint = constraint, estimand = estimand, method = method, ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  output <- if (method == "SBW" & grid.search) {
    args$method <- NULL
    # if(is.null(args$grid)) 
    if(!is.null(constraint)) args$grid <- constraint
    f.call <- as.call(c(list(as.name("sbw_grid_search")), argn))
    
    
    eval(f.call, envir = args) #do.call("sbw_grid_search", args)
  } else if (method == "RKHS" | method == "RKHS.dose") {
    if (grid.search & method == "RKHS.dose") {
      f.call <- as.call(c(list(as.name("RKHS_grid_search")), argn))
      eval(f.call, envir = args)
      # do.call("RKHS_grid_search", list(data = data, estimand = estimand, 
      #                                  method = method,
      #                                  ...))
    } else {
      f.call <- as.call(c(list(as.name("calc_weight_RKHS")), argn))
      eval(f.call, envir = args)
      # do.call("calc_weight_RKHS", list(data = data, estimand = estimand, 
      #                                  method = method,
      #                                  ...))
      
    }
  } else if (method == "NNM" ) {
    f.call <- as.call(c(list(as.name("calc_weight_NNM")), argn))
    eval(f.call, envir = args)
  } else if( method != "Logistic") {
    f.call <- as.call(c(list(as.name("calc_weight_bal")), argn))
    eval(f.call, envir = args)
    # do.call("calc_weight_bal", list(data, constraint,  estimand = estimand, 
    #                                 method = method,
    #                                 ...))
  } else {
    if(estimand == "cATE") estimand <- "ATE"
    calc_weight_glm (data, constraint, estimand, ...)
  }
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, ...)
  }
  output$estimand <- estimand
  class(output) <- "causalWeights"
  return(output)
}

calc_weight_NNM <- function(data, estimand = c("ATT", "ATC","ATE", "cATE"),
                            ...) {
  
  est <- match.arg(estimand)
  pd <- prep_data(data,...)
  z <- pd$z
  df <- pd$df
  
  if(any(colnames(df) == "y")) {
    df$y <- NULL
  }
  x1 <- as.matrix(df[z==1,,drop=FALSE])
  x0 <- as.matrix(df[z==0,,drop=FALSE])
  ns <- get_n(data, ...)
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  n  <- n0 + n1
  
  dots <- list(...)
  cost <- dots$cost
  if(is.null(dots$p)) dots$p <- 2
  if(is.null(dots$dist)) dots$dist <- "Lp"
  if(dots$dist == "RKHS" & is.null(dots$rkhs.args) & is.null(dots$cost)) {
    if(is.null(dots$opt.method)) dots$opt.method <- "stan"
    temp.est <- switch(estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE"
    )
    dots$rkhs.args <- RKHS_param_opt(x=data$get_x(), 
                                     z=data$get_z(), 
                                     y = data$get_y(),
                                     metric = "mahalanobis",
                                     is.dose = dots$is.dose,
                                     opt.method = dots$opt.method,
                                     estimand = temp.est,
                                     ...)
  }
  if(est == "cATE") {
    if(is.null(cost)) {
      cost <- cost_fun(x0, x1, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                      estimand = "ATC")
      if(dots$dist == "RKHS") {
        cost <- list(cost,
                     cost_fun(x0, x1, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                              estimand = "ATT"))
      }
    }
    if(dots$dist == "RKHS") {
      w0.tab <- table(apply(cost[[1]], 2, which.min))
      w1.tab <- table(apply(cost[[2]], 1, which.min))
      w0 <- w0.tab/n1
      w1 <- w1.tab/n0
    } else {
      w0.tab <- table(apply(cost, 2, which.min))
      w1.tab <- table(apply(cost, 1, which.min))
      w0 <- w0.tab/n1
      w1 <- w1.tab/n0
    }
    
  } else if (est == "ATE") {
    x <- as.matrix(df)
    if(is.null(cost)) {
      cost <- list(cost_fun(x0, x, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                            estimand = "ATE"),
                   cost_fun(x1, x, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                            estimand = "ATE")
      )
    }
    
    w0.tab <- tabulate(apply(cost[[1]], 2, which.min), nbins = n0)
    w1.tab <- tabulate(apply(cost[[2]], 2, which.min), nbins = n1)
    w0 <- w0.tab/n
    w1 <- w1.tab/n
  } else if (est == "ATT") {
    if(is.null(cost)) {
      cost <- cost_fun(x0, x1, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                       estimand = "ATT")
    }
    
    if(dots$dist == "RKHS") {
      w0.tab <- tabulate(apply(cost[[1]], 2, which.min), nbins = n0)
      # w1.tab <- table(apply(cost[[2]], 1, which.min))
    } else {
      w0.tab <- tabulate(apply(cost, 2, which.min), nbins = n0)
      # w1.tab <- table(apply(cost, 1, which.min))
    }
    w0 <- w0.tab/n1
    w1 <- rep(1/n1,n1)
   
  } else if (est == "ATC") {
    if(is.null(cost)) {
      cost <- cost_fun(x0, x1, power = dots$p, metric = dots$dist, rkhs.args = dots$rkhs.args,
                     estimand = "ATC")
    }
    if(dots$dist == "RKHS") {
      # w0.tab <- table(apply(cost[[1]], 2, which.min))
      w1.tab <- tabulate(apply(cost[[2]], 1, which.min), nbins = n1)
    } else {
      # w0.tab <- table(apply(cost, 2, which.min))
      w1.tab <- tabulate(apply(cost, 1, which.min), nbins = n1)
    }
    w0 <- rep(1/n0, n0)
    w1 <- w1.tab/n1
  }
  
  
  output <- list(w0 = w0, w1 = w1, gamma = NULL)
  # output <- convert_sol(sol, estimand, method, ns["n0"], ns["n1"])
  
  return(output)
}

calc_weight_glm <- function(data, constraint,  estimand = c("ATT", "ATC","ATE"),
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
  if (estimand == "ATT") {
    output$w1 <- rep(1/n1,n1)
    output$w0 <- pred[z==0]/(1 - pred[z==0]) * 1/n1
  } else if (estimand == "ATC") {
    output$w1 <- (1 - pred[z==1])/pred[z==1] * 1/n0
    output$w0 <- rep(1/n0,n0)
  } else if (estimand == "ATE") {
    output$w1 <- 1/pred[z==1] * 1/n
    output$w0 <- 1/(1-pred[z==0]) * 1/n
  }
  return(output)
}

calc_weight_bal <- function(data, constraint,  estimand = c("ATT", "ATC","ATE", "cATE", "feasible"), 
                            method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                            solver = c("cplex","gurobi", "mosek"),
                            ...) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  qp <- quadprog(data, constraint,  estimand, 
                 method,
                 ...)
  dots <- list(...)
  
  sol <- lapply(qp, function(q) QPsolver(q, solver = solver, ...))
  
  # sol <- switch(estimand,
  #               "ATT"  = QPsolver(qp, solver = solver,...),
  #               "ATC"  = QPsolver(qp, solver = solver,...),
  #               "cATE" = lapply(qp, function(q) QPsolver(q, solver = solver, ...)),
  #               "ATE"  = lapply(qp, function(q) QPsolver(q, solver = solver, ...)))
  
  # sol <- if( estimand == "ATE" & (method == "Wasserstein" | method == "Constrained Wasserstein") ) {
  #   lapply(1:2, function(i) QPsolver(qp[[i]], solver = solver, ...))
  # } else {
  #   QPsolver(qp, solver = solver,...) # normalize to have closer to sum 1
  # }
  ns <- get_n(data, ...)
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  output <- convert_sol(sol, estimand, method, ns["n0"], ns["n1"])
  
  return(output)
}

calc_weight_RKHS <- function(data, estimand = c("ATC", "ATT", "cATE", "ATE"), method = c("RKHS", "RKHS.dose"),
                             solver = c("gurobi","cplex","mosek"), opt.hyperparam = TRUE,
                             ...) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  opt.hyperparam <- isTRUE(opt.hyperparam)
  solver <- match.arg(solver)
  pd <- prep_data(data,...)
  
  if(estimand == "cATE") {
    stopifnot(method == "RKHS")
    args <- list(data = data, estimand = estimand, 
                 solver = solver, opt.hyperparam = opt.hyperparam,
                 ...)
    this.call <- match.call(expand.dots = TRUE)
    args$estimand <- "ATC"
    atc.call <- eval(this.call, envir = args)
    
    args$estimand <- "ATT"
    att.call <- eval(this.call, envir = args)
    
    output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
    output$w0 <- att.call$w0
    output$w1 <- atc.call$w1
    
    output$addl.args <- list("control" = list(theta = atc.call$addl.args$theta, 
                                              gamma = atc.call$addl.args$gamma,
                                              p = atc.call$addl.args$p,
                                              sigma_2 = atc.call$addl.args$sigma_2),
                             "treated" = list(theta = att.call$addl.args$theta, 
                                              gamma = att.call$addl.args$gamma,
                                              p = att.call$addl.args$p,
                                              sigma_2 = att.call$addl.args$sigma_2)
                             )
    
    return(output)
    
  }
  
  if(opt.hyperparam) {
    # pd <- prep_data(data,...)
    opt_args <- list(x=pd$df[,-1, drop = FALSE], y = pd$df$y, z = pd$z, power = 2:3, 
                     estimand = estimand, ...)
    opt_args <- opt_args[!duplicated(names(opt_args))]
    opt_argn <- lapply(names(opt_args), as.name)
    names(opt_argn) <- names(opt_args)
    f.call <- as.call(c(list(as.name("RKHS_param_opt")), opt_argn))
    param <- eval(f.call, envir = opt_args)
    rm(opt_args)
    
    args <- c(list(data = data,
                   constraint = NULL,
                   estimand = estimand,
                   method = method),
              param,
              list(...))
  } else {
    args <- list(data = data,
                 constraint = NULL,
                 estimand = estimand,
                 method = method,
                 ...)
  }
  args <- args[!duplicated(names(args))]
  qp <- do.call("quadprog",args)
  dots <- list(...)
  
  sol <- lapply(qp, function(q) QPsolver(q, solver = solver, ...)) # normalize to have closer to sum 1
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  
  ns <- get_n(data, ...)
  
  if(estimand == "ATC") {
    output$w0 <- rep(1/ns["n0"], ns["n0"])
    output$w1 <- renormalize(sol[[1]][pd$z == 1])
  } else if (estimand == "ATT") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- rep(1/ns["n1"], ns["n1"])
  } else if (estimand == "cATE") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- renormalize(sol[[2]][pd$z == 1])
  } else if (estimand == "ATE") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- renormalize(sol[[1]][pd$z == 1])
  }
  output$addl.args <- list(theta = args$theta, 
                           gamma = args$gamma,
                           p = args$p,
                           sigma_2 = args$sigma_2)
  return(output)
}

convert_sol <- function(sol, estimand, method, n0, n1) {
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  stopifnot(is.list(sol))
  
  if ( method %in% c("Wasserstein", "Constrained Wasserstein") ) {
    if (estimand == "ATC") {
      output$gamma <- matrix(sol[[1]], n0, n1)
      output$w0 <- rep.int(1/n0, n0)
      output$w1 <- colSums(output$gamma)
    } else if (estimand == "ATT") {
      output$gamma <- matrix(sol[[1]], n0, n1)
      output$w0 <- rowSums(output$gamma)
      output$w1 <- rep.int(1/n1, n1)
    } else if (estimand == "cATE") {
      output$w0 <- rowSums(matrix(sol[[2]], n0, n1))
      output$w1 <- colSums(matrix(sol[[1]], n0, n1))
    } else if (estimand == "ATE") {
      N <- n0 + n1
      output$w0 <-  rowSums(matrix(sol[[1]], n0, N))
      output$w1 <-  rowSums(matrix(sol[[2]], n1, N)) #note both are rowSums here
    }
  } else {
    if(estimand == "ATT") {
      output$w0 <- sol[[1]]
      output$w1 <- rep.int(1/n1,n1)
    } else if (estimand == "ATC") {
      output$w0 <- rep.int(1/n0,n0)
      output$w1 <- sol[[1]]
    } else if (estimand == "feasible") {
      output$w0 <- renormalize(sol[[1]][1:n0])
      output$w1 <- renormalize(sol[[1]][n0 + 1:n1])
    } else if (estimand == "cATE") {
      output$w0 <- renormalize(sol[[1]])
      output$w1 <- renormalize(sol[[2]])
    } else if (estimand == "ATE") {
      output$w0 <- renormalize(sol[[1]][1:n0])
      output$w1 <- renormalize(sol[[1]][n0 + 1:n1])
    }
  }
  
  return(output)
}

convert_ATE <- function(weight1, weight2, transport.matrix = FALSE,...) {
  list_weight <- list(weight1, weight2)
  check.vals <- sapply(list_weight, function(w) w$estimand)
  ATT.pres <- check.vals %in% "ATT"
  ATC.pres <- check.vals %in% "ATC"
  both <- (sum(c(ATT.pres, ATC.pres)) == 2)
  if(!both) stop("One set of weights must be from an ATT estimand and one must be from an ATC to combine")
  
  ATC.weight <- list_weight[[which(ATC.pres)]]
  ATT.weight <- list_weight[[which(ATT.pres)]]
  
  output <- list(w0 = NULL, w1 = NULL, gamma= NULL, estimand = NULL)
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
  output$estimand <- "ATE"
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

get_z.DataSim <- function(data,...) {
  return(data$get_z())
}

get_z.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  return(df$z)
}

get_n.DataSim <- function(data,...) {
  return(data$get_n())
}

get_n.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z==0), n1 = sum(df$z==1))
  return(ns)
}

setOldClass("DataSim")
setGeneric("get_z", function(data, ...) UseMethod("get_z"))
setMethod("get_z", "DataSim", get_z.DataSim)
setMethod("get_z", "data.frame", get_z.data.frame)
setGeneric("get_n", function(data, ...) UseMethod("get_n"))
setMethod("get_n", "DataSim", get_n.DataSim)
setMethod("get_n", "data.frame", get_n.data.frame)