setClass("causalWeights", slots = c(w0 = "numeric", w1 = "numeric", gamma = "matrix",estimand = "character",
                                    method = "character", args = "list"))
setClass("sampleWeights", slots = c(a = "numeric", b = "numeric", total = "numeric"))

calc_weight <- function(data, constraint=NULL,  estimand = c("ATE","ATT", "ATC","cATE", "feasible"), 
                        method = c("SBW", "RKHS", "RKHS.dose", "Wasserstein", "Constrained Wasserstein",
                                   "NNM",
                                   "Logistic", "None"),
                        formula = NULL,
                        transport.matrix = FALSE, grid.search = FALSE,
                        ...) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  grid.search <- isTRUE(grid.search)
  args <- list(data = data, constraint = constraint, estimand = estimand, method = method, 
               formula = formula, ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  output <- if (grid.search & method %in% c("SBW","RKHS.dose","Constrained Wasserstein","Wasserstein")) {
    # args$method <- NULL
    # if(is.null(args$grid)) 
    if (!is.null(constraint) & is.null(args$grid)) {
      args$grid <- constraint
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
    }
    
    grid.fun <- switch(method,
                       "SBW" = "sbw_grid_search",
                       "Constrained Wasserstein" = "wass_grid_search",
                       "Wasserstein" = "wass_grid_search",
                       "RKHS.dose" = "RKHS_grid_search"
                       )
    
    f.call <- as.call(c(list(as.name(grid.fun)), argn))
    
    eval(f.call, envir = args)
    
  } else if (method == "RKHS" | method == "RKHS.dose") {
    
      f.call <- as.call(c(list(as.name("calc_weight_RKHS")), argn))
      eval(f.call, envir = args)
      # do.call("calc_weight_RKHS", list(data = data, estimand = estimand, 
      #                                  method = method,
      #                                  ...))
      
  } else if (method == "NNM" ) {
    f.call <- as.call(c(list(as.name("calc_weight_NNM")), argn))
    eval(f.call, envir = args)
  }  else if (method == "None") {
    ns <- get_n(data, ... )
    n0 <- ns[1]
    n1 <- ns[2]
    
    list(w0 = rep(1/n0, n0),
                   w1 = rep(1/n1, n1),
                   gamma = NULL,
                   estimand = estimand,
                   method = "None",
                   args = list(NULL))
  } else if ( method != "Logistic") {
    f.call <- as.call(c(list(as.name("calc_weight_bal")), argn))
    eval(f.call, envir = args)
    # do.call("calc_weight_bal", list(data, constraint,  estimand = estimand, 
    #                                 method = method,
    #                                 ...))
  } else {
    if (estimand == "cATE") estimand <- "ATE"
    calc_weight_glm(data = data, constraint = constraint, estimand = estimand, formula = formula, ...)
  }
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, ...)
  }
  output$estimand <- estimand
  output$method <- method
  class(output) <- "causalWeights"
  return(output)
}

calc_weight_NNM <- function(data, estimand = c("ATE","ATT", "ATC", "cATE"),
                            transport.matrix = FALSE, sample_weight = NULL,
                            ...) {
  
  est <- match.arg(estimand)
  pd <- prep_data(data,...)
  z <- pd$z
  df <- pd$df
  
  if (any(colnames(df) == "y")) {
    y <- df$y
    df$y <- NULL
  }
  x  <- as.matrix(df)
  # x1 <- x[z == 1,, drop = FALSE]
  # x0 <- x[z == 0,, drop = FALSE]
  ns <- get_n(data, ...)
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  n  <- n0 + n1
  
  margmass <- get_sample_weight(sample_weight, z)
  
  dots <- list(...)
  cost <- dots$cost
  p <- dots$p
  if (is.null(p)) p <- 2
  if (is.null(dots$metric)) dots$metric <- "mahalanobis"
  if (dots$metric == "RKHS" & is.null(dots$rkhs.args) & is.null(dots$cost)) {
    if (is.null(dots$opt.method)) dots$opt.method <- "stan"
    temp.est <- switch(estimand,
                       "ATT" = "ATT",
                       "ATC" = "ATC",
                       "ATE"
    )
    if (is.null(dots$kernel)) dots$kernel <- "RBF"
    
    dots$rkhs.args <- RKHS_param_opt(x = x, 
                                     z = z, 
                                     y = y,
                                     metric = dots$metric,
                                     kernel = dots$kernel,
                                     is.dose = dots$is.dose,
                                     opt.method = dots$opt.method,
                                     estimand = temp.est,
                                     ...)
  }
  if (is.null(cost)) {
    if (est == "cATE") {
      cost <- list(cost_fun(x, z, power = p, metric = dots$metric, rkhs.args = dots$rkhs.args,
                            estimand = "ATC"),
                   cost_fun(x, z, power = p, metric = dots$metric, rkhs.args = dots$rkhs.args,
                            estimand = "ATT"))
    } else {
      cost <- cost_fun(x, z, power = p, metric = dots$metric, rkhs.args = dots$rkhs.args,
                       estimand = est)
    }
    
  }
  if (est == "cATE") {
    
    # if(dots$metric == "RKHS") {
    index0 <- factor(apply(cost[[1]]^p, 2, which.min), levels = 1:n0)
    index1 <- factor(apply(cost[[2]]^p, 1, which.min), levels = 1:n1)
    w0 <- tapply(margmass$b, INDEX = index0, FUN = "sum", default = 0)
    w1 <- tapply(margmass$a, INDEX = index1, FUN = "sum", default = 0)
      # w0.tab <- tabulate(apply(cost[[1]]^p, 2, which.min), nbins = n0)
      # w1.tab <- tabulate(apply(cost[[2]]^p, 1, which.min), nbins = n1)
      # w0 <- w0.tab / n1
      # w1 <- w0.tab / n0
    # } else {
    #   w0.tab <- tabulate(apply(cost^p, 2, which.min), nbins = n0)
    #   w1.tab <- tabulate(apply(cost^p, 1, which.min), nbins = n1)
    #   w0 <- w0.tab/n1
    #   w1 <- w1.tab/n0
    # }
    
  } else if (est == "ATE") {
    # x <- as.matrix(df)
    index0 <- factor(apply(cost[[1]]^p, 2, which.min), levels = 1:n0)
    index1 <- factor(apply(cost[[2]]^p, 2, which.min), levels = 1:n1)
    w0 <- tapply(margmass$total, INDEX = index0, FUN = "sum", default = 0)
    w1 <- tapply(margmass$total, INDEX = index1, FUN = "sum", default = 0)
    # w0.tab <- tabulate(apply(cost[[1]]^p, 2, which.min), nbins = n0)
    # w1.tab <- tabulate(apply(cost[[2]]^p, 2, which.min), nbins = n1)
    # w0 <- w0.tab / n
    # w1 <- w1.tab / n
  } else if (est == "ATT") {
    
   # if (dots$metric == "RKHS") {
   #    # w0.tab <- tabulate(apply(cost[[1]]^p, 2, which.min), nbins = n0)
   #   index <- factor(apply(cost[[1]]^p, 2, which.min), levels = 1:n0)
   #   w0 <-  tapply(margmass$b, INDEX = index, FUN = "sum", default = 0)
   #  } else {
      # w0.tab <- tabulate(apply(cost^p, 2, which.min), nbins = n0)
      index <- factor(apply(cost^p, 2, which.min), levels = 1:n0)
      w0 <-  tapply(margmass$b, INDEX = index, FUN = "sum", default = 0)
    # }
    # w0 <- w0.tab * / n1
    w1 <- margmass$b #rep(1/n1,n1)
   
  } else if (est == "ATC") {
    # if (dots$metric == "RKHS") {
    #   # w1.tab <- tabulate(apply(cost[[2]]^p, 1, which.min), nbins = n1)
    #   index <- factor(apply(cost[[2]]^p, 1, which.min), levels = 1:n1)
    #   w1 <-  tapply(margmass$a, INDEX = index, FUN = "sum", default = 0)
    # } else {
      # w1.tab <- tabulate(apply(cost^p, 1, which.min), nbins = n1)
      index <- factor(apply(cost^p, 1, which.min), levels = 1:n1)
      w1 <-  tapply(margmass$a, INDEX = index, FUN = "sum", default = 0)
    # }
    w0 <- margmass$a #rep(1/n0, n0)
    # w1 <- w1.tab / n0
  }
  
  addl.args <- list(power = p, 
                    metric = dots$metric)
  if (addl.args$metric == "RKHS") {
    addl.args$rkhs.args <- dots$rkhs.args
  }
  
  output <- list(w0 = c(w0), w1 = c(w1), 
                 gamma = NULL,
                 estimand = estimand,
                 method = "NNM", 
                 args = addl.args)
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, cost = cost, p = p, ...)
  }
  # output <- convert_sol(sol, estimand, method, ns["n0"], ns["n1"])
  
  return(output)
}

calc_weight_glm <- function(data, constraint,  estimand = c("ATE","ATT", "ATC"),
                            ...) {
  dots <- list(...)
  pd <- prep_data(data,...)
  z <- pd$z
  df <- pd$df
  
  n1 <- sum(z)
  n0 <- sum(1 - z)
  n  <- n1 + n0
  if (any(colnames(df) == "y")) {
    df$y <- NULL
  }
  if ( is.null(dots$formula) ) {
    dots$formula <- formula(z ~ .)
  }
  dots$formula <- form_all_squares(dots$formula, colnames(df))
  mod <- glm(dots$formula, data.frame(z = z, df), family = binomial(link = "logit"))
  pred <- predict(mod, type = "response")
  
  if (isTRUE(constraint > 0) & isTRUE(constraint < 1)) {
    Ks  <- sort(c(constraint, 1 - constraint))
    up  <- Ks[2]
    low <- Ks[1]
    
    pred[pred > up]  <- up
    pred[pred < low] <- low
  }
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL, estimand = estimand, method = "Logistic",
                 args = c(constraint = constraint , dots))
  if (estimand == "ATT") {
    output$w1 <- rep(1/n1,n1)
    output$w0 <- pred[z == 0]/(1 - pred[z == 0]) * 1/n1
  } else if (estimand == "ATC") {
    output$w1 <- (1 - pred[z == 1])/pred[z == 1] * 1/n0
    output$w0 <- rep(1/n0,n0)
  } else if (estimand == "ATE") {
    output$w1 <- 1/pred[z == 1] * 1/n
    output$w0 <- 1/(1 - pred[z == 0]) * 1/n
  }
  
  return(output)
}

calc_weight_bal <- function(data, constraint,  estimand = c("ATE","ATT", "ATC", "cATE", "feasible"), 
                            method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                            solver = c("gurobi","mosek","cplex"),
                            sample_weight = NULL,
                            ...) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  sample_weight <- get_sample_weight(sample_weight, get_z(data, ...))
  
  qp <- quadprog(data, constraint,  estimand, 
                 method, sample_weight,
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
  output <- convert_sol(sol, estimand, method, ns["n0"], ns["n1"], sample_weight)
  output$estimand <- estimand
  output$method <- method
  if(method %in% ot.methods()) {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    if(!is.null(dots$cost)) dots$cost <- NULL
    if(dots$metric == "RKHS") {
      if(is.null(dots$rkhs.args)) dots$rkhs.args <- list()
      if(is.null(dots$rkhs.args$opt.method)) dots$rkhs.args$opt.method <- "stan"
      if(is.null(dots$rkhs.args$kernel)) dots$rkhs.args$kernel <- "RBF"
    }
    dots$metric <- dots$metric
    dots$power <- dots$p
    dots$p <- NULL
  }
  output$args <- c(list(solver = solver, constraint = constraint),
                        dots)
  
  return(output)
}

calc_weight_RKHS <- function(data, estimand = c("ATE","ATC", "ATT", "cATE"), method = c("RKHS", "RKHS.dose"),
                             kernel = c("RBF", "polynomial"),
                             solver = c("gurobi","cplex","mosek"), opt.hyperparam = TRUE,
                             sample_weight = NULL,
                             ...) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  opt.hyperparam <- isTRUE(opt.hyperparam)
  solver <- match.arg(solver)
  kernel <- match.arg(kernel)
  pd <- prep_data(data,...)
  sample_weight <- get_sample_weight(sample_weight, z = pd$z)
  
  if (estimand == "cATE") {
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
    output$estimand <- estimand
    output$method <- "RKHS"
    output$args <- list("control" = list(theta = atc.call$args$theta, 
                                              gamma = atc.call$args$gamma,
                                              p = atc.call$args$p,
                                              sigma_2 = atc.call$args$sigma_2),
                             "treated" = list(theta = att.call$args$theta, 
                                              gamma = att.call$args$gamma,
                                              p = att.call$args$p,
                                              sigma_2 = att.call$args$sigma_2)
                             )
    
    return(output)
    
  }
  
  if (opt.hyperparam) {
    # pd <- prep_data(data,...)
    opt_args <- list(x = pd$df[,-1, drop = FALSE], y = pd$df$y, z = pd$z, power = 2:3, 
                     kernel = kernel,
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
                   kernel =  kernel,
                   method = method),
              param,
              list(...))
  } else {
    args <- list(data = data,
                 constraint = NULL,
                 estimand = estimand,
                 kernel = kernel,
                 method = method,
                 ...)
  }
  args <- args[!duplicated(names(args))]
  qp <- do.call("quadprog",args)
  dots <- list(...)
  
  sol <- lapply(qp, function(q) QPsolver(q, solver = solver, ...)) # normalize to have closer to sum 1
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL,
                 estimand = estimand,
                 method = "RKHS",
                 args = list(theta = args$theta, 
                      gamma = args$gamma,
                      p = args$p,
                      sigma_2 = args$sigma_2,
                      kernel = args$kernel,
                      metric = args$metric,
                      is.dose = args$is.dose,
                      is.standardized = args$is.standardized
                 ))
  
  ns <- get_n(data, ...)
  
  if (estimand == "ATC") {
    output$w0 <- sample_weight$a
    output$w1 <- renormalize(sol[[1]][pd$z == 1])
  } else if (estimand == "ATT") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- sample_weight$b
  } else if (estimand == "cATE") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- renormalize(sol[[2]][pd$z == 1])
  } else if (estimand == "ATE") {
    output$w0 <- renormalize(sol[[1]][pd$z == 0])
    output$w1 <- renormalize(sol[[1]][pd$z == 1])
  }
  return(output)
}


calc_weight_error <- function(n0 = NULL, n1 = NULL) {
  
  if (is.null(n0) || is.null(n1)) {
    return(list(w0 = NA_real_,
                w1 = NA_real_))
  }
  return(list(w0 = rep(NA_real_, n0),
              w1 = rep(NA_real_, n1)))
}

convert_sol <- function(sol, estimand, method, n0, n1, sample_weight) {
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  stopifnot(is.list(sol))
  stopifnot(inherits(sample_weight, "sampleWeights"))
  
  if ( method %in% c("Wasserstein", "Constrained Wasserstein") ) {
    if (estimand == "ATC") {
      output$gamma <- matrix(sol[[1]], n0, n1)
      output$w0 <- sample_weight$a
      output$w1 <- colSums(output$gamma)
    } else if (estimand == "ATT") {
      output$gamma <- matrix(sol[[1]], n0, n1)
      output$w0 <- rowSums(output$gamma)
      output$w1 <- sample_weight$b
    } else if (estimand == "cATE") {
      output$w0 <- rowSums(matrix(sol[[2]], n0, n1))
      output$w1 <- colSums(matrix(sol[[1]], n0, n1))
    } else if (estimand == "ATE") {
      N <- n0 + n1
      output$w0 <-  rowSums(matrix(sol[[1]], n0, N))
      output$w1 <-  rowSums(matrix(sol[[2]], n1, N)) #note both are rowSums here
    }
  } else {
    if (estimand == "ATT") {
      output$w0 <- sol[[1]]
      output$w1 <- sample_weight$b
    } else if (estimand == "ATC") {
      output$w0 <- sample_weight$a
      output$w1 <- sol[[1]]
    } else if (estimand == "feasible") {
      output$w0 <- renormalize(sol[[1]][1:n0])
      output$w1 <- renormalize(sol[[1]][n0 + 1:n1])
    } else if (estimand == "cATE") {
      output$w0 <- renormalize(sol[[1]])
      output$w1 <- renormalize(sol[[2]])
    } else if (estimand == "ATE") {
      output$w0 <- renormalize(sol[[1]])
      output$w1 <- renormalize(sol[[2]])
    }
  }
  
  return(output)
}

convert_ATE <- function(weight1, weight2, transport.matrix = FALSE, ...) {
  list_weight <- list(weight1, weight2)
  check.vals <- sapply(list_weight, function(w) w$estimand)
  ATT.pres <- check.vals %in% "ATT"
  ATC.pres <- check.vals %in% "ATC"
  both <- (sum(c(ATT.pres, ATC.pres)) == 2)
  if(!both) stop("One set of weights must be from an ATT estimand and one must be from an ATC to combine")
  if(weight1$method != weight2$method) warning("Methods of estimating weights don't agree!")
  
  ATC.weight <- list_weight[[which(ATC.pres)]]
  ATT.weight <- list_weight[[which(ATT.pres)]]
  
  output <- list(w0 = ATT.weight$w0, w1 = ATC.weight$w1, gamma= NULL, estimand = "ATE",
                 method = ATT.weight$method, args = ATT.weight$args)
  class(output) <- "causalWeights"

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
  
  return(output)
}

calc_gamma <- function(weights, ...) {
  if (!is.null(weights$gamma)) return(weights$gamma)
  dots <- list(...)
  n1 <- length(weights$w1)
  n0 <- length(weights$w0)
  if(!is.null(dots$cost) & !is.null(dots$p)) {
    if (length(dots$cost[[1]]) > 1) return(NULL)
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
      # transp_plan <- transport::transport(a, b, p = p, costm = cost)
      tplan <- transport::transport(a = a,
                                    b = b,
                                    # p = dots$p,
                                    costm = cost^dots$p)
      
      temp_gamma[cbind(tplan[[1]], tplan[[2]])] <- tplan[[3]]
      
    }
    gamma[nzero_row, nzero_col] <- temp_gamma
    # transport::transport.default
    
  } else {
    gamma <- NULL
  }
  return(gamma)
}

get_sample_weight <- function(sw,z) {
  if (inherits(sw, "sampleWeights")) {
    return(sw)
  } else if (is.list(sw)) {
    n1 <- sum(z == 1)
    n0 <- sum(z == 0)
    n <- n0 + n1
    swtotal <- rep(NA, n)
    if (any(names(sw) == "w0") & any(names(sw) == "w1")) {
      sw <- sw[c("w0","w1")]
    }
    if (any(names(sw) == "a") & any(names(sw) == "b")) {
      sw <- sw[c("a","b")]
    }
    swtotal[z == 1] <- sw[[2]] * n1
    swtotal[z == 0] <- sw[[1]] * n0
    outmass <- list(a = sw[[1]], b = sw[[2]],
                    total = renormalize(swtotal))
    stopifnot(sapply(outmass,length) %in% c(n0, n1, n0 + n1))
    
  } else if (is.numeric(sw)) {
    outmass <- list(a = renormalize(sw[z == 0]), b = renormalize(sw[z == 1]),
                total = renormalize(sw))
  } else if (is.null(sw)) {
    n1 <- sum(z == 1)
    n0 <- sum(z == 0)
    n <- n0 + n1
    outmass <- list(a = rep(1/n0,n0), b = rep(1/n1,n1),
                    total = rep(1/n, n))
  } else {
    stop("Sample weights of unknown form. must be vector length of data or list with weights for controls and then treated in separate slots (in that order)")
  }
  class(outmass) <- "sampleWeights"
  return(outmass)
}

get_z.DataSim <- function(data,...) {
  return(data$get_z())
}

get_z.data.frame <- function(data,...) {
  
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  
  if (is.null(tx_ind)) {
    stop("must specify treatment indicator 'treatment.indicator' either by name or column number")
  }
  tx.var <- if (is.character(tx_ind)) {
    match(tx_ind, colnames(data))
  } else {
    tx_ind
  }
  
  z <- as.integer(data[[tx.var]])
  attr(z, "treatment.indicator") <- "z"
  return(z)
}

get_z.matrix <- function(data,...) {
  
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  
  if (is.null(tx_ind)) {
    stop("must specify treatment indicator 'treatment.indicator' either by name or column number")
  }
  tx.var <- if (is.character(tx_ind)) {
    match(tx_ind, colnames(data))
  } else {
    tx_ind
  }
  
  z <- as.integer(data[,tx.var])
  attr(z, "treatment.indicator") <- "z"
  return(z)
}

get_n.DataSim <- function(data,...) {
  return(data$get_n())
}

get_n.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z == 0), n1 = sum(df$z == 1))
  return(ns)
}

get_n.matrix <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z == 0), n1 = sum(df$z == 1))
  return(ns)
}

get_p.DataSim <- function(data,...) {
  return(data$get_p())
}

get_p.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  if(!is.null(df$df$y)){
    return(ncol(df$df) - 1)
  } else {
    return(ncol(df$df))
  }
}

get_p.matrix <- function(data,...) {
  df <- prep_data(data,...)
  if(!is.null(df$df$y)){
    return(ncol(df$df) - 1)
  } else {
    return(ncol(df$df))
  }
}

get_subset.DataSim <- function(data,idx, ...) {
  stop("not implemented yet")
  data$set_index(idx)
  return(data)
}

get_subset.data.frame <- function(data, idx, ...) {
  return(data[idx,,drop = FALSE])
}

clear_subset.DataSim <- function(data,idx, ...) {
  stop("not implemented yet")
  data$clear_index
  # return(data)
}

clear_subset.data.frame <- function(data, idx, ...) {
  NULL
}

setOldClass("DataSim")
setOldClass("Hainmueller")
setOldClass("Sonabend2020")

setGeneric("get_z", function(data, ...) UseMethod("get_z"))
setMethod("get_z", "DataSim", get_z.DataSim)
setMethod("get_z", "data.frame", get_z.data.frame)
setMethod("get_z", "matrix", get_z.matrix)

setGeneric("get_n", function(data, ...) UseMethod("get_n"))
setMethod("get_n", "DataSim", get_n.DataSim)
setMethod("get_n", "data.frame", get_n.data.frame)
setMethod("get_n", "matrix", get_n.matrix)

setGeneric("get_p", function(data, ...) UseMethod("get_p"))
setMethod("get_p", "DataSim", get_p.DataSim)
setMethod("get_p", "data.frame", get_p.data.frame)
setMethod("get_p", "matrix", get_p.matrix)


setGeneric("get_subset", function(data, idx, ...) UseMethod("get_subset"))
setMethod("get_subset", "DataSim", get_subset.DataSim)
setMethod("get_subset", "data.frame", get_subset.data.frame)

setGeneric("clear_subset", function(data, idx, ...) UseMethod("clear_subset"))
setMethod("clear_subset", "DataSim", clear_subset.DataSim)
setMethod("clear_subset", "data.frame", clear_subset.data.frame)
