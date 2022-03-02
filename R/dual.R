#' Function to create balance functions
#'
#' @param balance a list with slots "formula", "balance.functions",
#' "balance.mean", or "balance.sd". One of the following combinations must be provided:
#' \itemize{
#' \item "formula" and "x"
#' \item "balance.functions"
#' }
#' @param x The data matrix for the source population.
#' @param target Either: 1) a list with slots target.mean and target.sd
#' or a data matrix of the target population. If a data matrix is provided
#' then the parameter `balance` should include a formula.
#'
#' @return a list with slots "balance.functions", "balance.delta",
#' "target.mean", "target.sd".
#' @keywords internal
balance.options <- function(balance, x, target) {
  
  #check if balance variable provided
  if (is.null(balance)) {
    return(list(NULL))
  }
  
  # if formula provided, will use that, otherwise uses all variables
  # in simple linear combination, ie X_1, X_2, ..., X_d
  if (is.null(balance$formula)) {
    balance$formula <- as.formula(~. + 0)
  } else {
    balance$formula <- as.formula(balance$formula)
    environment(balance$formula) <- environment()
  }
  
  # if the balance functions not provided (functions of covariates to balance)
  # will use formula to calculate from `x`. Otherwise, uses provided functions
  if (is.null(balance$balance.functions) ) {
    balance$balance.functions <- model.matrix(balance$formula, data.frame(x))
  }
  
  # get the desired means to target. also get the scale of the target data as
  # a measure of error. 
  if (is.null(balance$target.sd) | is.null(balance$target.mean)) target_m <- model.matrix(balance$formula, data.frame(target))
  if (is.null(balance$target.mean))  balance$target.mean <- matrix(colMeans(target_m),nrow = 1)
  if (is.null(balance$target.sd)) balance$target.sd <- matrixStats::colSds(target_m)
  
  # set up final list with desired quantities
  balance$balance.functions <- as.matrix(balance$balance.functions)
  balance$balance.delta <- as.double(balance$balance.constraints)
  balance$target.mean <- as.matrix(balance$target.mean)
  balance$target.sd <- as.double(balance$target.sd)
  if(length(balance$balance.delta) == 0) balance <- list(NULL)
  return(balance)
  
}

#' options for the Wasserstein optimization problem
#'
#' @param wasserstein parts of the output list can be provided a priori
#' @param x The source data matrix. Must be provided if wasserstein$cost is null.
#' @param target The target data matrix. 
#' Must be provided if wasserstein$cost is null.
#' @param sw The sample weights of each group. not necessary if provided in
#' `wasserstein` parameter.
#'
#' @return a list with slots "power", "cost","lambda","b". These are, respectively,
#' the power of the Wasserstein distance, the cost matrix, the penalty
#' parameter for regularization and the sample weights of the target population.
#' 
#' @keywords internal
wasserstein.options <- function(wasserstein, x, target, sw) {
  
  # the output list. can be provided as fully complete
  if (missing(wasserstein) || is.null(wasserstein)) {
    # wasserstein        <- list(metric = "mahalanobis")
    # wasserstein$power  <- 2
    # z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
    # wasserstein$cost   <- cost_fun(rbind(x, target), z = z, 
    #                                estimand = "ATT", power = 2, 
    #                                metric = wasserstein$metric)
    # wasserstein$lambda  <- 1.0
    # wasserstein$b <- sw$b
    wasserstein <- list()
  }
  
  # if the metric isn't provided, make it "mahalanobis", 
  # otherwise match argument
  if (is.null(wasserstein$metric)) {
    wasserstein$metric <- "mahalanobis"
  } else {
    wasserstein$metric <- match.arg(wasserstein$metric, 
                                    c("mahalanobis", "sdLp", "Lp"))
  }
  
  # set wasserstein power to default if null
  if (is.null(wasserstein$power)) wasserstein$power <- 2
  
  # calculate the cost if not provided
  if (is.null(wasserstein$cost)) {
    z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
    wasserstein$cost <- cost_fun(rbind(x, target),
                                 z = z,
                                 estimand = "ATT", 
                                 power = wasserstein$power, 
                                 metric = wasserstein$metric)
  }
  
  # sets multiplication factor of regularization
  if (is.null(wasserstein$lambda)) wasserstein$lambda <- 1.0
  
  # sets sample weights for target
  if (is.null(wasserstein$b)) {
    if (missing(sw) || is.null(sw)) {
      sw <- list(b = rep(1, nrow(target))/nrow(target))
    }
    wasserstein$b <- sw$b
  }
  
  # confirm data types of output list
  wasserstein$metric <- as.character(wasserstein$metric)
  wasserstein$power  <- as.double(wasserstein$power)
  wasserstein$cost   <- as.matrix(wasserstein$cost)^wasserstein$power
  wasserstein$lambda  <- as.double(wasserstein$lambda)
  wasserstein$b      <- as.double(wasserstein$b)
  
  stopifnot(wasserstein$power >= 1)
  stopifnot(wasserstein$lambda >= 0)
  
  return(wasserstein)
}

#' Marginal constraints for wasserstein
#'
#' @param marg.wass a list with at minimum slots "marginal.constraints". Can
#' also include "marginal.costs".
#' @param wass output from the 
#' [wasserstein.options][wasserstein.options()] function
#' @param x The source data matrix. Optional if include "marginal.costs" in 
#' `marg.wass`.
#' @param target The target data matrix. Optional if include "marginal.costs" in 
#' `marg.wass`.
#'
#' @return a list with slots "marginal.costs", "marginal.delta"
#' @keywords internal
marginal.wass.options <- function(marg.wass, wass, x, target) {
  
  # if missing marg.wass, then return NULL
  # function will assume that we don't want marginal constraints
  if (missing(marg.wass) || is.null(marg.wass)) {
    return(list(NULL))
    marg.wass        <- list()
  }
  
  # if the marginal constraints aren't provided, 
  # will assume that we don't want marginal constraints
  if (is.null(marg.wass$marginal.constraints)) return(list(NULL))
  
  # calculate the marginal costs for each covariate
  if (is.null(marg.wass$marginal.costs)) {
    z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
    
    marg.wass$marginal.costs <- lapply(1:ncol(x), function(i) 
      cost_fun(rbind(x[,i,drop = FALSE], 
                     target[,i,drop = FALSE]),
               z = z,
               estimand = "ATT", power = wass$power, 
               metric = wass$metric))
  }
  
  # raise costs to the power in the wasserstein list
  marg.wass$marginal.costs  <- lapply(marg.wass$marginal.costs, function(cc) cc^wass$power)
  
  # set up the marginal constraints vector
  marg.wass$marginal.delta  <- as.double(marg.wass$marginal.constraints^wass$power)
  if (length(marg.wass$marginal.delta) != length(marg.wass$marginal.costs)) {
    marg.wass$marginal.delta <- rep(marg.wass$marginal.delta, 
                                    length(marg.wass$marginal.costs))
  }
  
  # make sure deltas make sense
  stopifnot(all(marg.wass$marginal.delta >= 0))
  
  return(marg.wass)
}


#' Optimization controls
#'
#' @param control provided arguments
#' @param method The optimization method
#'
#' @return A list with slots depending on the optimization method
#' @keywords internal
control.options <- function(control, method) {
  
  #sets default controls for "oem" method
  if (method == "oem") {
    if (is.null(control) | !is.list(control)) {
      control <- list(scale.factor = numeric(0),
                      maxit = 500L,
                      tol = 1e-7,
                      irls.maxit = 100L,
                      irls.tol = 0.001,
                      groups = numeric(0),
                      group.weights = NULL
      )
    }
    if (is.null(control$scale.factor))  control$scale.factor <- numeric(0)
    if (is.null(control$maxit))         control$maxit <- 500L
    if (is.null(control$tol))           control$tol <- 1e-7
    if (is.null(control$irls.maxit))    control$irls.maxit <- 100L
    if (is.null(control$irls.tol))      control$irls.tol <- 0.001
    if (is.null(control$groups))        control$groups <- numeric(0)
    if (is.null(control$group.weights)) control$group.weights <- NULL
    
    control$scale.factor <- as.numeric(control$scale.factor)
    control$maxit <- as.integer(control$maxit)
    control$tol   <- as.double(control$tol)
    control$irls.maxit <- as.integer(control$irls.maxit)
    control$irls.tol   <- as.double(control$irls.tol)
    control$groups   <- as.numeric(control$groups)
    if (!is.null(control$group.weights)) control$group.weights   <- as.numeric(control$group.weights)
    
    control.names <- c("scale.factor","maxit",
                       "tol", "irls.maxit",
                       "irls.tol", "groups",
                       "group.weights")
    
  } else if (method == "lbfgs") { #args for lbfgs
    
    if (is.null(control) | !is.list(control)) {
      control <- list(trace = 0L,
                      factr = 1e7,
                      pgtol = 0,
                      abstol = 0,
                      reltol = 0,
                      lmm = 5,
                      maxit = 1000L,
                      info = FALSE
      )
    }
    if (is.null(control$trace))         control$trace <- 0
    if (is.null(control$factr))         control$factr <- 1e7
    if (is.null(control$pgtol))         control$pgtol <- 0
    if (is.null(control$abstol))        control$abstol <- 0
    if (is.null(control$reltol))        control$reltol <- 0
    if (is.null(control$lmm))           control$lmm <- 5
    if (is.null(control$maxit))         control$maxit <- 500L
    if (is.null(control$info))          control$info <- FALSE
    
    control$trace <- as.numeric(control$trace)
    control$factr <- as.numeric(control$factr)
    control$pgtol <- as.numeric(control$abstol)
    control$abstol <- as.numeric(control$abstol)
    control$reltol <- as.numeric(control$reltol)
    control$lmm   <-   as.integer(control$lmm)
    control$maxit <- as.integer(control$maxit)
    control$info   <- isTRUE(control$info)
    
    
    control.names <- c("trace","factr",
                       "pgtol", "abstol",
                       "reltol", "lmm",
                       "maxit", "info")
    
  } else {
    stop("Method ", method, " not found")
  }
  
  return(control[names(control) %in% control.names])
}



#cpp functions for speed

{
  
  otDualL2 <-  R6::R6Class("otDualL2",
                                 public = list(
                                   obj =  function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     diff <- (fmat + gmat - private$cost)
                                     
                                     obj <- sum(private$a * f) + sum(private$b * g) - 
                                       0.5 / private$lambda * sum((diff * (diff > 0))^2)
                                     return(-obj)
                                   },
                                   obj_self =  function(vars, a, cost, lambda) {
                                     n <- length(vars)
                                     
                                     fmat <- matrix(vars, n, n)
                                     gmat <- matrix(vars, n, n, byrow = TRUE)
                                     diff <- (fmat + gmat - cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     obj <- 2 * sum(a * vars) - 
                                       0.5 / lambda * sum(diff * diff)
                                     return(-obj)
                                   },
                                   grad = function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     diff <- (fmat + gmat - private$cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     f.grad <- private$a - c(rowSums(diff) / private$lambda)
                                     
                                     g.grad <- private$b - c(colSums(diff) / private$lambda)
                                     
                                     return(-c(f.grad, g.grad))
                                   },
                                   grad_self = function(vars, a, cost, lambda) {
                                     n <- length(vars)
                                     
                                     fmat <- matrix(vars, n, n)
                                     gmat <- matrix(vars, n, n, byrow = TRUE)
                                     diff <- (fmat + gmat - cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     grad <- 2 * a - 2 * colSums(diff) / lambda
                                     
                                     return(-grad)
                                   },
                                   forward = "function",
                                   h_star_inv = function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     eta <- (fmat + gmat - private$cost)
                                     return((eta * (eta > 0))/private$lambda)
                                   },
                                   get_a = function() {
                                     return(private$a)
                                   },
                                   get_b = function() {
                                     return(private$b)
                                   },
                                   get_f = function(vars) {
                                     vars[private$fidx]
                                   },
                                   get_g = function(vars) {
                                     vars[private$gidx]
                                   },
                                   get_lambda = function() {
                                     return(private$lambda)
                                   },
                                   get_weight = function(vars) {
                                     return( round_pi( self$h_star_inv(vars), private$a, private$b ))
                                   },
                                   get_dist = function(vars) {
                                     return(
                                       sum(private$cost * self$get_weight(vars))
                                     )
                                   },
                                   update_a = "function",
                                   
                                   initialize = function(
                                     x, y,
                                     a,
                                     b,
                                     p,
                                     lambda,
                                     debias = TRUE,
                                     solver = c("lbfgs","mosek","gurobi","quadprog"),
                                     cost,
                                     ...
                                   ) {
                                     private$lambda <- lambda
                                     stopifnot(private$lambda > 0)
                                     private$debias <- isTRUE(debias)
                                     private$solver <- match.arg(solver)
                                     
                                     
                                     if ( missing(cost) || is.null(cost)) {
                                       
                                       private$cost <- cost_calc_lp(x,y,p,"rowwise")^p
                                     } else {
                                       private$cost     <- as.matrix(cost^p)
                                     }
                                     
                                     if (private$debias) {
                                       private$cost_xx <- cost_calc_lp(x,x,p, "rowwise")^p
                                       private$cost_yy <- cost_calc_lp(y,y,p, "rowwise")^p
                                     }
                                     
                                     private$n     <- nrow(private$cost)
                                     private$m     <- ncol(private$cost)
                                     private$fidx <- 1:private$n
                                     private$gidx <- (private$n + 1):(private$n + private$m)
                                     private$a  <- a
                                     private$b  <- b
                                     private$nvars <- length(a) + length(b)
                                     private$options <- control.options(control = list(...), method = "lbfgs")
                                     private$init <- rep(0, private$nvars)
                                     
                                     self$forward <- switch(private$solver,
                                                            "lbfgs" = private$forward_lbfgs,
                                                            private$forward_solver)
                                     self$update_a <- switch(private$solver,
                                                             "lbfgs" = private$update_a_lbfgs,
                                                             private$update_a_solver)
                                     
                                     if (private$solver != "lbfgs") {
                                       
                                       ### X to Y ####
                                       private$qp_xy <- private$solver_setup(private$cost, private$a, private$b)
                                       
                                       if (private$debias) {
                                         private$qp_xx <- private$solver_setup(private$cost_xx, private$a, private$a)
                                         private$qp_yy <- private$solver_setup(private$cost_yy, private$b, private$b)
                                       }
                                     }
                                   }
                                 ),
                                 private = list(
                                   a = "numeric",
                                   b = "numeric",
                                   cost    = "matrix",
                                   cost_xx = "matrix",
                                   cost_yy = "matrix",
                                   debias  = "logical",
                                   f = "numeric",
                                   fidx = "integer",
                                   forward_lbfgs = function(...) {
                                     
                                     fit <- lbfgsb3c::lbfgsb3c(par = private$init,
                                                               fn = otDualL2_obj_,
                                                               gr = otDualL2_grad_,
                                                               lower = -Inf,
                                                               upper = Inf,
                                                               control = private$options,
                                                               a = private$a, 
                                                               b = private$b,
                                                               cost = private$cost, 
                                                               lambda = private$lambda
                                     )
                                     
                                     objective <- (-fit$value)
                                     private$f <- f <- self$get_f(fit$par)
                                     private$g <- g <- self$get_g(fit$par)
                                     
                                     if (private$debias) {
                                       fit_xx <- lbfgsb3c::lbfgsb3c(par = rep(1, private$n),
                                                                    fn = otDualL2_obj_self_,
                                                                    gr = otDualL2_grad_self_,
                                                                    lower = -Inf,
                                                                    upper = Inf,
                                                                    control = private$options,
                                                                    a = private$a, 
                                                                    cost = private$cost_xx, 
                                                                    lambda = private$lambda
                                       )
                                       fit_yy <- lbfgsb3c::lbfgsb3c(par = rep(1, private$m),
                                                                    fn = otDualL2_obj_self_,
                                                                    gr = otDualL2_grad_self_,
                                                                    lower = -Inf,
                                                                    upper = Inf,
                                                                    control = private$options,
                                                                    a = private$b, 
                                                                    cost = private$cost_yy, 
                                                                    lambda = private$lambda
                                       )
                                       
                                       private$p <- fit_xx$par
                                       private$q <- fit_yy$par
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       private$loss_yy <- (-fit_yy$value)
                                       
                                       objective <- objective + 0.5 * (fit_xx$value) - 0.5 * private$loss_yy #objectives are neg due to lbfgs requirement
                                       
                                     }
                                     
                                     
                                     if(fit$convergence != 0) warning(fit$message)
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   },
                                   forward_solver = function(...) {
                                     
                                     fit <- QPsolver(private$qp_xy, solver = private$solver, get.dual = TRUE)
                                     
                                     objective <- fit$value
                                     private$f <- f <- self$get_f(fit$dual[,1] - fit$dual[,2])
                                     private$g <- g <- self$get_g(fit$dual[,1] - fit$dual[,2])
                                     
                                     if (private$debias) {
                                       fit_xx <-  QPsolver(private$qp_xx, solver = private$solver, get.dual = TRUE)
                                       fit_yy <-  QPsolver(private$qp_yy, solver = private$solver, get.dual = TRUE)
                                       
                                       private$p <- self$get_f(fit_xx$dual[,1] - fit_xx$dual[,2])
                                       private$q <- fit_yy$dual[1:private$m + private$m,1] - fit_yy$dual[1:private$m + private$m,2]
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       private$loss_yy <- fit_yy$value
                                       
                                       objective <- objective - 0.5 * fit_xx$value - 0.5 * private$loss_yy 
                                       
                                     }
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   },
                                   g = "numeric",
                                   gidx = "integer",
                                   init = "numeric",
                                   lambda = "numeric",
                                   loss_xy = "numeric",
                                   loss_xx = "numeric",
                                   loss_yy = "numeric",
                                   m = "integer",
                                   n = "integer",
                                   nvars = "integer",
                                   options = "list",
                                   p = "numeric",
                                   q = "numeric",
                                   qp_xy = "list",
                                   qp_xx = "list",
                                   qp_yy = "list",
                                   solver = "character",
                                   solver_setup = function(cost, a, b) {
                                     n0 <- length(a)
                                     n1 <- length(b)
                                     qp <- list(obj = list(L = c(as.numeric(cost))),
                                                LC = NULL)
                                     
                                     sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                                                          j = 1:(n0*n1),
                                                                          x = rep(1, n1 * n0),
                                                                          dims = c(1, n0*n1))
                                     col_const <- vec_to_col_constraints(n0,n1)
                                     row_const <- vec_to_row_constraints(n0,n1)
                                     
                                     qp$LC$A <- rbind(row_const,
                                                                 col_const,
                                                                 sum_const)
                                     qp$LC$uc <- qp$LC$lc <- c(a,
                                                                 b,
                                                                 1)
                                     qp$bounds <- list(lb = rep(0, n0*n1), ub = rep(Inf, n0*n1))
                                     soc <- if(private$solver == "mosek") {TRUE} else {FALSE}
                                     qp <- qp_pen(qp, n0, n1, a, b, "L2", private$lambda, soc = soc, divergence = FALSE)
                                     return(qp)
                                   },
                                   update_a_solver = function(a, ...) {
                                     private$a <- a
                                     
                                     private$qp_xy$LC$lc[1:private$n] <- private$qp_xy$LC$uc[1:private$n] <- a
                                     
                                     
                                     fit <- QPsolver(private$qp_xy, solver = private$solver, get.dual = TRUE)
                                     
                                     objective <- fit$value
                                     private$f <- f <- self$get_f(fit$dual[,1] - fit$dual[,2])
                                     private$g <- g <- self$get_g(fit$dual[,1] - fit$dual[,2])
                                     
                                     if (private$debias) {
                                       private$qp_xx$LC$lc[1:(2*private$n)] <- private$qp_xx$LC$uc[1:(2*private$n)] <- a
                                       fit_xx <-  QPsolver(private$qp_xx, solver = private$solver, get.dual = TRUE)
                                       
                                       private$p <- self$get_f(fit_xx$dual[,1] - fit_xx$dual[,2])
                                       
                                       
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       objective <- objective - 0.5 * (fit_xx$value) - 0.5 * private$loss_yy #objectives are neg due to lbfgs requirement
                                       
                                     }
                                     
                                     
                                     
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   },
                                   update_a_lbfgs = function(a, ...) {
                                     private$a <- a
                                     
                                     fit <- lbfgsb3c::lbfgsb3c(par = rep(0, private$n + private$m),
                                                               fn = otDualL2_obj_,
                                                               gr = otDualL2_grad_,
                                                               lower = -Inf,
                                                               upper = Inf,
                                                               control = private$options,
                                                               a = private$a, 
                                                               b = private$b,
                                                               cost = private$cost, 
                                                               lambda = private$lambda
                                                               
                                     )
                                     
                                     objective <- (-fit$value)
                                     private$f <- f <- self$get_f(fit$par)
                                     private$g <- g <- self$get_g(fit$par)
                                     if(fit$convergence != 0) warning(fit$message)
                                     
                                     if (private$debias) {
                                       fit_xx <- lbfgsb3c::lbfgsb3c(par = private$p,
                                                                    fn = otDualL2_obj_self_,
                                                                    gr = otDualL2_grad_self_,
                                                                    lower = -Inf,
                                                                    upper = Inf,
                                                                    control = private$options,
                                                                    a = private$a, 
                                                                    cost = private$cost_xx, 
                                                                    lambda = private$lambda
                                       )
                                       
                                       if(fit_xx$convergence != 0) {
                                         warning(fit_xx$message)
                                         # private$p <- rep(0, private$n)
                                         # fit_xx$value <- objective * 2
                                       } 
                                       private$p <- fit_xx$par
                                       
                                       
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       objective <- objective + 0.5 * (fit_xx$value) - 0.5 * private$loss_yy #objectives are neg due to lbfgs requirement
                                       
                                     }
                                     
                                     
                                     
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   }
                                 )
  )
  
}

# R version of L2 OT for testing
{
  
  otDualL2_RONLY <-  R6::R6Class("otDualL2_RONLY",
                                 public = list(
                                   obj =  function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     diff <- (fmat + gmat - private$cost)
                                     
                                     obj <- sum(private$a * f) + sum(private$b * g) - 
                                       0.5 / private$lambda * sum((diff * (diff > 0))^2)
                                     return(-obj)
                                   },
                                   obj_self =  function(vars, a, cost, lambda) {
                                     n <- length(vars)
                                     
                                     fmat <- matrix(vars, n, n)
                                     gmat <- matrix(vars, n, n, byrow = TRUE)
                                     diff <- (fmat + gmat - cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     obj <- 2 * sum(a * vars) - 
                                       0.5 / lambda * sum(diff * diff)
                                     return(-obj)
                                   },
                                   grad = function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     diff <- (fmat + gmat - private$cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     f.grad <- private$a - c(rowSums(diff) / private$lambda)
                                     
                                     g.grad <- private$b - c(colSums(diff) / private$lambda)
                                     
                                     return(-c(f.grad, g.grad))
                                   },
                                   grad_self = function(vars, a, cost, lambda) {
                                     n <- length(vars)
                                     
                                     fmat <- matrix(vars, n, n)
                                     gmat <- matrix(vars, n, n, byrow = TRUE)
                                     diff <- (fmat + gmat - cost)
                                     
                                     diff <- diff * (diff > 0)
                                     
                                     grad <- 2 * a - 2 * colSums(diff) / lambda
                                     
                                     return(-grad)
                                   },
                                   forward = function(...) {
                                     
                                     fit <- lbfgsb3c::lbfgsb3c(par = private$init,
                                                               fn = self$obj,
                                                               gr = self$grad,
                                                               lower = -Inf,
                                                               upper = Inf,
                                                               control = private$options
                                     )
                                     
                                     objective <- (-fit$value)
                                     private$f <- f <- self$get_f(fit$par)
                                     private$g <- g <- self$get_g(fit$par)
                                     
                                     if (private$debias) {
                                       fit_xx <- lbfgsb3c::lbfgsb3c(par = rep(1, private$n),
                                                                 fn = self$obj_self,
                                                                 gr = self$grad_self,
                                                                 lower = -Inf,
                                                                 upper = Inf,
                                                                 control = private$options,
                                                                 a = private$a, 
                                                                 cost = private$cost_xx, 
                                                                 lambda = private$lambda
                                       )
                                       fit_yy <- lbfgsb3c::lbfgsb3c(par = rep(1, private$m),
                                                                    fn = self$obj_self,
                                                                    gr = self$grad_self,
                                                                    lower = -Inf,
                                                                    upper = Inf,
                                                                    control = private$options,
                                                                    a = private$b, 
                                                                    cost = private$cost_yy, 
                                                                    lambda = private$lambda
                                       )
                                       
                                       private$p <- fit_xx$par
                                       private$q <- fit_yy$par
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       private$loss_yy <- (-fit_yy$value)
                                       
                                       objective <- objective + 0.5 * (fit_xx$value) - 0.5 * private$loss_yy #objectives are neg due to lbfgs requirement
                                       
                                     }
                                     
                                     
                                     if(fit$convergence != 0) warning(fit$message)
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   },
                                   h_star_inv = function(vars) {
                                     f <- self$get_f(vars)
                                     g <- self$get_g(vars)
                                     
                                     fmat <- matrix(f, private$n, private$m)
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     eta <- (fmat + gmat - private$cost)
                                     return((eta * (eta > 0))/private$lambda)
                                   },
                                   get_a = function() {
                                     return(private$a)
                                   },
                                   get_b = function() {
                                     return(private$b)
                                   },
                                   get_f = function(vars) {
                                     vars[private$fidx]
                                   },
                                   get_g = function(vars) {
                                     vars[private$gidx]
                                   },
                                   get_lambda = function() {
                                     return(private$lambda)
                                   },
                                   get_weight = function(vars) {
                                     return( round_pi( self$h_star_inv(vars), private$a, private$b ))
                                   },
                                   get_dist = function(vars) {
                                     return(
                                       sum(private$cost * self$get_weight(vars))
                                     )
                                   },
                                   update_a = function(a, ...) {
                                     private$a <- a
                                     
                                     fit <- lbfgsb3c::lbfgsb3c(par = c(private$f, private$g),
                                                               fn = self$obj,
                                                               gr = self$grad,
                                                               lower = -Inf,
                                                               upper = Inf,
                                                               control = private$options
                                     )
                                     
                                     objective <- (-fit$value)
                                     private$f <- f <- self$get_f(fit$par)
                                     private$g <- g <- self$get_g(fit$par)
                                     
                                     if (private$debias) {
                                       fit_xx <- lbfgsb3c::lbfgsb3c(par = private$p,
                                                                    fn = self$obj_self,
                                                                    gr = self$grad_self,
                                                                    lower = -Inf,
                                                                    upper = Inf,
                                                                    control = private$options,
                                                                    a = private$a, 
                                                                    cost = private$cost_xx, 
                                                                    lambda = private$lambda
                                       )
                                       
                                       
                                       private$p <- fit_xx$par
                                       
                                       f <- f - private$p
                                       g <- g - private$q
                                       
                                       objective <- objective + 0.5 * (fit_xx$value) - 0.5 * private$loss_yy #objectives are neg due to lbfgs requirement
                                       
                                     }
                                     
                                     
                                     if(fit$convergence != 0) warning(fit$message)
                                     
                                     return(list (f = f,
                                                  g = g,
                                                  loss = objective))
                                   },
                                   initialize = function(
                                                         x, y,
                                                         a,
                                                         b,
                                                         p,
                                                         lambda,
                                                         debias = TRUE,
                                                         cost,
                                                         ...
                                   ) {
                                     private$lambda <- lambda
                                     stopifnot(private$lambda > 0)
                                     private$debias <- isTRUE(debias)
                                     
                                     if( missing(cost) || is.null(cost)) {
                                       
                                       private$cost <- cost_calc_lp(x,y,p,"rowwise")^p
                                     } else {
                                       private$cost     <- as.matrix(cost^p)
                                     }
                                     
                                     if(debias) {
                                       private$cost_xx <- cost_calc_lp(x,x,p, "rowwise")^p
                                       private$cost_yy <- cost_calc_lp(y,y,p, "rowwise")^p
                                     }
                                     
                                     private$n     <- nrow(private$cost)
                                     private$m     <- ncol(private$cost)
                                     private$fidx <- 1:private$n
                                     private$gidx <- (private$n + 1):(private$n + private$m)
                                     private$a  <- a
                                     private$b  <- b
                                     private$nvars <- length(a) + length(b)
                                     private$options <- control.options(control = list(...), method = "lbfgs")
                                     private$init <- rep(0, private$nvars)
                                   }
                                 ),
                                 private = list(
                                   a = "numeric",
                                   b = "numeric",
                                   cost    = "matrix",
                                   cost_xx = "matrix",
                                   cost_yy = "matrix",
                                   debias  = "logical",
                                   f = "numeric",
                                   fidx = "integer",
                                   g = "numeric",
                                   gidx = "integer",
                                   init = "numeric",
                                   lambda = "numeric",
                                   loss_xy = "numeric",
                                   loss_xx = "numeric",
                                   loss_yy = "numeric",
                                   m = "integer",
                                   n = "integer",
                                   nvars = "integer",
                                   options = "list",
                                   p = "numeric",
                                   q = "numeric"
                                 )
  )
  
}



#### COT Dual R version for Testing ####
{
  
  cotDual_RONLY <-  R6::R6Class("cotDual_RONLY",
                                 public = list(
                                   obj =  function(vars) {
                                     param <- private$param_map(vars)
                                     
                                     g <- param$g
                                     xi <- param$xi
                                     
                                     eta <- private$eta_calc(g, xi, private$cost, private$lambda)
                                     
                                     obj <- sum(g * private$b) - sum(sign(xi) * private$delta) - sum(private$bal_tar * xi) - sum(private$omega(eta))
                                     
                                     
                                     return(-obj)
                                   },
                                   grad = function(vars) {
                                     
                                     param <- private$param_map(vars)
                                     xi <- param$xi
                                     g <- param$g
                                     
                                     
                                     eta <- private$eta_calc(g, xi, private$cost, private$lambda)
                                     
                                     omega.final <- private$omega_deriv(eta)
                                     
                                     xi_grad_omega <- crossprod(private$bal_mat, rowSums(omega.final))
                                     
                                     
                                     xi_grad_pos <- -sum((xi>0) * private$delta) - private$bal_targ - xi_grad_omega
                                     xi_grad_neg <- -sum((xi<0) * private$delta) + private$bal_targ + xi_grad_omega
                                     
                                     g_grad <- b - colSums(omega.final)
                                     
                                     grad <- c(g_grad, xi_grad_pos, xi_grad_neg)
                                     
                                     
                                     return(-grad)
                                   },
                                   forward = function(...) {
                                     
                                     # for (i in 1:private$options$maxit) {
                                     #   v <- self$grad(vars)
                                     #   d <- 
                                     # }
                                     
                                     fit <- lbfgsb3c::lbfgsb3c(par = private$init,
                                                               fn = self$obj,
                                                               gr = self$grad,
                                                               lower = private$bounds[,1],
                                                               upper = private$bounds[,2],
                                                               control = private$options
                                     )
                                     
                                     objective <- (-fit$value)
                                     weight <- self$h_star_inv(fit$par)
                                     a <- rowSums(weight)
                                     dual <- fit$par
                                     
                                     if (is.null(fit$convergence) || fit$convergence != 0) warning(fit$message)
                                     
                                     return(list (a = a,
                                                  gamma = weight,
                                                  dual = dual,
                                                  loss = objective))
                                   },
                                   h_star_inv = function(vars) {
                                     param <- private$param_map(vars)
                                     g <- param$g
                                     xi <- param$xi
                                     eta <- private$eta_calc(g, xi, private$cost, private$lambda)
                                     
                                     if(private$penalty == "entropy") {
                                       return(exp(eta))
                                     } else if (private$penalty == "L2") {
                                       return(eta * (eta > 0))
                                     }
                                     
                                   },
                                   get_b = function() {
                                     return(private$b)
                                   },
                                   get_g = function(vars) {
                                     vars[private$gidx]
                                   },
                                   get_lambda = function() {
                                     return(private$lambda)
                                   },
                                   get_xi = function(vars) {
                                     xipos <- vars[private$xp_idx]
                                     xineg <- vars[private$xn_idx]
                                     return(xipos - xineg)
                                   },
                                   get_dist = function(vars) {
                                     return(
                                       sum(private$cost * self$get_weight(vars))
                                     )
                                   },
                                   initialize = function(
                                     x, y,
                                     b,
                                     p,
                                     lambda,
                                     cost,
                                     balance.formula,
                                     delta,
                                     penalty = c("L2", "entropy"),
                                     ...
                                   ) {
                                     private$lambda <- lambda
                                     stopifnot(private$lambda > 0)
                                     
                                     private$penalty <- match.arg(penalty, c("L2","entropy"))
                                     
                                     if( missing(cost) || is.null(cost)) {
                                       
                                       private$cost <- cost_calc_lp(x,y,p,"rowwise")^p
                                     } else {
                                       private$cost     <- as.matrix(cost^p)
                                     }
                                     
                                     
                                     
                                     private$n     <- nrow(private$cost)
                                     private$m     <- ncol(private$cost)
                                     
                                     
                                     x <- data.frame(x)
                                     y <- data.frame(y)
                                     
                                     bform <- form_all_squares(balance.formula, colnames(x))
                                     
                                     bform.temp <- as.character(bform[length(bform)])
                                     bform <- as.formula(paste0("~ 0 +", bform.temp))
                                     
                                     targ.mat <- model.matrix(bform, y)
                                     private$bal_targ <- colMeans(targ.mat)
                                     private$bal_mat  <- model.matrix(bform, x)
                                     
                                     private$delta <- matrixStats::colSds(targ.mat) * delta
                                     
                                     private$k <- ncol(private$bal_mat)
                                     
                                     private$gidx <- 1:private$m
                                     private$xp_idx <- (private$m+1):(private$m + private$k)
                                     private$xn_idx <- (private$m + private$k + 1):(private$m + private$k * 2)
                                     
                                     
                                     private$b  <- b
                                     private$nvars <- private$m + private$k * 2
                                     
                                     private$bounds <- matrix(Inf, private$nvars, 2)
                                     private$bounds[private$gidx,1] <- -Inf
                                     private$bounds[private$xp_idx,1] <- 0
                                     private$bounds[private$xn_idx,1] <- 0
                                     
                                     private$init <- rep(0, private$nvars)
                                     
                                     private$options <- control.options(control = list(...), method = "lbfgs")
                                    
                                   }
                                 ),
                                 private = list(
                                   b = "numeric",
                                   bal_mat = "matrix",
                                   bal_targ = "numeric",
                                   bounds = "matrix",
                                   cost    = "matrix",
                                   delta = "numeric",
                                   eta_calc = function(g, xi, cost, lambda) {
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     ximat <- matrix(private$bal_mat %*% xi, private$n, private$m)
                                     eta <- (gmat - ximat - cost)/lambda
                                     return(eta)
                                   },
                                   g = "numeric",
                                   gidx = "integer",
                                   init = "numeric",
                                   k = "integer",
                                   lambda = "numeric",
                                   loss_xy = "numeric",
                                   loss_xx = "numeric",
                                   loss_yy = "numeric",
                                   m = "integer",
                                   n = "integer",
                                   nvars = "integer",
                                   omega = function(eta) { # eta already includes lambda
                                     if (private$penalty == "entropy") {
                                       return(private$lambda * exp(eta))
                                     } else if (private$penalty == "L2") {
                                       diff <- eta * (eta > 0)
                                       return(0.5 * diff * diff)
                                     }
                                   },
                                   omega_deriv = function(eta) {
                                     if (private$penalty == "entropy") {
                                       return(exp(eta))
                                     } else if (private$penalty == "L2") {
                                       diff <- eta * (eta > 0)
                                       return(diff)
                                     }
                                   },
                                   options = "list",
                                   penalty = "character",
                                   param_map = function(vars) {
                                     g <- self$get_g(vars)
                                     xi <- self$get_xi(vars)
                                     
                                     gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                                     ximat <- matrix(private$bal_mat %*% xi, private$n, private$m)
                                     
                                     eta <- gmat - ximat  - private$cost
                                     
                                     omega.out <- private$omega_deriv(eta)
                                     
                                     check <- abs(private$bal_targ + crossprod(private$bal_mat, rowSums(omega.out))) > private$delta
                                     
                                     xi  <- xi * as.numeric(check)
                                     return(list(g = g,
                                                 xi = xi))
                                   },
                                   xp_idx = "numeric",
                                   xn_idx = "numeric"
                                 )
  )
  
}


#### COT Dual ####
cotDual <- R6::R6Class("cotDual",
                       public = list(
                         obj =  cotDual_obj_,
                         grad = cotDual_grad_,
                         r_obj = function(vars) {
                           lambda <- private$lambda
                           
                           pf.obj <- sum(private$pf * vars)
                           eta <- as.numeric(private$Q %*% vars - private$cost)
                           
                           
                           omega.obj <- if(private$penalty == "L2") {
                             diff <- eta * (eta>0)
                             sum(diff * diff) * 0.5 / lambda
                           } else if (private$penalty == "entropy") {
                             log_sum_exp(eta/lambda) * lambda
                           }
                           
                           return(-(pf.obj - omega.obj))
                         },
                         r_grad = function(vars) {
                           lambda <- private$lambda
                           
                           pf.grad <- private$pf
                           eta <- as.numeric(private$Q %*% vars - private$cost)/lambda
                           
                           
                           omega.grad <- if(private$penalty == "L2") {
                             diff <- eta * (eta>0)
                             as.numeric(Matrix::crossprod(private$Q, diff))
                           } else if (private$penalty == "entropy") {
                             as.numeric(Matrix::crossprod(private$Q, exp(eta - log_sum_exp(eta))))
                           }
                           
                           return(-(pf.grad - omega.grad))
                         },
                         forward = function(init = NULL,...) {
                           
                           if (missing(init) || is.null(init)) init <- self$init()
                           
                           fit <- lbfgsb3c::lbfgsb3c(par = init,
                                                    fn = self$obj,
                                                    gr = self$grad,
                                                    lower = private$bounds[,1],
                                                    upper = private$bounds[,2],
                                                    control = private$options,
                                                    QQ = private$Q,
                                                    cost_ = private$cost,
                                                    pf_ = private$pf,
                                                    lambda = private$lambda,
                                                    penalty = private$penalty
                           )
                           
                           if(is.null(fit$convergence) || fit$convergence != 0) warning(fit$message)
                           
                           return(fit$par)
                           
                         },
                         h_star_inv = function(vars) {
                           eta <- (private$Q %*% vars - private$cost) / private$lambda
                           if(private$penalty == "entropy") {
                             return(exp(eta - log_sum_exp(eta)))
                           } else if (private$penalty == "L2") {
                             return(eta * (eta > 0))
                           }
                         },
                         get_weight = function(vars) {
                           return(matrix(self$h_star_inv(vars), private$n, private$m) )
                         },
                         get_bounds = function(vars) {
                           l_bd <- length(private$balance.delta)
                           bounds <- if (private$balfun) {
                             rbind(cbind(rep(-Inf, private$m ),
                                         rep(Inf ,  private$m)),
                                   cbind(rep(0, l_bd ),
                                         rep(Inf ,  l_bd)))
                           } else {
                             cbind(rep(-Inf, private$m ),
                                   rep(Inf ,  private$m))
                           }
                           return(bounds)
                         },
                         get_nvars = function() {
                           return(private$nvars)
                         },
                         init = function() {
                           rep(0, private$nvars)
                         },
                         penalty.factor = function() {
                           return(private$pf)
                         },
                         initialize = function(
                           b,
                           cost,
                           lambda,
                           penalty = c("L2","entropy"),
                           # marginal.costs = NULL,
                           # marginal.delta = NULL,
                           balance.functions = NULL,
                           balance.delta = NULL,
                           target.mean = NULL,
                           target.sd = NULL,
                           control = NULL,
                           ...
                         ) {
                           private$lambda <- lambda
                           private$cost <- c(cost)
                           private$n    <- nrow(cost)
                           private$m    <- ncol(cost)
                           
                           private$penalty <- match.arg(penalty)
                           private$dual.idx <- 1:private$m
                           cur.idx <- private$m + 1
                           
                           private$b <- b
                           
                           
                           
                           # if (!is.null(marginal.costs) && 
                           #     !is.null(marginal.delta) ) {
                           #   private$margins <- TRUE
                           #   private$marg.idx <- cur.idx:(cur.idx + length(marginal.costs) - 1)
                           #   cur.idx <- cur.idx + length(marginal.costs)
                           #   private$marginal.costs <- sapply(marginal.costs, c)
                           #   private$marginal.delta <- marginal.delta
                           #   stopifnot(length(private$marg.idx) == length(private$marginal.delta))
                           #   
                           # } else {
                           #   private$margins <- FALSE
                           # } 
                           
                           if (!is.null(balance.functions) & 
                               !is.null(balance.delta)) {
                             private$balfun <- TRUE
                             private$bal.idx <- cur.idx:(cur.idx + 2 * ncol(balance.functions) - 1)
                             balance.var <- matrixStats::colVars( balance.functions )
                             private$balance.functions <- Matrix::crossprod(vec_to_row_constraints(private$n, private$m),
                                                                            balance.functions)
                             # scale <- if (missing(target.sd) | is.null(target.sd)) {
                             #   sqrt(private$n/(private$n + private$m)  * matrixStats::colVars( balance.functions ) +
                             #          private$m/(private$n + private$m)  * matrixStats::colVars( target.mean ))
                             # } else {
                             #   sqrt(private$n/(private$n + private$m) * balance.var +
                             #          private$m/(private$n + private$m) * target.sd^2) 
                             # }
                             scale <- if (missing(target.sd) | is.null(target.sd)) { 
                               matrixStats::colSds( target.mean )
                             } else {
                               target.sd
                             }
                             if ( is.matrix(target.mean)) {
                               if (nrow(target.mean) > 1) target.mean <- colMeans(target.mean)
                             }
                             private$balance.delta <- c(-balance.delta * scale - target.mean, 
                                                        -balance.delta * scale + target.mean)
                             private$target.mean <- target.mean
                             stopifnot(length(private$bal.idx) == length(private$balance.delta))
                           } else {
                             private$balfun <- FALSE
                           }
                           
                           # fun.num <- if (private$margins & private$balfun) {
                           #   4L
                           # } else if (private$balfun) {
                           #   3L
                           # } else if (private$margins) {
                           #   2L
                           # } else {
                           #   1L
                           # }
                           # private$p.grad.fun = switch(fun.num,
                           #                             "1" = private$p.grad.fun.c,
                           #                             "2" = private$p.grad.fun.cm,
                           #                             "3" = private$p.grad.fun.cb,
                           #                             "4" = private$p.grad.fun.cmb)
                           # 
                           # private$p.fun = switch(fun.num,
                           #                        "1" = private$p.fun.c,
                           #                        "2" = private$p.fun.cm,
                           #                        "3" = private$p.fun.cb,
                           #                        "4" = private$p.fun.cmb)
                           # 
                           # private$Q  <- switch(fun.num,
                           #                      "1" = Matrix::t(vec_to_col_constraints_csparse(private$n,
                           #                                                                     private$m)),
                           #                      "2" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                           #                                                                           private$m)),
                           #                                  -private$marginal.costs),
                           #                      "3" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                           #                                                                           private$m)),
                           #                                  -private$balance.functions,
                           #                                  private$balance.functions),
                           #                      "4" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                           #                                                                           private$m)),
                           #                                  -private$marginal.costs,
                           #                                  -private$balance.functions,
                           #                                  private$balance.functions))
                           # private$nvars = switch(fun.num,
                           #                        "1" = private$m,
                           #                        "2" = private$m + ncol(private$marginal.costs),
                           #                        "3" = private$m + 2 * ncol(private$balance.functions),
                           #                        "4" = private$m +
                           #                          ncol(private$marginal.costs) + 
                           #                          2 * ncol(private$balance.functions))
                           # private$pf <- switch(fun.num,
                           #                      "1" = c(-private$b, private$b),
                           #                      "2" = c(-private$b, private$b,
                           #                              private$marginal.delta),
                           #                      "3" = c(-private$b, private$b,  
                           #                              private$balance.delta),
                           #                      "4" = c(-private$b, private$b,
                           #                              private$marginal.delta,
                           #                              private$balance.delta
                           #                      )
                           # )
                           
                           fun.num <- if (private$balfun) {
                             2L
                           } else {
                             1L
                           }
                           
                           private$p.grad.fun = switch(fun.num,
                                                       "1" = private$p.grad.fun.c,
                                                       "2" = private$p.grad.fun.cb)
                           
                           private$p.fun = switch(fun.num,
                                                  "1" = private$p.fun.c,
                                                  "2" = private$p.fun.cb)
                           
                           private$Q  <- switch(fun.num,
                                                "1" = Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                               private$m)),
                                                "2" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                                     private$m)),
                                                            -private$balance.functions, # this is the positive portion
                                                            private$balance.functions)) # this is the negative portion
                           private$nvars = switch(fun.num,
                                                  "1" = private$m,
                                                  "2" = private$m + 2 * ncol(private$balance.functions))
                           private$pf <- switch(fun.num,
                                                "1" = private$b,
                                                "2" = c(private$b,  
                                                        private$balance.delta)
                           )
                           
                           
                           private$bounds <- self$get_bounds()
                           
                           private$options <- control.options(control, method = "lbfgs")
                           
                         }
                       ),
                       private = list(
                         b  = "numeric",
                         bal.idx = "integer",
                         bounds = "matrix",
                         n = "integer",
                         m = "integer",
                         nvars = "integer",
                         cost = "numeric",
                         dual.idx = "integer",
                         lambda = "numeric",
                         margins = "logical",
                         balfun  = "logical",
                         marginal.costs = "list",
                         marginal.delta = "numeric",
                         balance.delta = "numeric",
                         balance.functions = "matrix",
                         Q = "matrix",
                         pf = "numeric",
                         target.mean = "numeric",
                         p.grad.fun = "function",
                         p.fun = "function",
                         p.grad.fun.c = function(vars) {
                           return(private$b)
                         },
                         p.grad.fun.cb = function(vars) {
                           beta_b <- vars[private$bal.idx]
                           return( c(private$b,
                                     -private$balance.delta )
                           )
                         }, 
                         p.fun.c = function(vars) {
                           return(sum(private$b * vars))
                         },
                         p.fun.cb = function(vars) {
                           g <- vars[private$dual.idx]
                           beta_b <- vars[private$bal.idx]
                           return( sum(private$b * g) -
                                     sum(beta_b * private$balance.delta)
                           )
                         },
                         penalty = "character",
                         options = "list"
                       )
)






#### SBW Dual ####
sbwDual <- R6::R6Class("sbwDual",
                       public = list(
                         obj =  sbwDual_obj_,
                         grad = sbwDual_grad_,
                         forward = function(init = NULL,...) {
                           
                           if (missing(init) || is.null(init)) init <- self$init()
                           
                           fit <- lbfgsb3c::lbfgsb3c(par = init,
                                                    fn = self$obj,
                                                    gr = self$grad,
                                                    lower = private$bounds[,1],
                                                    upper = private$bounds[,2],
                                                    control = private$options,
                                                    QQ = private$balance.functions,
                                                    pf_ = private$delta,
                                                    penalty = private$penalty
                           )
                           
                           if(is.null(fit$convergence) || fit$convergence != 0) warning(fit$message)
                           
                           return(fit$par)
                           
                         },
                         # obj = function(vars) {
                         #   return(sum(private$delta * private$scale * abs(vars)) +
                         #            sum(self$h_star_inv(vars)^2) + 
                         #            sum(private$target.functions * vars))
                         # },
                         # grad = function(vars) {
                         #   lambda <- private$delta * private$scale
                         #   # vars <- vars * as.numeric(abs(vars) > lambda)
                         #   return(sign(vars) * lambda +
                         #            c(private$target.mean) -
                         #            c(crossprod(private$balance.functions, self$h_star_inv(vars))))
                         # },
                         h_star_inv = function(vars) {
                           if(private$penalty == "L2") {
                             eta <- -0.5 * private$balance.functions %*% vars + 1.0/private$n
                             return(eta * as.numeric(eta > 0))
                           } else if (private$penalty == "entropy") {
                             eta <- -private$balance.functions %*% vars + log(1/private$n)
                             return(exp(eta - log_sum_exp(eta )))
                           }
                           
                         },
                         get_weight = function(vars) {
                           return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                         },
                         get_bounds = function(vars) {
                           bounds <- cbind(rep(0,    private$nvars ),
                                           rep(Inf , private$nvars ))
                           return(bounds)
                         },
                         get_nvars = function() {
                           return(private$nvars)
                         },
                         init = function() {
                           rep(0, private$nvars)
                         },
                         penalty.factor = function() {
                           return(private$pf)
                         },
                         initialize = function(
                           penalty = c("L2","entropy"),
                           balance.functions = NULL,
                           balance.delta = NULL,
                           target.mean = NULL,
                           target.sd = NULL,
                           control = NULL,
                           ...
                         ) {
                           
                           private$penalty <- match.arg(penalty)
                           private$balance.functions     <- as.matrix(balance.functions)
                           
                           private$balance.functions <- cbind(private$balance.functions, #positive vars
                                                             -private$balance.functions) #negative vars
                           private$n     <- nrow(private$balance.functions)
                           
                           if(missing(target.sd) || is.null(target.sd)) target.sd <- matrixStats::colSds(target.mean)
                           private$target.mean <- colMeans(as.matrix(target.mean))
                           
                           private$delta <-  c(-target.sd * balance.delta - private$target.mean, #positive vars
                                               -target.sd * balance.delta + private$target.mean) #negative vars
                           
                           private$nvars <- ncol(private$balance.functions)
                           
                           private$bounds <- self$get_bounds()
                           
                           private$options <- control.options(control, method = "lbfgs")
                           
                         }
                       ),
                       private = list(
                         balance.delta = "numeric",
                         balance.functions = "matrix",
                         bounds = "matrix",
                         delta = "numeric",
                         n = "integer",
                         nvars = "integer",
                         options = "list",
                         penalty = "character",
                         target.mean = "numeric"
                       )
)





#### Dual prep and run function ####

dual_opt <- function(x, target, 
                     init = NULL,
                     sample_weights = NULL, 
                     method = c("Wasserstein", "SBW"),
                     penalty = c("L2","entropy"),
                     wasserstein = list(metric = c("mahalanobis"),
                                        power = 2,
                                        cost = NULL,
                                        lambda = 1.0),
                     balance = list(balance.functions = NULL,
                                    formula = NULL,
                                    balance.constraints = NULL),
                     control = list(maxit = 2e3),
                     ...
) {
  
  method <- match.arg(method)
  penalty <- match.arg(penalty)
  
  sw  <- get_sample_weight(sample_weights, c(rep(0, nrow(x)), rep(1, nrow(target))))
  if (method == "Wasserstein" ) {
    wasserstein <- wasserstein.options(wasserstein, x, target, sw)
  } else {
    wasserstein <- list()
  }
  
  balance <- balance.options(balance, x, target)
  
  if (method == "SBW" & is.null(balance$balance.functions)) stop("Balance functions must be specified for SBW")
  if (method == "Wasserstein" & is.null(wasserstein$cost)) stop("Wasserstein cost must be provided or data to calculate the cost, with given metric and power")
  
  opt.class <- switch(method,
                      "Wasserstein" = cotDual,
                      "SBW" = sbwDual)
  
  optimizer <- opt.class$new(b = wasserstein$b,
                             lambda = wasserstein$lambda, 
                             cost = wasserstein$cost, 
                             balance.functions = balance$balance.functions,
                             balance.delta = balance$balance.delta,
                             target.mean = balance$target.mean,
                             target.sd = balance$target.sd,
                             penalty = penalty,
                             control = control)
  
  beta <- optimizer$forward(init)
  
  
  weight <- optimizer$get_weight(beta)
  
  if (grepl("Wasserstein", method)) {
    n <- nrow(wasserstein$cost)
    m <- ncol(wasserstein$cost)
    
    # unconst_weight <- Matrix::Matrix(data = unconst_weight,
    #                                  nrow = n, ncol = m)
    weight <- if(is.matrix(weight)) {
      Matrix::Matrix(data = weight)
    } else {
      Matrix::Matrix(data = weight,
                     nrow = n, ncol = m)
    }
    
  } 
  
  return(list(weight = weight, dual = beta,
              optimizer = optimizer))
}

#### Utility functions ####

#' Optimization controls
#'
#' @param control provided arguments
#' @param method The optimization method
#'
#' @return A list with slots depending on the optimization method
#' @keywords internal
control.options <- function(control, method) {
  
  #sets default controls for "oem" method
  if (method == "oem") {
    if (is.null(control) | !is.list(control)) {
      control <- list(scale.factor = numeric(0),
                      maxit = 500L,
                      tol = 1e-7,
                      irls.maxit = 100L,
                      irls.tol = 0.001,
                      groups = numeric(0),
                      group.weights = NULL
      )
    }
    if (is.null(control$scale.factor))  control$scale.factor <- numeric(0)
    if (is.null(control$maxit))         control$maxit <- 500L
    if (is.null(control$tol))           control$tol <- 1e-7
    if (is.null(control$irls.maxit))    control$irls.maxit <- 100L
    if (is.null(control$irls.tol))      control$irls.tol <- 0.001
    if (is.null(control$groups))        control$groups <- numeric(0)
    if (is.null(control$group.weights)) control$group.weights <- NULL
    
    control$scale.factor <- as.numeric(control$scale.factor)
    control$maxit <- as.integer(control$maxit)
    control$tol   <- as.double(control$tol)
    control$irls.maxit <- as.integer(control$irls.maxit)
    control$irls.tol   <- as.double(control$irls.tol)
    control$groups   <- as.numeric(control$groups)
    if (!is.null(control$group.weights)) control$group.weights   <- as.numeric(control$group.weights)
    
    control.names <- c("scale.factor","maxit",
                       "tol", "irls.maxit",
                       "irls.tol", "groups",
                       "group.weights")
    
  } else if (method == "lbfgs") { #args for lbfgs
    
    if (missing(control) || is.null(control) || !is.list(control) ) {
      control <- list(trace = 0L,
                      factr = 1e7,
                      pgtol = 0,
                      abstol = 0,
                      reltol = 0,
                      lmm = 5,
                      maxit = 1000L,
                      info = FALSE
      )
    }
    if (is.null(control$trace))         control$trace <- 0
    if (is.null(control$factr))         control$factr <- 1e7
    if (is.null(control$pgtol))         control$pgtol <- 0
    if (is.null(control$abstol))        control$abstol <- 0
    if (is.null(control$reltol))        control$reltol <- 0
    if (is.null(control$lmm))           control$lmm <- 5
    if (is.null(control$maxit))         control$maxit <- 1000L
    if (is.null(control$info))          control$info <- FALSE
    
    control$trace <- as.numeric(control$trace)
    control$factr <- as.numeric(control$factr)
    control$pgtol <- as.numeric(control$abstol)
    control$abstol <- as.numeric(control$abstol)
    control$reltol <- as.numeric(control$reltol)
    control$lmm   <-   as.integer(control$lmm)
    control$maxit <- as.integer(control$maxit)
    control$info   <- isTRUE(control$info)
    
    
    control.names <- c("trace","factr",
                       "pgtol", "abstol",
                       "reltol", "lmm",
                       "maxit", "info")
    
  } else {
    stop("Method ", method, " not found")
  }
  
  return(control[names(control) %in% control.names])
}

#' Function to create balance functions
#'
#' @param balance a list with slots "formula", "balance.functions",
#' "balance.mean", or "balance.sd". One of the following combinations must be provided:
#' \itemize{
#' \item "formula" and "x"
#' \item "balance.functions"
#' }
#' @param x The data matrix for the source population.
#' @param target Either: 1) a list with slots target.mean and target.sd
#' or a data matrix of the target population. If a data matrix is provided
#' then the parameter `balance` should include a formula.
#'
#' @return a list with slots "balance.functions", "balance.delta",
#' "target.mean", "target.sd".
#' @keywords internal
balance.options <- function(balance, x, target) {
  
  #check if balance variable provided
  if (is.null(balance)) {
    return(list(NULL))
  }
  
  if(is.list(balance)) {
    if((is.null(balance$formula) || is.na(balance$formula)) && 
       (is.null(balance$balance.functions) || is.na(balance$balance.functions) ) &&
       (is.null(balance$balance.constraints) || is.na(balance$balance.constraints))) return(balance)
  }
  
  # if formula provided, will use that, otherwise uses all variables
  # in simple linear combination, ie X_1, X_2, ..., X_d
  if (is.null(balance$formula) || is.na(balance$formula)) {
    return(list(NULL)) #balance$formula <- as.formula(~. + 0)
  } else {
    form <- form_all_squares(balance$formula, colnames(x))
    form.temp <- as.character(form[length(form)])
    balance$formula <- as.formula(paste0("~ 0 +", form.temp))
    environment(balance$formula) <- environment()
  }
  
  # if the balance functions not provided (functions of covariates to balance)
  # will use formula to calculate from `x`. Otherwise, uses provided functions
  if (is.null(balance$balance.functions) || is.na(balance$balance$functions)) {
    balance$balance.functions <- model.matrix(balance$formula, data.frame(x))
  }
  
  # get the desired means to target. also get the scale of the target data as
  # a measure of error. 
  if (is.null(balance$target.sd) | is.null(balance$target.mean)) target_m <- model.matrix(balance$formula, data.frame(target))
  if (is.null(balance$target.mean))  balance$target.mean <- matrix(colMeans(target_m),nrow = 1)
  if (is.null(balance$target.sd)) balance$target.sd <- matrixStats::colSds(target_m)
  
  # set up final list with desired quantities
  balance$balance.functions <- as.matrix(balance$balance.functions)
  balance$balance.delta <- as.double(balance$balance.constraints)
  balance$target.mean <- as.matrix(balance$target.mean)
  balance$target.sd <- as.double(balance$target.sd)
  if(length(balance$balance.delta) == 0) {
    warning("Balance formula given but not a constraint. Not using balancing functions!")
    balance <- list(NULL)
  }
  return(balance)
  
}

#' options for the Wasserstein optimization problem
#'
#' @param method If method is "Wasserstein", will fill out full list. Otherwise just returns generic parts.
#' @param wasserstein parts of the output list can be provided a priori
#' @param x The source data matrix. Must be provided if wasserstein$cost is null.
#' @param target The target data matrix. 
#' Must be provided if wasserstein$cost is null.
#' @param sw The sample weights of each group. not necessary if provided in
#' `wasserstein` parameter.
#'
#' @return a list with slots "power", "cost","lambda","b". These are, respectively,
#' the power of the Wasserstein distance, the cost matrix, the penalty
#' parameter for regularization and the sample weights of the target population.
#' 
#' @keywords internal
wasserstein.options <- function( wasserstein, x, target, sw) {
  
  # the output list. can be provided as fully complete
  if (missing(wasserstein) || is.null(wasserstein)) {
    # wasserstein        <- list(metric = "mahalanobis")
    # wasserstein$power  <- 2
    # z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
    # wasserstein$cost   <- cost_fun(rbind(x, target), z = z, 
    #                                estimand = "ATT", power = 2, 
    #                                metric = wasserstein$metric)
    # wasserstein$lambda  <- 1.0
    # wasserstein$b <- sw$b
    wasserstein <- list()
  }
  
    # if the metric isn't provided, make it "mahalanobis", 
    # otherwise match argument
    if (is.null(wasserstein$metric)) {
      wasserstein$metric <- "mahalanobis"
    } else {
      wasserstein$metric <- match.arg(wasserstein$metric, 
                                      c("mahalanobis", "sdLp", "Lp"))
    }
    
    # set wasserstein power to default if null
    if (is.null(wasserstein$power)) wasserstein$power <- 2
    
    # calculate the cost if not provided
    if (is.null(wasserstein$cost)) {
      z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
      wasserstein$cost <- cost_fun(rbind(x, target),
                                   z = z,
                                   estimand = "ATT", 
                                   power = wasserstein$power, 
                                   metric = wasserstein$metric)
    }
    
  # sets multiplication factor of regularization
  if (is.null(wasserstein$lambda)) wasserstein$lambda <- 1.0
  
  # sets sample weights for target
  if (is.null(wasserstein$b)) {
    if (missing(sw) || is.null(sw)) {
      sw <- list(b = rep(1, nrow(target))/nrow(target))
    }
    wasserstein$b <- sw$b
  }
  
  # confirm data types of output list
  wasserstein$metric <- as.character(wasserstein$metric)
  wasserstein$power  <- as.double(wasserstein$power)
  wasserstein$cost   <- as.matrix(wasserstein$cost)^wasserstein$power
  wasserstein$lambda  <- as.double(wasserstein$lambda)
  wasserstein$b      <- as.double(wasserstein$b)
  
  stopifnot(wasserstein$power >= 1)
  stopifnot(wasserstein$lambda > 0)
  
  return(wasserstein)
}



#### backup cotDuall2 ####
cotDualL2_2 <- R6::R6Class("cotDualL2",
                           public = list(
                             obj =  function(vars) {
                               
                               return(
                                 cotDualL2_obj_(vars, QQ = private$Q,
                                                cost_ = private$cost, 
                                                pf_ = as.double(private$p.fun(vars)),
                                                lambda = as.double(private$lambda)))
                               
                               
                               diff <- (private$Q %*% vars - private$cost)
                               diff <- diff * (diff > 0)
                               return(-(private$p.fun(vars) - 0.5 * sum(diff^2) / private$lambda))
                             },
                             grad = function(vars) {
                               
                               return(
                                 cotDualL2_grad_(vars_ = vars, 
                                                 QQ = private$Q,
                                                 cost_ = private$cost, 
                                                 pf_ = private$p.grad.fun(vars),
                                                 lambda = private$lambda)
                               )
                               
                               diff <- (private$Q %*% vars - private$cost)
                               diff <- diff * (diff > 0)
                               return(as.numeric(-(private$p.grad.fun(vars) - 
                                                     Matrix::crossprod(private$Q, 
                                                                       diff) / private$lambda)))
                             },
                             h_star_inv = function(vars) {
                               eta <- (private$Q %*% vars - private$cost) / private$lambda
                               return(eta * (eta > 0))
                             },
                             get_weight = function(vars) {
                               return(matrix(self$h_star_inv(vars), private$n, private$m) )
                             },
                             bounds = function(vars) {
                               l_md <- length(private$marginal.delta)
                               l_bd <- length(private$balance.delta)
                               bounds <- if (private$margins & private$balfun) {
                                 rbind(
                                   cbind(rep(-Inf, private$m ),
                                         rep(Inf ,  private$m)),
                                   cbind(rep(0, l_md),
                                         rep(Inf, l_md)),
                                   cbind(rep(0, l_bd ),
                                         rep(Inf ,  l_bd))
                                 )
                               } else if (private$margins & !private$balfun) {
                                 rbind(
                                   cbind(rep(-Inf, private$m ),
                                         rep(Inf ,  private$m)),
                                   cbind(rep(0, l_md),
                                         rep(Inf, l_md)))
                                 
                               } else if (!private$margins & private$balfun) {
                                 rbind(cbind(rep(-Inf, private$m ),
                                            rep(Inf ,  private$m)),
                                       cbind(rep(0, l_bd ),
                                             rep(Inf ,  l_bd)))
                               } else {
                                 cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m))
                               }
                               return(bounds)
                             },
                             get_nvars = function() {
                               return(private$nvars)
                             },
                             init = function() {
                               rep(max(private$cost), private$nvars)
                             },
                             get_xtx = function(vars) {
                               # diff <- (private$Q %*% vars - private$cost)
                               # pos <- (diff > 0)
                               return(Matrix::crossprod(private$Q))
                             },
                             get_xty = function() {
                               NULL
                             },
                             get_hessian = function(vars) {
                               # diff <- (private$Q %*% vars - private$cost)
                               # pos <- as(diff > 0, "sparseVector")
                               return(self$get_xtx(vars)/(-private$lambda))
                             },
                             get_grad_var = function(vars) {
                               # ignore negative sign and penalty functions since variance of constants is 0 and negative sign is squared
                               diff <- (private$Q %*% vars - private$cost)
                               # diff  <- diff * (diff > 0)
                               # pos.prob <- mean(diff > 0)
                               return(Matrix::crossprod(private$Q * diff)/private$lambda)
                             },
                             get_pos_prob = function(vars) {
                               # ignore negative sign and penalty functions since variance of constants is 0 and negative sign is squared
                               diff <- (private$Q %*% vars - private$cost)
                               # diff  <- diff * (diff > 0)
                               pos.prob <- mean(as.numeric(diff) > 0)
                               return(pos.prob)
                             },
                             penalty.factor = function() {
                               return(private$pf)
                             },
                             initialize = function(lambda, cost, 
                                                   b,
                                                   marginal.costs = NULL,
                                                   marginal.delta = NULL,
                                                   balance.functions = NULL,
                                                   balance.delta = NULL,
                                                   target.mean = NULL,
                                                   target.sd = NULL,
                                                   ...
                             ) {
                               private$lambda <- lambda
                               private$cost <- c(cost)
                               private$n    <- nrow(cost)
                               private$m    <- ncol(cost)
                               private$dual.idx <- 1:private$m
                               cur.idx <- private$m + 1
                               
                               private$b <- b
                               
                               if (!is.null(marginal.costs) && 
                                   !is.null(marginal.delta) ) {
                                 private$margins <- TRUE
                                 private$marg.idx <- cur.idx:(cur.idx + length(marginal.costs) - 1)
                                 cur.idx <- cur.idx + length(marginal.costs)
                                 private$marginal.costs <- sapply(marginal.costs, c)
                                 private$marginal.delta <- marginal.delta
                                 stopifnot(length(private$marg.idx) == length(private$marginal.delta))
                                 
                               } else {
                                 private$margins <- FALSE
                               } 
                               
                               if (!is.null(balance.functions) & 
                                   !is.null(balance.delta)) {
                                 private$balfun <- TRUE
                                 private$bal.idx <- cur.idx:(cur.idx + 2 * ncol(balance.functions) - 1)
                                 balance.var <- matrixStats::colVars( balance.functions )
                                 private$balance.functions <- Matrix::crossprod(vec_to_row_constraints(private$n, private$m),
                                                                                balance.functions)
                                 scale <- if (missing(target.sd) | is.null(target.sd)) {
                                   sqrt(private$n/(private$n + private$m)  * matrixStats::colVars( balance.functions ) +
                                          private$m/(private$n + private$m)  * matrixStats::colVars( target.mean ))
                                 } else {
                                   sqrt(private$n/(private$n + private$m) * balance.var +
                                          private$m/(private$n + private$m) * target.sd^2) 
                                 }
                                 if ( is.matrix(target.mean)) {
                                   if (nrow(target.mean) > 1) target.mean <- colMeans(target.mean)
                                 }
                                 private$balance.delta <- c(balance.delta * scale + target.mean, 
                                                            -balance.delta * scale - target.mean)
                                 private$target.mean <- target.mean
                                 stopifnot(length(private$bal.idx) == length(private$balance.delta))
                               } else {
                                 private$balfun <- FALSE
                               }
                               
                               fun.num <- if (private$margins & private$balfun) {
                                 4L
                               } else if (private$balfun) {
                                 3L
                               } else if (private$margins) {
                                 2L
                               } else {
                                 1L
                               }
                               
                               private$p.grad.fun = switch(fun.num,
                                                           "1" = private$p.grad.fun.c,
                                                           "2" = private$p.grad.fun.cm,
                                                           "3" = private$p.grad.fun.cb,
                                                           "4" = private$p.grad.fun.cmb)
                               
                               private$p.fun = switch(fun.num,
                                                      "1" = private$p.fun.c,
                                                      "2" = private$p.fun.cm,
                                                      "3" = private$p.fun.cb,
                                                      "4" = private$p.fun.cmb)
                               
                               private$Q  <- switch(fun.num,
                                                    "1" = Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                                   private$m)),
                                                    "2" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                                         private$m)),
                                                                -private$marginal.costs),
                                                    "3" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                                         private$m)),
                                                                -private$balance.functions,
                                                                private$balance.functions),
                                                    "4" = cbind(Matrix::t(vec_to_col_constraints_csparse(private$n,
                                                                                                         private$m)),
                                                                -private$marginal.costs,
                                                                -private$balance.functions,
                                                                private$balance.functions))
                               private$nvars = switch(fun.num,
                                                      "1" = private$m,
                                                      "2" = private$m + ncol(private$marginal.costs),
                                                      "3" = private$m + 2 * ncol(private$balance.functions),
                                                      "4" = private$m +
                                                        ncol(private$marginal.costs) + 
                                                        2 * ncol(private$balance.functions))
                               private$pf <- switch(fun.num,
                                                    "1" = c(-private$b, private$b),
                                                    "2" = c(-private$b, private$b,
                                                            private$marginal.delta),
                                                    "3" = c(-private$b, private$b,  
                                                            private$balance.delta),
                                                    "4" = c(-private$b, private$b,
                                                            private$marginal.delta,
                                                            private$balance.delta
                                                    )
                               )
                               
                             }
                           ),
                           private = list(
                             b  = "numeric",
                             n = "integer",
                             m = "integer",
                             nvars = "integer",
                             cost = "numeric",
                             lambda = "numeric",
                             margins = "logical",
                             balfun  = "logical",
                             marginal.costs = "list",
                             marginal.delta = "numeric",
                             balance.delta = "numeric",
                             balance.functions = "matrix",
                             Q = "matrix",
                             pf = "numeric",
                             dual.idx = "integer",
                             marg.idx = "integer",
                             bal.idx = "integer",
                             target.mean = "numeric",
                             p.grad.fun = "function",
                             p.fun = "function",
                             p.grad.fun.c = function(vars) {
                               return(private$b)
                             },
                             p.grad.fun.cm = function(vars) {
                               return( c(private$b, 
                                         -private$marginal.delta)
                               )
                             },
                             p.grad.fun.cb = function(vars) {
                               beta_b <- vars[private$bal.idx]
                               return( c(private$b,
                                         -private$balance.delta )
                               )
                             }, 
                             p.grad.fun.cmb = function(vars) {
                               beta_b <- vars[private$bal.idx]
                               return( c(private$b,
                                         -private$marginal.delta,
                                         -private$balance.delta)
                               )
                             },
                             p.fun.c = function(vars) {
                               return(sum(private$b * vars))
                             },
                             p.fun.cm = function(vars) {
                               g <- vars[private$dual.idx]
                               beta_m <- vars[private$marg.idx]
                               return( sum(private$b * g) -
                                         sum(private$marginal.delta * beta_m)
                               )
                             },
                             p.fun.cb = function(vars) {
                               g <- vars[private$dual.idx]
                               beta_b <- vars[private$bal.idx]
                               return( sum(private$b * g) -
                                         sum(beta_b * private$balance.delta)
                               )
                             }, 
                             p.fun.cmb = function(vars) {
                               g <- vars[private$dual.idx]
                               beta_m <- vars[private$marg.idx]
                               beta_b <- vars[private$bal.idx]
                               return( sum(private$b * g) -
                                         sum(private$marginal.delta * beta_m) -
                                         sum(beta_b * private$balance.delta)
                               )
                             }
                             # , h_star_inv.c = function(vars) {
                             #   eta <- rep(private$vars, each = private$n) - private$cost
                             #   return(eta * as.numeric(eta > 0) / private$lambda)
                             # },
                             # h_star_inv.cm = function(vars) {
                             #   g <- vars[private$dual.idx]
                             #   beta_m <- vars[private$marg.idx]
                             #   eta <- rep(g, each = private$n) +  
                             #     private$marginal.costs %*% beta_m - 
                             #     private$cost
                             #   return(eta * as.numeric(eta > 0) / private$delta)
                             # },
                             # h_star_inv.cb = function(vars) {
                             #   g <- vars[private$dual.idx]
                             #   beta_m <- vars[private$marg.idx]
                             #   eta <- rep(g, each = private$n) +  
                             #     private$marginal.costs %*% beta_m - 
                             #     private$cost
                             #   return(eta * as.numeric(eta > 0) / private$delta)
                             # },
                             # h_star_inv.cmb = function(vars) {
                             #   g <- vars[private$dual.idx]
                             #   beta_m <- vars[private$marg.idx]
                             #   beta_b <- vars[private$bal.idx]
                             #   eta <- rep(g, each = private$n) +  
                             #     private$marginal.costs %*% beta_m +
                             #     Matrix::crossprod(private$balance.functions, beta_b) - 
                             #     private$cost
                             #   return(eta * as.numeric(eta > 0) / private$delta)
                             # }
                           )
)

#### cot_dual opt ####

cot_dual_opt <- function(x, target, 
                         init = NULL,
                         sample_weights = NULL, 
                         method = c("SBW", "Wasserstein"),
                         penalty = c("KL","Entropy","L2"),
                         wasserstein = list(metric = c("mahalanobis"),
                                            power = 2,
                                            cost = NULL,
                                            lambda = 1.0),
                         balance = list(balance.functions = NULL,
                                        formula = NULL,
                                        balance.constraints = NULL),
                         marginal.wasserstein = list(marginal.costs = NULL,
                                                     marginal.constraints = NULL),
                         control = list(maxit = 1e4)
) {
  
  method <- match.arg(method)
  
  sw  <- get_sample_weight(sample_weights, c(rep(0, nrow(x)), rep(1, nrow(target))))
  if (method == "Wasserstein" | method == "Constrained Wasserstein") {
    wasserstein <- wasserstein.options(wasserstein, x, target, sw)
    marginal.wasserstein <- marginal.wass.options(marginal.wasserstein,
                                                  wasserstein, x, target)
  }
  
  balance <- balance.options(balance, x, target)
  
  if (method == "SBW" & is.null(balance$balance.functions)) stop("Balance functions must be specified for SBW")
  if (method == "Wasserstein" & is.null(wasserstein$cost)) stop("Wasserstein list must be provided")
  if (method == "Constrained Wasserstein" & is.null(wasserstein$cost)) stop("Wasserstein list must be provided")
  
  opt.class <- switch(method,
                      "Wasserstein" = cotDualL2_2)
  
  optimizer <- opt.class$new(lambda = wasserstein$lambda, cost = wasserstein$cost, 
                             b = wasserstein$b,
                             marginal.costs = marginal.wasserstein$marginal.costs,
                             marginal.delta = marginal.wasserstein$marginal.delta,
                             balance.functions = balance$balance.functions,
                             balance.delta = balance$balance.delta,
                             target.mean = balance$target.mean,
                             target.sd = balance$target.sd)
  
  control <- control.options(control, method = "lbfgs")
  if (is.null(init)) init <- optimizer$init()
  bounds <- optimizer$bounds()
  
  fit <- lbfgsb3c::lbfgsb3c(par = init,
                            fn = optimizer$obj,
                            gr = optimizer$grad,
                            lower = bounds[,1],
                            upper = bounds[,2],
                            control = control
  )
  if (is.null(fit$convergence) | isFALSE(fit$convergence == 0)) warning(fit$message)
  beta <- c(fit$par)
  
  
  weight <- optimizer$get_weight(beta)
  # unconst_weight <- as.numeric(optimizer$h_star_inv(beta))
  # # unconst_weight <- unconst_weight * as.numeric(unconst_weight > 0)
  # weight <- simplex_proj(unconst_weight)
  
  if (grepl("Wasserstein", method)) {
    n <- nrow(wasserstein$cost)
    m <- ncol(wasserstein$cost)
    
    # unconst_weight <- Matrix::Matrix(data = unconst_weight,
    #                                  nrow = n, ncol = m)
    weight <- if(is.matrix(weight)) {
                Matrix::Matrix(data = weight)
    } else {
      Matrix::Matrix(data = weight,
                     nrow = n, ncol = m)
    }
                             
  }
  
  return(list(weight = weight, beta = beta,
              # unconstrained = list(weight = unconst_weight),
              fit = fit,
              optimizer = optimizer))
}

#### deriv of divergence ####
# function(f,g,a,b)