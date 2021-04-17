sbw_dual <- function(x, target, constraint) {
  
  n  <- nrow(x)
  x_m <- colMeans(x)
  x_v <- matrixStats::colVars(x)
  
  if (!is.null(nrow(target)) ) {
    if (nrow(target) > 1 )  {
      S <- sqrt(0.5 * x_v + 0.5 * matrixStats::colVars(target))
      target <- colMeans(target)
    }
  } else {
    S <- sqrt(x_v)
  }
  
  if ( constraint == 0) {
    # equivalent to (x'x)^(-1) (x' 1_n - target)
    QR <- qr(x)
    R <- qr.R(QR)
    
    beta <- 2 * (backsolve(R, forwardsolve(l = R, target,
                                                         transpose = TRUE))
                  - qr.coef(QR, rep(1,ncol(x))) )
    
  } else {
    # use lasso in oem that can take xtx and xty
    
    XtY <- (target - x_m) #* 1/S)
    # XtX <- crossprod(scale(x, center = FALSE, scale = S))
    XtX <- crossprod(x)/2
    
    fit <- oem::oem.xtx(xtx = XtX, xty = XtY, family = "gaussian",
                        penalty = "lasso", lambda = constraint
                        , penalty.factor = S
    )
    beta <- fit$beta$lasso
  }
  
 
  
  unconst_wt  <- sbw_dual_to_wt(x, lambda = beta, n = n)
  
  return(list(weight = simplex_proj(unconst_wt), 
              unconstrained_weight = unconst_wt,
              lambda = beta))
}

sbw_dual_to_wt <- function(x, lambda, n) {
  eta <- (x %*% lambda )
  return(-(eta)/2.0 + 1.0 / n)
}


cot_const_dual <- function(x, target, sample_weights = NULL, constraint, metric, power, 
                           add.margins = FALSE, marg.constraints = NULL,
                           formula = NULL, balance.constraints = NULL) {
  n <- nrow(x)
  m <- nrow(target)
  
  z <- c( rep(0, n), 
          rep(1, m))
  Q_mat <- rbind(rep(1, n*m),
                 vec_to_col_constraints(n,m),
                 cost_fun(x = rbind(x, target), z = z, 
                          metric = metric, ground_power = power, estimand = "ATT")^power
  )
  
  sw    <- get_sample_weight(sample_weights, z)
  d <- c(1, sw$b, constraint)
  
  if (add.margins) {
    combined.mat <- rbind(x, target)
    Q_mat <- do.call("rbind", c(list(c(Q_mat)),
                                lapply( 1:ncol(x), function(d)
                                  c(cost_fun(x = combined.mat[,d, drop = FALSE], z = z, 
                                           metric = metric, ground_power = power, estimand = "ATT")^power
                                  ))))
    if (length(marg.constraints) == 1) marg.constraints <- rep(marg.constraints, ncol(x))
    stopifnot(length(marg.constraints) == ncol(x))
    d <- c(d, marg.constraints)
  }
  
  
  if (!is.null(form)  & !is.null(balance.constraints)) {
    form    <- form_all_squares(formula, colnames(target))
    targ_mm <- model.matrix(form, data = data.frame(target))
    x_mm    <-  model.matrix(form, data = data.frame(x))
    S       <- sqrt(0.5 * matrixStats::weightedVar(targ_mm, w = sw$b) + 
                      0.5 * matrixStats::weightedVar(x_mm, w = sw$a))
    t_m     <- c(rep(0, nrow(Q_mat)), matrixStats::colWeightedMeans(targ_mm, w = sw$b))
    # x_m     <- matrixStats::colWeightedMeans(x_mm, w = sw$a)
    Q_mat   <- rbind(Q_mat,
                     Matrix::crossprod(x_mm, vec_to_col_constraints(n , m)))
    d  <- c(d , balance.constraints * S)
    
    XtY <- rowMeans(Q_mat) - t_m
  } else {
    XtY <- rowMeans(Q_mat)
  }
  XtX <- crossprod(Q_mat)/2
  
  fit <- oem::oem.xtx(xtx = XtX, xty = XtY, family = "gaussian",
                      penalty = "lasso", lambda = 1, penalty.factor = d)
  
  beta <- fit$beta$lasso
  
  unconst_wt  <- cot_const_dual_to_wt(x = x, lambda = beta, n = ncol(Q_mat))
  
  return(list(weight = simplex_proj(unconst_wt), 
              unconstrained_weight = unconst_wt,
              lambda = beta))
}

cot_const_dual_to_wt <- function(x, lambda, n) {
  eta <- (x %*% lambda )
  return((eta)/2.0 + 1.0 / n)
}

cot_dual <- function(x, target, sample_weights = NULL, constraint, metric, power, 
                     add.margins = FALSE, marg.constraints = NULL,
                     formula = NULL, balance.constraints = NULL) {
  
}

optProblem <- R6::R6Class("optProblem", 
                          public = list(
                            get_xtx = function() NULL,
                            get_xty = function() NULL
                            # bounds = "function",
                            # obj  = "function",
                            # grad = "function",
                            # get_weight = "function",
                            # intialize = function(delta) {
                            #   private$delta <- delta
                            #   invisible(self)
                            # }
                          ),
                          private = list(
                            # "n" = integer
                          )
                          
)

sbwDualL2 <- R6::R6Class("sbwDualL2",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                              return(sum(private$delta * private$scale * abs(vars)) +
                               sum(self$h_star_inv(vars)^2) + 
                               sum(private$target.functions * vars))
                           },
                           grad = function(vars) {
                             lambda <- private$delta * private$scale
                             # vars <- vars * as.numeric(abs(vars) > lambda)
                             return(sign(vars) * lambda +
                                c(private$target.mean) -
                               c(crossprod(private$balance.functions, self$h_star_inv(vars))))
                           },
                           h_star_inv = function(vars) {
                             eta <- -0.5 * private$balance.functions %*% vars + 1.0/private$n
                             return(eta * as.numeric(eta > 0))
                           },
                           get_weight = function(vars) {
                             return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                          },
                          bounds = function() {
                            cbind(rep(-Inf, ncol(private$balance.functions)),
                                  rep(Inf ,  ncol(private$balance.functions)))
                          },
                          init = function() {
                            rep(0, private$nvars)
                          },
                          get_xtx = function() {
                            return(0.5 * crossprod(private$balance.functions))
                          },
                          get_xty = function() {
                            return(
                              c(colMeans(private$balance.functions) - private$target.mean)
                            )
                          },
                          penalty.factor = function() {
                            return(
                              private$delta * private$scale
                            )
                          },
                         initialize = function(balance.delta,
                                               balance.functions,
                                               target.mean = NULL,
                                               target.sd = NULL,
                                               ...
                         ) {
                           private$delta <- balance.delta
                           private$balance.functions     <- as.matrix(balance.functions)
                           private$n     <- nrow(private$balance.functions)
                           private$target.mean <- colMeans(as.matrix(target.mean))
                           
                           private$scale <- if (missing(target.sd) | is.null(target.sd)) {
                             sqrt(0.5 * matrixStats::colVars( private$balance.functions ) +
                             0.5 * matrixStats::colVars( target.mean ))
                           } else {
                            sqrt(0.5 * matrixStats::colVars( private$balance.functions ) +
                                 0.5 * target.sd^2) 
                           }
                           
                           private$nvars <- ncol(private$balance.functions)
                         }
                        ),
                      private = list(
                        n = "integer",
                        nvars = "integer",
                        scale = "numeric",
                        balance.functions = "matrix",
                        target.mean = "numeric",
                        delta = "numeric"
                      )
)

cotDualL2 <- R6::R6Class("cotDualL2",
                         inherit = optProblem,
                         public = list(
                             obj =  function(vars) {
                               return(-(private$p.fun(vars) - 0.5 * sum(self$h_star_inv(vars)^2)))
                             },
                             grad = function(vars) {
                               
                               return(as.numeric(-(private$p.grad.fun(vars) - Matrix::crossprod(private$Q, 
                                                                    self$h_star_inv(vars)))))
                             },
                             h_star_inv = function(vars) {
                               eta <- (private$Q %*% vars - private$cost) / private$delta
                               return(eta * as.numeric(eta > 0))
                             },
                             get_weight = function(vars) {
                               return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
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
                                   cbind(rep(-Inf, l_bd ),
                                         rep(Inf ,  l_bd))
                                   )
                               } else if (private$margins & !private$balfun) {
                                 rbind(
                                   cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m)),
                                   cbind(rep(0, l_md),
                                         rep(Inf, l_md)))
                                 
                               } else if (!private$margins & private$balfun) {
                                 cbind(rep(-Inf, private$m + l_bd),
                                       rep(Inf ,  private$m + l_bd))
                               } else {
                                 cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m))
                               }
                               return(bounds)
                             },
                             get_nvars = function(){
                               return(private$nvars)
                             },
                             init = function() {
                               rep(0, private$nvars)
                             },
                             get_xtx = function() {
                               neg.mass.const <- -private$Q[,private$dual.idx]
                               cmb.Q <- cbind(neg.mass.const,
                                              private$Q)
                               return(Matrix::crossprod(cmb.Q) / private$delta)
                             },
                             get_xty = function() {
                               neg.mass.const <- -private$Q[,private$dual.idx]
                               cmb.Q <- cbind(neg.mass.const,
                                              private$Q)
                               QtC  <- as.numeric(Matrix::crossprod(cmb.Q, 
                                                                   private$cost)) / private$delta
                               l_t  <- length(private$target.mean)
                               A    <- c(rep(0, length(QtC) - l_t), private$target.mean)
                               return(QtC + A)
                             },
                             # get_X = function() {
                             #   neg.mass.const <- -private$Q[,private$dual.idx]
                             #   cmb.Q <- cbind(neg.mass.const,
                             #                  private$Q)
                             #   A <- cbind(
                             #     matrix(0, nrow = length(private$target.mean),
                             #          ncol = ncol(cmb.Q) - length(private$target.mean)),
                             #     Matrix::Diagonal(length(private$target.mean),
                             #                      x = private$target.mean))
                             #   return(rbind(cmb.Q / sqrt(private$delta),
                             #                 A) )
                             # },
                             # get_Y = function() {
                             #   return(c(private$cost / sqrt(private$delta), 
                             #            rep(1, length(private$target.mean))))
                             # },
                             penalty.factor = function() {
                               return(private$pf)
                             },
                           initialize = function(delta, cost, 
                                                 b,
                                                 marginal.costs = NULL,
                                                 marginal.delta = NULL,
                                                 balance.functions = NULL,
                                                 balance.delta = NULL,
                                                 target.mean = NULL,
                                                 target.sd = NULL,
                                                 ...
                                                 ) {
                             private$delta <- delta
                             private$cost <- c(cost)
                             private$n    <- nrow(cost)
                             private$m    <- ncol(cost)
                             private$dual.idx <- 1:private$m
                             cur.idx <- private$m + 1
                             
                             private$b <- b
                             
                             if (!is.null(marginal.costs) & 
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
                               private$bal.idx <- cur.idx:(cur.idx + ncol(balance.functions)-1)
                               private$balance.functions <- Matrix::crossprod(vec_to_row_constraints(private$n, private$m),
                                                                              balance.functions)
                               scale <- if (missing(target.sd) | is.null(target.sd)) {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        0.5 * matrixStats::colVars( target.mean ))
                               } else {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        target.sd^2) 
                               }
                               private$balance.delta <- balance.delta * scale
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
                                                  "1" = Matrix::t(vec_to_col_constraints(private$n,
                                                                                             private$m)),
                                                  "2" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                             private$m)),
                                                            private$marginal.costs),
                                                  "3" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                    private$balance.functions),
                                                  "4" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                    private$marginal.costs,
                                                                    private$balance.functions))
                             private$nvars = switch(fun.num,
                                                    "1" = private$m,
                                                    "2" = private$m + ncol(private$marginal.costs),
                                                    "3" = private$m + ncol(private$balance.functions),
                                                    "4" = private$m +
                                                      ncol(private$marginal.costs) + 
                                                      ncol(private$balance.functions))
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
                           delta = "numeric",
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
                                     -sign(beta_b) * private$balance.delta -
                                       private$target.mean)
                             )
                           }, 
                           p.grad.fun.cmb = function(vars) {
                             beta_b <- vars[private$bal.idx]
                             return( c(private$b,
                                       -private$marginal.delta,
                                       -sign(beta_b) * private$balance.delta -
                                       private$target.mean)
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
                                      sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           }, 
                           p.fun.cmb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             beta_b <- vars[private$bal.idx]
                             return( sum(private$b * g) -
                                       sum(private$marginal.delta * beta_m) -
                                       sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           },
                           h_star_inv.c = function(vars) {
                             eta <- rep(private$vars, each = private$n) - private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cm = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             eta <- rep(g, each = private$n) +  
                                private$marginal.costs %*% beta_m - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             eta <- rep(g, each = private$n) +  
                               private$marginal.costs %*% beta_m - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cmb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             beta_b <- vars[private$bal.idx]
                             eta <- rep(g, each = private$n) +  
                               private$marginal.costs %*% beta_m +
                               Matrix::crossprod(private$balance.functions, beta_b) - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           }
                         )
)

cotConstDualL2 <- R6::R6Class("cotConstDualL2",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                             return(-(private$p.fun(vars) - sum(
                               self$h_star_inv(vars)^2
                               # (1 / (private$n * private$m) - 0.5 * private$Q %*% vars)^2
                               )))
                           },
                           grad = function(vars) {
                             
                             return(as.numeric(-(private$p.grad.fun(vars) + Matrix::crossprod(private$Q, 
                                                                           self$h_star_inv(vars)
                                                                           # 1 / (private$n * private$m) - 0.5 * private$Q %*% vars
                                                                           ))))
                           },
                           h_star_inv = "function",
                           # h_star_inv = function(vars) {
                           #   eta <-  1 / (private$n * private$m) - 0.5 * private$Q %*% vars
                           #   return(eta)
                           #   return(eta * as.numeric(eta > 0))
                           # },
                           get_weight = function(vars) {
                             return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                           },
                           bounds = function(vars) {
                             l_md <- length(private$marginal.delta)
                             l_bd <- length(private$balance.delta)
                             bounds <- if (private$margins & private$balfun) {
                               rbind(
                                 cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m)),
                                 cbind(0, Inf),
                                 cbind(rep(0, l_md),
                                       rep(Inf, l_md)),
                                 cbind(rep(-Inf, l_bd ),
                                       rep(Inf ,  l_bd))
                               )
                             } else if (private$margins & !private$balfun) {
                               rbind(
                                 cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m)),
                                 cbind(0, Inf),
                                 cbind(rep(0, l_md),
                                       rep(Inf, l_md)))
                               
                             } else if (!private$margins & private$balfun) {
                               rbind(cbind(rep(-Inf, private$m ),
                                           rep(Inf ,  private$m)),
                                     cbind(0, Inf),
                                     cbind(rep(-Inf, l_bd),
                                     rep(Inf ,   l_bd)))
                             } else {
                               rbind(cbind(rep(-Inf, private$m ),
                                           rep(Inf ,  private$m)),
                               cbind(0, Inf))
                             }
                             return(bounds)
                           },
                           get_nvars = function(){
                             return(private$nvars)
                           },
                           init = function() {
                             rep(0, private$nvars)
                           },
                           get_xtx = function() {
                             neg.mass.const <- -private$Q[,private$dual.idx]
                             cmb.Q <- cbind(neg.mass.const,
                                            private$Q)
                             return(Matrix::crossprod(cmb.Q) * 0.5)
                           },
                           get_xty = function() {
                             neg.mass.const <- -private$Q[,private$dual.idx]
                             cmb.Q <- cbind(neg.mass.const,
                                            private$Q)
                             QtO  <- as.numeric(Matrix::colMeans(cmb.Q))
                             l_t  <- length(private$target.mean)
                             A    <- c(rep(0, length(QtO) - l_t), 
                                       private$target.mean)
                             return( A - QtO )
                           },
                           penalty.factor = function() {
                             return(private$pf)
                           },
                           initialize = function(delta, cost, 
                                                 b,
                                                 marginal.costs = NULL,
                                                 marginal.delta = NULL,
                                                 balance.functions = NULL,
                                                 balance.delta = NULL,
                                                 target.mean = NULL,
                                                 target.sd = NULL,
                                                 ...
                           ) {
                             private$delta <- delta
                             private$cost <- c(cost)
                             private$n    <- nrow(cost)
                             private$m    <- ncol(cost)
                             private$dual.idx <- 1:private$m
                             private$cwass.idx <- private$m + 1
                             cur.idx <- private$m + 2
                             
                             private$b <- b
                             
                             if (!is.null(marginal.costs) & 
                                 !is.null(marginal.delta) ) {
                               private$margins <- TRUE
                               private$marg.idx <- cur.idx:(cur.idx + length(marginal.costs) - 1)
                               cur.idx <- cur.idx + length(marginal.costs) 
                               private$marginal.costs <- sapply(marginal.costs, c)
                               private$marginal.delta <- marginal.delta
                             } else {
                               private$margins <- FALSE
                             } 
                             
                             if (!is.null(balance.functions) & 
                                 !is.null(balance.delta)) {
                               private$balfun <- TRUE
                               private$bal.idx <- cur.idx:(cur.idx + ncol(balance.functions) - 1)
                               private$balance.functions <- Matrix::crossprod(vec_to_row_constraints(private$n, private$m),
                                                                              balance.functions)
                               scale <- if (missing(target.sd) | is.null(target.sd)) {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        0.5 * matrixStats::colVars( target.mean ))
                               } else {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        target.sd^2) 
                               }
                               private$balance.delta <- balance.delta * scale
                               private$target.mean <- target.mean
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
                                                  "1" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                              private$cost),
                                                  "2" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                    private$cost,
                                                                    private$marginal.costs),
                                                  "3" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                    private$cost,
                                                                    private$balance.functions),
                                                  "4" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                     private$m)),
                                                                    private$cost,
                                                                    private$marginal.costs,
                                                                    private$balance.functions))
                             private$nvars = switch(fun.num,
                                                    "1" = private$m + 1,
                                                    "2" = private$m + 1 + ncol(private$marginal.costs),
                                                    "3" = private$m + 1 + ncol(private$balance.functions),
                                                    "4" = private$m + 1 +
                                                      ncol(private$marginal.costs) + 
                                                      ncol(private$balance.functions))
                             
                             private$pf <- switch(fun.num,
                                                  "1" = c(-private$b, private$b, 
                                                          private$delta),
                                                  "2" = c(-private$b, private$b, private$delta, 
                                                          private$marginal.delta),
                                                  "3" = c(-private$b, private$b, private$delta, 
                                                          private$balance.delta),
                                                  "4" = c(-private$b, private$b, private$delta, 
                                                          private$marginal.delta,
                                                          private$balance.delta
                                                          )
                                                  )
                             
                             self$h_star_inv <- switch(fun.num,
                                                       "1" = private$h_star_inv.c,
                                                       "2" = private$h_star_inv.cm,
                                                       "3" = private$h_star_inv.cb,
                                                       "4" = private$h_star_inv.cmb
                                                       )
                             
                           }
                         ),
                         private = list(
                           b  = "numeric",
                           n = "integer",
                           m = "integer",
                           nvars = "integer",
                           cost = "numeric",
                           delta = "numeric",
                           margins = "logical",
                           balfun  = "logical",
                           marginal.costs = "list",
                           marginal.delta = "numeric",
                           balance.delta = "numeric",
                           balance.functions = "matrix",
                           Q = "matrix",
                           pf = "numeric",
                           dual.idx = "integer",
                           cwass.idx = "integer",
                           marg.idx = "integer",
                           bal.idx = "integer",
                           target.mean = "numeric",
                           p.grad.fun = "function",
                           p.fun = "function",
                           p.grad.fun.c = function(vars) {
                             return(c(-private$b, -private$delta))
                           },
                           p.grad.fun.cm = function(vars) {
                             return( c(-private$b,
                                       -private$delta,
                                       -private$marginal.delta)
                             )
                           },
                           p.grad.fun.cb = function(vars) {
                             beta_b <- vars[private$bal.idx]
                             return( c(-private$b , 
                                       -private$delta,
                                       -sign(beta_b) * private$balance.delta -
                                       private$target.mean)
                             )
                           }, 
                           p.grad.fun.cmb = function(vars) {
                             beta_b <- vars[private$bal.idx]
                             return( c( -private$b,
                                     -private$delta,
                                     -private$marginal.delta ,
                                     -sign(beta_b) * private$balance.delta -
                                       private$target.mean)
                             )
                           },
                           p.fun.c = function(vars) {
                             g <- vars[private$dual.idx]
                             e <- vars[private$cwass.idx]
                             return(sum(-private$b * g) - private$delta * e)
                           },
                           p.fun.cm = function(vars) {
                             g <- vars[private$dual.idx]
                             e <- vars[private$cwass.idx]
                             beta_m <- vars[private$marg.idx]
                             return( sum(-private$b * g) - private$delta * e - 
                                       sum(private$marginal.delta * beta_m)
                             )
                           },
                           p.fun.cb = function(vars) {
                             g <- vars[private$dual.idx]
                             e <- vars[private$cwass.idx]
                             beta_b <- vars[private$bal.idx]
                             return( sum(-private$b * g) - private$delta * e -
                                       sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           }, 
                           p.fun.cmb = function(vars) {
                             g <- vars[private$dual.idx]
                             e <- vars[private$cwass.idx]
                             beta_m <- vars[private$marg.idx]
                             beta_b <- vars[private$bal.idx]
                             return( sum(-private$b * g) - private$delta * e -
                                       sum(private$marginal.delta * beta_m) -
                                       sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           },
                           h_star_inv.c = function(vars) {
                             eta <- -private$Q %*% vars * 0.5 + 
                               1 / (private$n * private$m)
                             return(eta)
                             return(eta * as.numeric(eta > 0))
                           },
                           h_star_inv.cm = function(vars) {
                             eta <- -private$Q %*% vars * 0.5 + 
                               1 / (private$n * private$m)
                             return(eta)
                             return(eta * as.numeric(eta > 0))
                           },
                           h_star_inv.cb = function(vars) {
                             eta <-  -private$Q %*% vars * 0.5 + 
                               1 / (private$n * private$m)
                             return(eta)
                             return(eta * as.numeric(eta > 0))
                           },
                           h_star_inv.cmb = function(vars) {
                             eta <-  -private$Q %*% vars * 0.5 + 
                               1 / (private$n * private$m)
                             return(eta)
                             return(eta * as.numeric(eta > 0))
                           }
                           
                         )
)

balance.options <- function(balance, x, target) {
  if (is.null(balance)) {
    return(list(NULL))
  #   balance <- list()
  #   balance$balance.functions <- model.matrix(~., data.frame(x))
  #   target_m <- model.matrix(~., data.frame(target))
  #   balance$target.sd <- matrixStats::colSds(target_m)
  #   balance$target.mean <- matrix(colMeans(target_m),nrow = 1)
  #   balance$balance.delta <- 0.2
  }
  
  if (is.null(balance$formula)) {
    balance$formula <- as.formula(~. + 0)
  } else {
    balance$formula <- as.formula(balance$formula)
    environment(balance$formula) <- environment()
  }
  
  if (is.null(balance$balance.functions) ) {
    balance$balance.functions <- model.matrix(balance$formula, data.frame(x))
  }
  
  if (is.null(balance$target.sd) | is.null(balance$target.mean)) target_m <- model.matrix(balance$formula, data.frame(target))
  if (is.null(balance$target.mean))  balance$target.mean <- matrix(colMeans(target_m),nrow = 1)
  if (is.null(balance$target.sd)) balance$target.sd <- matrixStats::colSds(target_m)
  # if (is.null(balance$balance.delta)) balance$balance.delta <- 0.2
  
  balance$balance.functions <- as.matrix(balance$balance.functions)
  balance$balance.delta <- as.double(balance$balance.constraints)
  balance$target.mean <- as.matrix(balance$target.mean)
  balance$target.sd <- as.double(balance$target.sd)
  if(length(balance$balance.delta) == 0) balance <- list(NULL)
  return(balance)
  
}

wasserstein.options <- function(wasserstein, x, target, sw) {
    if (is.null(wasserstein) | missing(wasserstein)) {
      # return(list(NULL))
      wasserstein        <- list(metric = "mahalanobis")
      wasserstein$power  <- 2
      z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
      wasserstein$cost   <- cost_fun(rbind(x, target), z = z, 
                                     estimand = "ATT", power = 2, 
                                     metric = "mahalanobis")
      wasserstein$lambda  <- 1.0
      wasserstein$b <- sw$b
    }
    if (is.null(wasserstein$metric)) wasserstein$metric <- match.arg(wasserstein$metric, c("mahalanobis", "sdLp", "Lp"))
    if (is.null(wasserstein$power)) wasserstein$power <- 2
    if (is.null(wasserstein$cost)) {
      z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
      wasserstein$cost <- cost_fun(rbind(x, target),
                                   z = z,
                                  estimand = "ATT", power = 2, 
                                  metric = "mahalanobis")
    }
    if (is.null(wasserstein$lambda)) wasserstein$lambda <- 1.0
    if (is.null(wasserstein$b)) wasserstein$b <- sw$b
    
    wasserstein$metric <- match.arg(wasserstein$metric, c("mahalanobis", "sdLp", "Lp"))
    wasserstein$power  <- as.double(wasserstein$power)
    wasserstein$cost   <- as.matrix(wasserstein$cost)^wasserstein$power
    wasserstein$lambda  <- as.double(wasserstein$lambda)
    wasserstein$b      <- as.double(wasserstein$b)
    
    stopifnot(wasserstein$power > 0)
    stopifnot(wasserstein$lambdaa > 0)
    
    return(wasserstein)
}

marginal.wass.options <- function(marg.wass, wass, x, target) {
  if (is.null(marg.wass) | missing(marg.wass)) {
    return(list(NULL))
    marg.wass        <- list()
  }
  if (is.null(marg.wass$marginal.constraints)) return(list(NULL))
  if (is.null(marg.wass$marginal.costs)) {
    z <- c(rep(0, nrow(x)), rep(1, nrow(target)))
           
    marg.wass$marginal.costs <- lapply(1:ncol(x), function(i) 
                                      cost_fun(rbind(x[,i,drop = FALSE], 
                                               target[,i,drop = FALSE]),
                                               z = z,
                                          estimand = "ATT", power = wass$power, 
                                          metric = wass$metric))
  }
  # if (is.null(marg.wass$marginal.delta)) {
  #   marg.wass$marginal.delta <- rep(1.0^wass$power, length(marg.wass$marginal.costs))
  # }
  
  marg.wass$marginal.costs  <- lapply(marg.wass$marginal.costs, function(cc) cc^wass$power)
  marg.wass$marginal.delta  <- as.double(marg.wass$marginal.constraints^wass$power)
  if (length(marg.wass$marginal.delta) != length(marg.wass$marginal.costs)) marg.wass$marginal.delta <- rep(marg.wass$marginal.delta, 
                                                                                                           length(marg.wass$marginal.costs))
  
  stopifnot(all(marg.wass$marginal.delta > 0))
  
  return(marg.wass)
}

control.options <- function(control, method) {
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
    
  } else if (method == "lbfgs") {
    
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

dual_opt <- function(x, target, 
                     init = NULL,
                     sample_weights = NULL, 
                     method = c("SBW", "Wasserstein",
                                "Constrained Wasserstein"),
                     wasserstein = list(metric = c("mahalanobis"),
                                        power = 2,
                                        cost = NULL,
                                        delta = 1.0),
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
         "SBW" = sbwDualL2,
         "Constrained Wasserstein" = cotConstDualL2,
         "Wasserstein" = cotDualL2)
  
  optimizer <- opt.class$new(delta = wasserstein$delta, cost = wasserstein$cost, 
          b = wasserstein$b,
          marginal.costs = marginal.wasserstein$marginal.costs,
          marginal.delta = marginal.wasserstein$marginal.delta,
          balance.functions = balance$balance.functions,
          balance.delta = balance$balance.delta,
          target.mean = balance$target.mean,
          target.sd = balance$target.sd)
  
  if (!is.null(balance$balance.functions) & !is.null(balance$balance.delta)) {
    # if (method != "Wasserstein") {
      control <- control.options(control, method = "oem")
      XtX <- as.matrix(optimizer$get_xtx())
      XtY <- as.matrix(optimizer$get_xty())
      
      fit <- oem::oem.xtx(xtx = XtX, 
                          xty = XtY,
                   family = "gaussian",
                   penalty = "lasso",
                   lambda = 1, nlambda = 1,
                   lambda.min.ratio = NULL,
                   alpha = 1, gamma = 3, tau = 0.5, #not used,
                   groups = control$groups,
                   scale.factor = control$scale.factor,
                   penalty.factor = optimizer$penalty.factor(),
                   group.weights = control$group.weights,
                   maxit = control$maxit,
                   tol = control$tol,
                   irls.maxit = control$irls.maxit,
                   irls.tol = control$irls.tol
                   )
      beta <- c(fit$beta[[1]])
      if (grepl("Wasserstein", method)) {
        neg.dual.idx <- optimizer$.__enclos_env__$private$dual.idx
        pos.dual.idx <- neg.dual.idx + length(neg.dual.idx)
        beta_neg <- beta[neg.dual.idx]
        beta_pos <- beta[pos.dual.idx]
        beta <- c(beta_pos - beta_neg, beta[-c(neg.dual.idx, pos.dual.idx)])
      }
    # } 
    # else {
    #   X <- optimizer$get_X()
    #   Y <- optimizer$get_Y()
    #   bounds <- optimizer$bounds()
    #   neg.dual.idx <- optimizer$.__enclos_env__$private$dual.idx
    #   pos.dual.idx <- neg.dual.idx + length(neg.dual.idx)
    #   bounds[neg.dual.idx,1] <- 0
    #   bounds <- rbind(bounds[neg.dual.idx,,drop = FALSE],
    #                   bounds)
    #   
    #   fit <- glmnet::glmnet(x = X, y = Y, family = "gaussian",
    #                  lower.limits = bounds[,1],
    #                  upper.limits = bounds[,2],
    #                  maxit = control$maxit, relax = FALSE,
    #                  penalty.factor = optimizer$penalty.factor(),
    #                  intercept = FALSE, lambda = 1)
    #   beta <- c(fit$beta)
    #   beta_neg <- beta[neg.dual.idx]
    #   beta_pos <- beta[pos.dual.idx]
    #   beta <- c(beta_pos - beta_neg, beta[-c(neg.dual.idx, pos.dual.idx)])
    #   
    # }
  } else {
    control <- control.options(control, method = "lbfgs")
    if (is.null(init)) init <- optimizer$init()
    bounds <- optimizer$bounds()
    
    fit <- lbfgsb3c::lbfgsb3(par = init,
                             fn = optimizer$obj,
                             gr = optimizer$grad,
                             lower = bounds[,1],
                             upper = bounds[,2],
                             control = control
    )
    if (is.null(fit$convergence) | isFALSE(fit$convergence == 0)) warning(fit$message)
    beta <- c(fit$par)
  }
  
  
  
  unconst_weight <- as.numeric(optimizer$h_star_inv(beta))
  unconst_weight <- unconst_weight * as.numeric(unconst_weight > 0)
  weight <- simplex_proj(unconst_weight)
  
  if (grepl("Wasserstein", method)) {
    n <- nrow(wasserstein$cost)
    m <- ncol(wasserstein$cost)
    
    unconst_weight <- Matrix::Matrix(data = unconst_weight,
                                     nrow = n, ncol = m)
    weight <- Matrix::Matrix(data = weight,
                                     nrow = n, ncol = m)
  }
  
  return(list(weight = weight, beta = beta,
              unconstrained = list(weight = unconst_weight),
              fit = fit))
}

convergence <- function(old, new) {
  
}

# simplex_opt_lbfgs <- function(optimizer, control) {
#   control <- control.options(control, method = "lbfgs")
#   if (is.null(init)) init <- optimizer$init()
#   bounds <- optimizer$bounds()
#   lfit <- simp_pred <- cur_param <- old_param <- list(NULL)
#   
#   old_param[[1]] <- init
#   
#   X <- optimizer$get_Q()
#   
#   for (i in 1:control$maxit) {
#     lfit[[1]] <- lbfgsb3c::lbfgsb3(par = old_param[[1]],
#                     fn = optimizer$obj,
#                     gr = optimizer$grad,
#                     lower = bounds[,1],
#                     upper = bounds[,2],
#                     control = control
#   )$par
#     simp_pred[[1]] <- optimizer$get_weight(lfit[[1]])
#     cur_param[[1]] <- lm(optimizer$get_y_pgd(simp_pred[[1]]) ~ X + 0)
#     if (convergence(old_param[[1]], cur_param[[1]])) {
#       break
#     } else {
#       old_param[[1]] <- cur_param[[1]]
#     }
#   }
#   
#   return(cur_param[[1]])
# }
# 
# convergence <- function()

fullDualL2 <- R6::R6Class("fullDualL2",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                             f <- vars[1:private$n0]
                             g <- vars[-c(1:private$n0)]
                             
                             return(sum(private$alpha * f) + private$beta * g)
                           },
                           grad = function(vars) {
                             
                             return(as.numeric(-(private$p.grad.fun(vars) - Matrix::crossprod(private$Q, 
                                                                                              self$h_star_inv(vars)))))
                           },
                           h_star_inv = function(vars) {
                             eta <- (private$Q %*% vars - private$cost) / private$delta
                             return(eta * as.numeric(eta > 0))
                           },
                           get_weight = function(vars) {
                             return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
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
                                 cbind(rep(-Inf, l_bd ),
                                       rep(Inf ,  l_bd))
                               )
                             } else if (private$margins & !private$balfun) {
                               rbind(
                                 cbind(rep(-Inf, private$m ),
                                       rep(Inf ,  private$m)),
                                 cbind(rep(0, l_md),
                                       rep(Inf, l_md)))
                               
                             } else if (!private$margins & private$balfun) {
                               cbind(rep(-Inf, private$m + l_bd),
                                     rep(Inf ,  private$m + l_bd))
                             } else {
                               cbind(rep(-Inf, private$m ),
                                     rep(Inf ,  private$m))
                             }
                             return(bounds)
                           },
                           get_nvars = function(){
                             return(private$nvars)
                           },
                           init = function() {
                             rep(0, private$nvars)
                           },
                           get_xtx = function() {
                             neg.mass.const <- -private$Q[,private$dual.idx]
                             cmb.Q <- cbind(neg.mass.const,
                                            private$Q)
                             return(Matrix::crossprod(cmb.Q) / private$delta)
                           },
                           get_xty = function() {
                             neg.mass.const <- -private$Q[,private$dual.idx]
                             cmb.Q <- cbind(neg.mass.const,
                                            private$Q)
                             QtC  <- as.numeric(Matrix::crossprod(cmb.Q, 
                                                                  private$cost)) / private$delta
                             l_t  <- length(private$target.mean)
                             A    <- c(rep(0, length(QtC) - l_t), private$target.mean)
                             return(QtC + A)
                           },
                           # get_X = function() {
                           #   neg.mass.const <- -private$Q[,private$dual.idx]
                           #   cmb.Q <- cbind(neg.mass.const,
                           #                  private$Q)
                           #   A <- cbind(
                           #     matrix(0, nrow = length(private$target.mean),
                           #          ncol = ncol(cmb.Q) - length(private$target.mean)),
                           #     Matrix::Diagonal(length(private$target.mean),
                           #                      x = private$target.mean))
                           #   return(rbind(cmb.Q / sqrt(private$delta),
                           #                 A) )
                           # },
                           # get_Y = function() {
                           #   return(c(private$cost / sqrt(private$delta), 
                           #            rep(1, length(private$target.mean))))
                           # },
                           penalty.factor = function() {
                             return(private$pf)
                           },
                           initialize = function(delta, cost, 
                                                 b,
                                                 marginal.costs = NULL,
                                                 marginal.delta = NULL,
                                                 balance.functions = NULL,
                                                 balance.delta = NULL,
                                                 target.mean = NULL,
                                                 target.sd = NULL,
                                                 ...
                           ) {
                             private$delta <- delta
                             private$cost <- c(cost)
                             private$n    <- nrow(cost)
                             private$m    <- ncol(cost)
                             private$dual.idx <- 1:private$m
                             cur.idx <- private$m + 1
                             
                             private$b <- b
                             
                             if (!is.null(marginal.costs) & 
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
                               private$bal.idx <- cur.idx:(cur.idx + ncol(balance.functions)-1)
                               private$balance.functions <- Matrix::crossprod(vec_to_row_constraints(private$n, private$m),
                                                                              balance.functions)
                               scale <- if (missing(target.sd) | is.null(target.sd)) {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        0.5 * matrixStats::colVars( target.mean ))
                               } else {
                                 sqrt(0.5 * matrixStats::colVars( balance.functions ) +
                                        target.sd^2) 
                               }
                               private$balance.delta <- balance.delta * scale
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
                                                  "1" = Matrix::t(vec_to_col_constraints(private$n,
                                                                                         private$m)),
                                                  "2" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                               private$m)),
                                                              private$marginal.costs),
                                                  "3" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                               private$m)),
                                                              private$balance.functions),
                                                  "4" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                               private$m)),
                                                              private$marginal.costs,
                                                              private$balance.functions))
                             private$nvars = switch(fun.num,
                                                    "1" = private$m,
                                                    "2" = private$m + ncol(private$marginal.costs),
                                                    "3" = private$m + ncol(private$balance.functions),
                                                    "4" = private$m +
                                                      ncol(private$marginal.costs) + 
                                                      ncol(private$balance.functions))
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
                           delta = "numeric",
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
                                       -sign(beta_b) * private$balance.delta -
                                         private$target.mean)
                             )
                           }, 
                           p.grad.fun.cmb = function(vars) {
                             beta_b <- vars[private$bal.idx]
                             return( c(private$b,
                                       -private$marginal.delta,
                                       -sign(beta_b) * private$balance.delta -
                                         private$target.mean)
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
                                       sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           }, 
                           p.fun.cmb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             beta_b <- vars[private$bal.idx]
                             return( sum(private$b * g) -
                                       sum(private$marginal.delta * beta_m) -
                                       sum(abs(beta_b) * private$balance.delta) -
                                       sum(private$target.mean * beta_b )
                             )
                           },
                           h_star_inv.c = function(vars) {
                             eta <- rep(private$vars, each = private$n) - private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cm = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             eta <- rep(g, each = private$n) +  
                               private$marginal.costs %*% beta_m - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             eta <- rep(g, each = private$n) +  
                               private$marginal.costs %*% beta_m - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           },
                           h_star_inv.cmb = function(vars) {
                             g <- vars[private$dual.idx]
                             beta_m <- vars[private$marg.idx]
                             beta_b <- vars[private$bal.idx]
                             eta <- rep(g, each = private$n) +  
                               private$marginal.costs %*% beta_m +
                               Matrix::crossprod(private$balance.functions, beta_b) - 
                               private$cost
                             return(eta * as.numeric(eta > 0) / private$delta)
                           }
                         )
)

otDualL2 <-  R6::R6Class("otDualL2",
                               inherit = optProblem,
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
                                   return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                                 },
                                 get_dist = function(vars) {
                                   f <- self$get_f(vars)
                                   g <- self$get_g(vars)
                                   
                                   return(
                                     sum(private$a * f) + sum(private$b * g)
                                   )
                                 },
                                 init = function() {
                                   rep(max(private$cost), private$nvars)
                                 },
                                 update_a = function(a) {
                                   private$a <- a
                                 },
                                 initialize = function(lambda,
                                                       cost,
                                                       p,
                                                       a,
                                                       b,
                                                       ...
                                 ) {
                                   private$lambda <- lambda
                                   private$cost     <- as.matrix(cost^p)
                                   private$n     <- nrow(private$cost)
                                   private$m     <- ncol(private$cost)
                                   private$fidx <- 1:private$n
                                   private$gidx <- (private$n + 1):(private$n + private$m)
                                   private$a  <- a
                                   private$b  <- b
                                   private$nvars <- length(a) + length(b)
                                 }
                               ),
                               private = list(
                                 a = "numeric",
                                 b = "numeric",
                                 cost = "matrix",
                                 fidx = "integer",
                                 gidx = "integer",
                                 lambda = "numeric",
                                 m = "integer",
                                 n = "integer",
                                 nvars = "integer"
                               )
)

otDualL2_lambda <-  R6::R6Class("otDualL2",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                             f <- self$get_f(vars)
                             g <- self$get_g(vars)
                             private$lambda <- exp(vars[private$nvars])
                             
                             fmat <- matrix(f, private$n, private$m)
                             gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                             diff <- (fmat + gmat - private$cost)
                             
                             obj <- sum(private$a * f) + sum(private$b * g) - 
                               0.5 / private$lambda * sum((diff * (diff > 0))^2)
                             return(-obj)
                           },
                           grad = function(vars) {
                             f <- self$get_f(vars)
                             g <- self$get_g(vars)
                             private$lambda <- exp(vars[private$nvars])
                             
                             fmat <- matrix(f, private$n, private$m)
                             gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                             diff <- (fmat + gmat - private$cost)
                             
                             diff <- diff * (diff > 0)
                             
                             f.grad <- private$a - c(rowSums(diff) / private$lambda)
                             
                             g.grad <- private$b - c(colSums(diff) / private$lambda)
                             
                             return(-c(f.grad, g.grad, 0.5 * sum(diff^2)/private$lambda^2))
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
                             return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                           },
                           get_dist = function(vars) {
                             f <- self$get_f(vars)
                             g <- self$get_g(vars)
                             
                             return(
                               sum(private$a * f) + sum(private$b * g)
                             )
                           },
                           init = function() {
                             c(rep(max(private$cost), private$nvars-1), -5)
                           },
                           update_a = function(a) {
                             private$a <- a
                           },
                           initialize = function(lambda,
                                                 cost,
                                                 p,
                                                 a,
                                                 b,
                                                 ...
                           ) {
                             private$lambda <- lambda
                             private$cost     <- as.matrix(cost^p)
                             private$n     <- nrow(private$cost)
                             private$m     <- ncol(private$cost)
                             private$fidx <- 1:private$n
                             private$gidx <- (private$n + 1):(private$n + private$m)
                             private$a  <- a
                             private$b  <- b
                             private$nvars <- length(a) + length(b) + 1
                           }
                         ),
                         private = list(
                           a = "numeric",
                           b = "numeric",
                           cost = "matrix",
                           fidx = "integer",
                           gidx = "integer",
                           lambda = "numeric",
                           m = "integer",
                           n = "integer",
                           nvars = "integer"
                         )
)

otDualL2_self <-  R6::R6Class("otDualL2",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                             f <- self$get_f(vars)

                             fmat <- matrix(f, private$n, private$n)
                             gmat <- matrix(f, private$n, private$n, byrow = TRUE)
                             
                             diff <- (fmat + gmat - private$cost)
                             
                             obj <- 2 * sum(private$a * f) - 
                               0.5 / private$lambda * sum((diff * (diff > 0))^2)
                             return(-obj)
                           },
                           grad = function(vars) {
                             f <- self$get_f(vars)

                             fmat <- matrix(f, private$n, private$n)
                             gmat <- matrix(f, private$n, private$n, byrow = TRUE)
                             diff <- ( fmat + gmat - private$cost)
                             
                             diff <- diff * (diff > 0)
                             f.grad <- private$a - c(rowSums(diff) / private$lambda)
                             g.grad <- private$a - c(colSums(diff) / private$lambda)
                             
                             return(-(f.grad + g.grad))
                           },
                           h_star_inv = function(vars) {
                             f <- self$get_f(vars)

                             fmat <- matrix(f, private$n, private$m)
                             eta <- (2 * fmat - private$cost)
                             return((eta * (eta > 0))/private$lambda)
                           },
                           get_a = function() {
                             return(private$a)
                           },
                           get_b = function() {
                             return(private$a)
                           },
                           get_f = function(vars) {
                             vars[private$fidx]
                           },
                           get_g = function(vars) {
                             vars[private$fidx]
                           },
                           get_lambda = function() {
                             return(private$lambda)
                           },
                           get_weight = function(vars) {
                             return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                           },
                           get_dist = function(vars) {
                             f <- self$get_f(vars)
                             
                             return(
                               2 * sum(private$a * f)
                             )
                           },
                           init = function() {
                             rep(max(private$cost), private$nvars)
                           },
                           update_a = function(a) {
                             private$a <- a
                           },
                           initialize = function(lambda,
                                                 cost,
                                                 p,
                                                 a,
                                                 b,
                                                 ...
                           ) {
                             private$lambda <- lambda
                             private$cost     <- as.matrix(cost^p)
                             private$n     <- nrow(private$cost)
                             # private$m     <- ncol(private$cost)
                             private$fidx <- 1:private$n
                             # private$gidx <- (private$n + 1):(private$n + private$m)
                             private$a  <- a
                             # private$b  <- b
                             private$nvars <- length(a)
                           }
                         ),
                         private = list(
                           a = "numeric",
                           b = "numeric",
                           cost = "matrix",
                           fidx = "integer",
                           gidx = "integer",
                           lambda = "numeric",
                           m = "integer",
                           n = "integer",
                           nvars = "integer"
                         )
)

otDualL2_self_lambda <-  R6::R6Class("otDualL2",
                              inherit = optProblem,
                              public = list(
                                obj =  function(vars) {
                                  f <- self$get_f(vars)
                                  private$lambda <- exp(vars[length(vars)])
                                  
                                  fmat <- matrix(f, private$n, private$n)
                                  gmat <- matrix(f, private$n, private$n, byrow = TRUE)
                                  
                                  diff <- (fmat + gmat - private$cost)
                                  
                                  obj <- 2 * sum(private$a * f) - 
                                    0.5 / private$lambda * sum((diff * (diff > 0))^2)
                                  return(-obj)
                                },
                                grad = function(vars) {
                                  f <- self$get_f(vars)
                                  private$lambda <- exp(vars[length(vars)])
                                  
                                  fmat <- matrix(f, private$n, private$n)
                                  gmat <- matrix(f, private$n, private$n, byrow = TRUE)
                                  
                                  diff <- ( fmat + gmat - private$cost)
                                  
                                  diff <- diff * (diff > 0)
                                  f.grad <- private$a - c(rowSums(diff) / private$lambda)
                                  g.grad <- private$a - c(colSums(diff) / private$lambda)
                                  
                                  return(-c(f.grad + g.grad, sum(diff^2)/private$lambda^2))
                                },
                                h_star_inv = function(vars) {
                                  f <- self$get_f(vars)

                                  fmat <- matrix(f, private$n, private$m)
                                  eta <- (2 * fmat - private$cost)
                                  return((eta * (eta > 0))/private$lambda)
                                },
                                get_a = function() {
                                  return(private$a)
                                },
                                get_b = function() {
                                  return(private$a)
                                },
                                get_lambda = function() {
                                  return(private$lambda)
                                },
                                get_f = function(vars) {
                                  vars[private$fidx]
                                },
                                get_g = function(vars) {
                                  vars[private$fidx]
                                },
                                get_weight = function(vars) {
                                  return( simplex_proj( as.numeric(self$h_star_inv(vars) ) ))
                                },
                                get_dist = function(vars) {
                                  f <- self$get_f(vars)
                                  
                                  return(
                                    2 * sum(private$a * f)
                                  )
                                },
                                init = function() {
                                  rep(max(private$cost), private$nvars)
                                },
                                update_a = function(a) {
                                  private$a <- a
                                },
                                initialize = function(lambda,
                                                      cost,
                                                      p,
                                                      a,
                                                      b,
                                                      ...
                                ) {
                                  private$lambda <- lambda
                                  private$cost     <- as.matrix(cost^p)
                                  private$n     <- nrow(private$cost)
                                  # private$m     <- ncol(private$cost)
                                  private$fidx <- 1:private$n
                                  # private$gidx <- (private$n + 1):(private$n + private$m)
                                  private$a  <- a
                                  # private$b  <- b
                                  private$nvars <- length(a) + 1
                                }
                              ),
                              private = list(
                                a = "numeric",
                                b = "numeric",
                                cost = "matrix",
                                fidx = "integer",
                                gidx = "integer",
                                lambda = "numeric",
                                m = "integer",
                                n = "integer",
                                nvars = "integer"
                              )
)


otDualEntropy <-  R6::R6Class("otDualEntropy",
                         inherit = optProblem,
                         public = list(
                           obj =  function(vars) {
                             fmat <- matrix(private$f, private$n, private$m)
                             gmat <- matrix(private$g, private$n, private$m, byrow = TRUE)
                             
                             E    <- exp((fmat + gmat - private$cost)/private$lambda) - 1
                             aXb <-  matrix(a, private$n, private$m) * matrix(b, private$n, private$m, byrow = TRUE)
                             
                             obj <- sum(private$a * f) + sum(private$b * g) - 
                               private$lamba * (sum(aXb * E) - sum(aXb))
                             return(-obj)
                           },
                           grad = function(vars) {
                             
                             fmat <- matrix(private$f, private$n, private$m)
                             gmat <- matrix(private$g, private$n, private$m, byrow = TRUE)

                             E    <- exp((fmat + gmat - private$cost)/private$lambda) - 1
                             
                             # TODO: right gradient
                             f.grad <- private$a - c(rowSums(diff) / private$lambda)
                             g.grad <- private$b - c(colSums(diff) / private$lambda)
                             
                             return(-c(f.grad, g.grad))
                           },
                           h_star_inv = function(vars) {
                             f <- self$get_f(vars)
                             g <- self$get_g(vars)
                             
                             fmat <- matrix(f, private$n, private$m)
                             gmat <- matrix(g, private$n, private$m, byrow = TRUE)
                             eta <- (fmat + gmat - private$cost)/private$lambda
                             return((exp(eta)) * matrix(a, private$n, private$m) * matrix(b, private$n, private$m, byrow = TRUE))
                           },
                           get_a = function() {
                             return(private$a)
                           },
                           get_b = function() {
                             return(private$b)
                           },
                           get_f = function(vars) {
                             private$f
                           },
                           get_g = function(vars) {
                             private$g
                           },
                           get_lambda = function() {
                             return(private$lambda)
                           },
                           get_weight = function(vars) {
                             return( as.numeric(self$h_star_inv(vars) ) )
                           },
                           get_dist = function(vars) {
                             
                             return(
                               sum(private$a * private$f) + sum(private$b * prviate$g)
                             )
                           },
                           init = function() {
                             rep(max(private$cost), private$nvars)
                           },
                           update_a = function(a) {
                             private$a <- a
                           },
                           initialize = function(lambda,
                                                 cost,
                                                 p,
                                                 a,
                                                 b,
                                                 x0,
                                                 x1,
                                                 reach = NULL,
                                                 diameter = NULL,
                                                 scaling = 0.5, truncate = 5,
                                                 metric = "sdLp", kernel = NULL,
                                                 cluster_scale=NULL, 
                                                 debias=TRUE, 
                                                 verbose=FALSE, backend='auto',
                                                 ...
                           ) {
                             private$lambda <- lambda
                             private$cost     <- as.matrix(cost^p)
                             private$n     <- nrow(private$cost)
                             private$m     <- ncol(private$cost)
                             private$fidx <- 1:private$n
                             private$gidx <- (private$n + 1):(private$n + private$m)
                             private$a  <- a
                             private$b  <- b
                             private$nvars <- length(a) + length(b)
                             private$x <- x0
                             private$y <- x1
                             private$sinkhorn_args <- list(p = p,
                                                           reach = reach,
                                                           diameter = diameter,
                                                           scaling = scaling,
                                                           truncate = truncate,
                                                           metric = metric,
                                                           kernel = kernel,
                                                           cluster_scale = cluster_scale,
                                                           debias = debias,
                                                           verbose = verbose,
                                                           backend = backend
                                                           )
                           },
                           solve = function(){
                             sol <- sinkhorn_geom(x, y, a, b, power = private$sinkhorn_args$p, 
                                           blur = private$lambda, reach = private$sinkhorn_args$reach, 
                                           diameter = private$sinkhorn_args$diameter,
                                           scaling = private$sinkhorn_args$scaling, 
                                           truncate = private$sinkhorn_args$truncate,
                                           metric = "sdLp", kernel = private$sinkhorn_args$kernel,
                                           cluster_scale=private$sinkhorn_args$cluster_scale, 
                                           debias=private$sinkhorn_args$debias, 
                                           verbose=private$sinkhorn_args$verbose, 
                                           backend=private$sinkhorn_args$backend)
                             private$f <- sol$f
                             private$g <- sol$g
                           }
                         ),
                         private = list(
                           a = "numeric",
                           b = "numeric",
                           cost = "matrix",
                           f = "numeric",
                           g = "numeric",
                           lambda = "numeric",
                           m = "integer",
                           n = "integer",
                           nvars = "integer",
                           sinkhorn_args = "list",
                           x = "matrix",
                           y = "matrix"
                         )
)  

  
otDualOpt <- function(x, z, p, metric, lambda, penalty = "L2", optimizer = NULL, sample_weight = NULL, self = FALSE, control = NULL, ...) {
  
  if (is.null(optimizer)) {
    sw <- get_sample_weight(sample_weight, z)
    a <- sw$a
    b <- sw$b
    cost <- cost_fun(x = x, z = z, p = p, metric = metric, estimand = "ATT")
    if (!self ) {
      optimizer <- switch(penalty,
                        "L2" = otDualL2$new(lambda,
                                              cost,
                                              p,
                                              a,
                                              b),
                        "entropy" = otDualEntropy$init(lambda = lambda,
                                                       cost = cost,
                                                       a = a,
                                                       b = b,
                                                       p = p,
                                                       x0 = x[z==0,,drop = FALSE],
                                                       x1 = x[z==1,,drop = FALSE],
                                                       metric = metric,
                                                       ...))
    } else {
      optimizer <- switch(penalty,
                          "L2" = otDualL2_self$new(lambda,
                                              cost,
                                              p,
                                              a),
                          "entropy" = otDualEntropy$init())
    }
  }
  
  if (penalty == "L2") {
    fit <- lbfgsb3c::lbfgsb3(par = optimizer$init(),
                             fn = optimizer$obj,
                             gr = optimizer$grad,
                             lower = -Inf,
                             upper = Inf,
                             control = control
    )
    f <- optimizer$get_f(fit$par)
    g <- optimizer$get_g(fit$par)
    lambda <- optimizer$get_lambda()
  } else {
    optimizer$solve()
    f <- optimizer$get_f()
    g <- optimizer$get_g()
    lambda <- optimizer$get_lambda()
    fit <- NULL
  }
  
  return(list(f = f, g = g, lambda =  lambda, optimizer = optimizer, fit = fit))
  
}
  
cotDualL2_2 <- R6::R6Class("cotDualL2",
                           inherit = optProblem,
                           public = list(
                             obj =  function(vars) {
                               diff <- (private$Q %*% vars - private$cost)
                               diff <- diff * (diff > 0)
                               return(-(private$p.fun(vars) - 0.5 * sum(diff^2) / private$lambda))
                             },
                             grad = function(vars) {
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
                                 rbind(bind(rep(-Inf, private$m ),
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
                               diff <- (private$Q %*% vars - private$cost)
                               pos <- (diff > 0)
                               return(Matrix::crossprod(private$Q * pos))
                             },
                             get_xty = function() {
                               neg.mass.const <- -private$Q[,private$dual.idx]
                               cmb.Q <- cbind(neg.mass.const,
                                              private$Q)
                               QtC  <- as.numeric(Matrix::crossprod(cmb.Q, 
                                                                    private$cost)) / private$delta
                               l_t  <- length(private$target.mean)
                               A    <- c(rep(0, length(QtC) - l_t), private$target.mean)
                               return(QtC + A)
                             },
                             get_hessian = function(vars) {
                               diff <- (private$Q %*% vars - private$cost)
                               pos <- as(diff > 0, "sparseVector")
                               return(-Matrix::crossprod(private$Q * pos)/private$lambda)
                             },
                             get_grad_var = function(vars) {
                               # ignore negative sign and penalty functions since variance of constants is 0 and negative sign is squared
                               diff <- (private$Q %*% vars - private$cost)
                               diff  <- diff * (diff > 0)
                               return(Matrix::crossprod(private$Q * as(diff, "sparseVector"))/private$lambda)
                             },
                             # get_X = function() {
                             #   neg.mass.const <- -private$Q[,private$dual.idx]
                             #   cmb.Q <- cbind(neg.mass.const,
                             #                  private$Q)
                             #   A <- cbind(
                             #     matrix(0, nrow = length(private$target.mean),
                             #          ncol = ncol(cmb.Q) - length(private$target.mean)),
                             #     Matrix::Diagonal(length(private$target.mean),
                             #                      x = private$target.mean))
                             #   return(rbind(cmb.Q / sqrt(private$delta),
                             #                 A) )
                             # },
                             # get_Y = function() {
                             #   return(c(private$cost / sqrt(private$delta), 
                             #            rep(1, length(private$target.mean))))
                             # },
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
                                                    "1" = Matrix::t(vec_to_col_constraints(private$n,
                                                                                           private$m)),
                                                    "2" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                 private$m)),
                                                                -private$marginal.costs),
                                                    "3" = cbind(Matrix::t(vec_to_col_constraints(private$n,
                                                                                                 private$m)),
                                                                -private$balance.functions,
                                                                private$balance.functions),
                                                    "4" = cbind(Matrix::t(vec_to_col_constraints(private$n,
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


cot_dual_opt <- function(x, target, 
                         init = NULL,
                         sample_weights = NULL, 
                         method = c("SBW", "Wasserstein"),
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
  
  fit <- lbfgsb3c::lbfgsb3(par = init,
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
    weight <- Matrix::Matrix(data = weight,
                             nrow = n, ncol = m)
  }
  
  return(list(weight = weight, beta = beta,
              # unconstrained = list(weight = unconst_weight),
              fit = fit,
              optimizer = optimizer))
}

