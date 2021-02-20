# conditional gradient descent
cg <- function(optimizer, verbose = TRUE) {
  stopifnot(inherits(optimizer, "cgOptimizer")) #needs to be R6 method
  
  optimizer$solve_G()
  if (verbose) pb <- txtProgressBar(min = 0, max = floor(optimizer$get_niter()/10), style = 3)
  for (i in 1:optimizer$get_niter()) {
    optimizer$solve_param()
    optimizer$solve_S()
    optimizer$step(i)
    if (optimizer$converged()) break
    if (verbose && i %% 10 == 0) setTxtProgressBar(pb, i/10)
  }
  if (verbose) close(pb)
  return(invisible(optimizer))
}

cgOptimizer <- R6::R6Class("cgOptimizer", 
                          public = list(
                            cold_start = function(lambda) {
                              self$warm_start(lambda)
                              private$qp$obj$L[private$cost_idx] <- c(private$cost)
                              private$run_warm_start <- FALSE
                            },
                            converged = function() {
                              private$f_val <- self$f()
                              f_val_diff <- abs(private$f_val - private$f_val_old)
                              # diff_param <- (private$param - private$param_old)
                              # paramnorm <-  sum(diff_param^2)/sum(private$param_old^2)
                              # if(is.nan(paramnorm) & sum(diff_param^2) == 0 ) paramnorm <- 0
                              # diff_G <- (private$G - private$G_old)
                              # Gnorm <- sum(diff_G^2)/sum(private$G_old^2)
                              # if(is.nan(Gnorm) & sum(diff_G^2) == 0 ) Gnorm <- 0
                              
                              converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol) 
                              # ||
                              #   isTRUE(f_val_diff < private$tol) 
                              # || (isTRUE( sqrt(paramnorm) < private$tol) &&
                              #   isTRUE( sqrt(Gnorm) < private$tol) )
                              # private$G_old <- private$G
                              # private$param_old <- private$param
                              private$f_val_old <- private$f_val
                              return(converged)
                            },
                            f = function() {
                              return(private$f_user(private$X1, private$X2, 
                                                    private$Y1, private$Y2, 
                                                    private$cost, 
                                                    private$G, 
                                                    private$param, 
                                                    private$b, private$p,
                                                    private$lambda))
                            },
                            df = function() {
                              return(private$df_user(private$X1, private$X2, 
                                                     private$Y1, private$Y2, 
                                                    private$cost, 
                                                    private$G, 
                                                    private$param,
                                                    private$b, private$p,
                                                    private$lambda))
                            },
                            get_G = function() {
                              return(private$G)
                            },
                            get_param = function() {
                              return(private$param)
                            },
                            get_niter = function() {
                              return(private$niter)
                            },
                            return_cw = function() {
                              out <- list(w0 = rowSums(private$G),
                                          w1 = colSums(private$G),
                                          gamma = private$G,
                                          estimand = "ATT",
                                          method = "Wasserstein",
                                          args = list(
                                                      constraint = list(joint = private$lambda[1],
                                                                        penalty = private$lambda[3],
                                                                        margins = private$penalty_list$margins),
                                                      power = private$power,
                                                      metric = private$metric,
                                                      niter = private$niter,
                                                      cur_iter = private$cur_iter,
                                                      add.margins = private$add.margins,
                                                      joint.mapping = private$add.mapping,
                                                      conditional.gradient = TRUE))
                              class(out) <- "causalWeights"
                              return(out)
                            },
                            solve_G = function() {
                              if (!private$run_warm_start) {
                                private$G <- matrix(private$solver(private$qp), private$n1, private$n2)
                              }
                              # private$G_old <- private$G
                              # private$param_old <- private$param
                              private$f_val_old <- self$f()
                            },
                            solve_param = function() {
                              return( private$solve_param_user(private$X1, private$X2, 
                                                               private$Y1, private$Y2, 
                                                                private$cost, private$G,
                                                                private$param,
                                                                private$lambda) )
                            },
                            solve_S = function() {
                              private$qp$obj$L[private$cost_idx] <- c(self$df())
                              private$S <- private$solver(private$qp)
                            },
                            step = function(i) {
                              private$cur_iter <- i
                              deltaG <- private$S - private$G
                              f_val <- self$f()
                              df_val <- private$qp$obj$L[private$cost_idx]
                              
                              f <- function(x, ...) {
                                g <- matrix(x[1:length(private$G)], private$n1, private$n2)
                                private$f_user(private$X1, private$X2, 
                                                private$Y1, private$Y2, 
                                                private$cost, 
                                                g, 
                                                private$param,
                                                private$b, private$p,
                                                private$lambda)
                              }
                              # search_res <- rje::armijo(
                              #   fun = f,
                              #   x = private$G,
                              #   dx = deltaG,
                              #   beta = 3,
                              #   sigma = 1e-10,
                              #   grad = df_val,
                              #   maximise = FALSE,
                              #   searchup = TRUE,
                              #   adj.start = 1
                              # )
                              # private$G <- private$G + search_res$adj * deltaG
                              derphi0 <- sum(deltaG * df_val)
                              if (-derphi0 >  private$tol) {
                                search_res <-  scalar_search_armijo(phi = f,
                                                                    phi0 = f_val,
                                                                    derphi0 = derphi0,
                                                                    x = c(private$G),
                                                                    dx = c(deltaG),
                                                                    c1 = 1e-4, alpha0 = .99, amin = 0)
                                if ( !is.null(search_res$alpha)) {
                                  search_res$alpha = min(1, max(search_res$alpha,0))
                                  private$G <- private$G + search_res$alpha * deltaG
                                }
                                
                              }
                              
                              
                            }, #function(iter) return(2/(iter + 1)),
                            warm_start = function(lambda = list(joint = NULL, parameters = NULL, penalty = NULL, margins = NULL)) {
                              private$lambda <- c(lambda$joint, 0, lambda$penalty)
                              private$qp$obj$L[grep("pen", names(private$qp$obj$L))] <- lambda$penalty
                              private$qp$LC$vals[grep("margins", names(private$qp$LC$vals))] <- lambda$margins
                              private$run_warm_start <- TRUE
                            }
                          ),
                          private = list(
                            "a" = "numeric",
                            "add.mapping" = "logical",
                            "add.margins" = "logical",
                            "b" = "numeric",
                            "cost" = "numeric",
                            "cost_idx" = "numeric",
                            "cur_iter" = "integer",
                            "df_user" = "function",
                            "f_user" = "function",
                            "f_val" = "numeric",
                            "f_val_old" = "numeric",
                            "G" = "numeric",
                            "G_old" = "numeric",
                            "lambda" = "numeric",
                            "n1" = "numeric",
                            "n2" = "numeric",
                            "niter" = "numeric",
                            "p" = "numeric",
                            "param" = "numeric",
                            "param_old" = "numeric",
                            "penalty" = "character",
                            "penalty_list" = "list",
                            "qp" = "list",
                            "run_warm_start" = "logical",
                            "S" = "numeric",
                            "solve_param_user" = "function",
                            "solver" = "function",
                            "tol" = "numeric",
                            "X1" = "matrix",
                            "X2" = "matrix",
                            "Y1" = "matrix",
                            "Y2" = "matrix"
                          )
                          
)

wassCGD <- R6::R6Class("wassCGD",
  inherit = cgOptimizer,
  public = list(
    initialize = function(X1, X2, 
                          cost,
                          qp_constraint = list(joint = 1,
                                               margins = 1,
                                               penalty = 1),
                          qp_solver = c("mosek", "gurobi", 
                                        "cplex"),
                          qp_penalty = "L2",
                          lambda = list(joint = 1,
                                        penalty = 1
                                        ),
                          add.mapping = FALSE,
                          add.margins = FALSE,
                          penalty = "L2",
                          metric = dist.metrics(),
                          power = 2,
                          niter = 1000,
                          tol = 1e-7,
                          sample_weight = NULL,
                          ...
    ) {
      metric    <- match.arg(metric, dist.metrics())
      private$add.margins <- isTRUE(add.margins)
      private$add.mapping <- isTRUE(add.mapping)
      private$penalty <- match.arg(penalty, c("L2","none", "entropy",
                                              "variance"))
      qp_penalty <- match.arg(qp_penalty, c("L2","none", "entropy",
                                            "variance"))
      if (missing(niter) || length(niter) == 0) niter <- 1000
      if (missing(tol) || length(tol) == 0) tol <= 1e-7 
      
      
      private$n1 <- nrow(X1)
      private$n2 <- nrow(X2)
      
      x <- rbind(X1,X2)
      z <- c(rep(0, private$n1), rep(1, private$n2))
      
      sw  <- get_sample_weight(sample_weight, z)
      private$a <- sw$a
      private$b <- sw$b
      
      if (missing(cost) || length(cost) == 0) {
        private$cost <- cost_fun(x = x, z = z, power = power, metric = metric,
                                 estimand = "ATT")^power
        if (add.margins) {
          costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
                                                              estimand = "ATT")),
                      list(private$cost^(1/power)))
        } else {
          costqp <- private$cost^(1/power)
        }
      } else {
        if (add.margins) {
          costqp <- cost
          private$cost   <- cost[[length(cost)]]^(power)
        } else {
          private$cost <- cost^(power)
        }
      }
      
      if (add.mapping) private$cost <- private$cost / max(private$cost)
      
      private$penalty_list <- list(margins = qp_constraint$margins)
      
      private$qp <- qp_wass(x, z, K = qp_constraint, p = power,
                                                   penalty = qp_penalty,
                                                   cost = costqp,
                                                   add.margins = add.margins,
                                                   add.joint = TRUE, joint.mapping = FALSE,
                                                   neg.weights = FALSE, bf = NULL,
                                                   sample_weight = NULL, soc = FALSE
                           )
                           
      stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
      # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
      private$cost_idx <- grep("cost", names(private$qp$obj$L))
      solver <- switch(qp_solver,
                       "cplex" = cplex_solver,
                       "gurobi" = gurobi_solver,
                       "mosek" = mosek_solver)
      private$solver <- function(qp){
        sol <- solver(qp)
        return(matrix(sol, private$n1, private$n2))
      }
      private$tol <- tol
      private$niter <- as.numeric(niter)
      
      f.df.fun <- if (add.mapping) {
          jm_nooutcome(method = "Wasserstein", penalty = private$penalty)
        } else {
          nojm_nooutcome(method = "Wasserstein", penalty = private$penalty)
        }
      
      private$f_user <- f.df.fun$f_user
      #   function(X1, X2,
      #                            Y1, Y2,
      #                            cost, G,
      #                            param,
      #                            lambda) {
      #   sum((X2 - sweep(crossprod(G, X1), MARGIN = 1, FUN = "/", STAT = private$b,
      #                   check.margin = FALSE))^2 * private$b / (private$p)) + 
      #     lambda[1] * sum(cost * G) +
      #     lambda[2] * sum(param^2)/private$p^2 +
      #     lambda[3] * 0.5 * sum(G^2)
      #   
      # }
      private$df_user <-  f.df.fun$df_user
      #   function(X1, X2,
      #                              Y1, Y2,
      #                              cost, G,
      #                              param,
      #                              lambda) {
      #   2 * tcrossprod(X1, X2 - sweep(crossprod(G, X1), MARGIN = 1, FUN = "/", STAT = private$b,
      #                                 check.margin = FALSE))  * private$b / (private$p) + 
      #     lambda[1] * cost +
      #     lambda[3] * G
      # }
      private$solve_param_user <- function(...) {
      }
      
      private$X1 <- X1
      private$X2 <- X2
      
      private$p  <- ncol(X1)
      private$n1 <- nrow(X1)
      private$n2 <- nrow(X2)
      
      private$param <- rep(0, private$p)
      
      private$lambda <- c(cost = lambda$joint, coefficients = 0.0,
                          gamma = lambda$penalty)
      
      private$run_warm_start <- FALSE
      
    }
  )
                                              
                                              
)

jointMapAndWeights <- R6::R6Class("jointMapAndWeights",
                              inherit = cgOptimizer,
                              public = list(
                                initialize = function(data,
                                                      map = c("linear", "gp"),
                                                      qp_method = c("Wasserstein",
                                                                    "Constrained Wasserstein"),
                                                      qp_constraint = list(joint = 1,
                                                                           margins = 1),
                                                      lambda_c = 1,
                                                      lambda_b = 1,
                                                      lambda_g = 1,
                                                      metric = dist.metrics(),
                                                      power = 2,
                                                      qp_solver = c("mosek", "gurobi", 
                                                                    "cplex"),
                                                      niter = 1000,
                                                      tol = 1e-7,
                                                      sample_weight = NULL,
                                                      penalty = "L2",
                                                      add.margins = FALSE,
                                                      ...
                                ) {
                                  qp_method <- match.arg(qp_method)
                                  metric    <- match.arg(metric, dist.metrics())
                                  private$add.margins <- isTRUE(add.margins)
                                  private$penalty <- match.arg(penalty, c("L2","none", "entropy",
                                                                          "variane"))
                                  if (length(niter) == 0) niter <- 1000
                                  if (length(tol) == 0) tol <= 1e-7 
                                  
                                  pd  <- prep_data(data, ...)
                                  
                                  y <- pd$df$y
                                  x <- as.matrix(pd$df[,-match("y", colnames(pd$df)), drop = FALSE])
                                  z <- pd$z
                                  
                                  sw  <- get_sample_weight(sample_weight, z)
                                  private$a <- sw$a
                                  private$b <- sw$b
                                  
                                  X1 <- x[z == 0, , drop = FALSE]
                                  X2 <- x[z == 1, , drop = FALSE]
                                  
                                  if (!is.null(y)) {
                                    Y1 <- y[z == 0]
                                    Y2 <- y[z == 1]
                                  } else {
                                    Y1 <- Y2 <- NULL
                                  }
                                  
                                  private$n1 <- nrow(X1)
                                  private$n2 <- nrow(X2)
                                  
                                  private$cost <- cost_fun(x = x, z = z, power = power, metric = metric,
                                                           estimand = "ATT")^power
                                  private$cost <- private$cost / max(private$cost)
                                  
                                  if (add.margins) {
                                    costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
                                                                      estimand = "ATT")^power),
                                                      list(private$cost))
                                  } else {
                                    costqp <- private$cost
                                  }
                                  
                                  private$qp <- switch(qp_method,
                                                       "Wasserstein" = qp_wass(x, z, K = qp_constraint, p = power,
                                                                               penalty = penalty,
                                                                               cost = costqp,
                                                                               add.margins = add.margins,
                                                                               add.joint = TRUE, joint.mapping = FALSE,
                                                                               neg.weights = FALSE, bf = NULL,
                                                                               sample_weight = NULL, soc = FALSE
                                                                               ),
                                                       "Constrained Wasserstein" = qp_wass_const(x, z, K = qp_constraint, p = power,
                                                                                           penalty = penalty,
                                                                                           cost = private$cost,
                                                                                           add.margins = add.margins,
                                                                                           add.joint = TRUE, joint.mapping = FALSE,
                                                                                           neg.weights = FALSE, bf = NULL,
                                                                                           sample_weight = NULL, soc = FALSE
                                                       )
                                                       
                                  )
                                  stopifnot(all(names(private$qp) %in% c("obj", "LC","bounds","nvar")))
                                  # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
                                  
                                  solver <- switch(qp_solver,
                                                "cplex" = cplex_solver,
                                                "gurobi" = gurobi_solver,
                                                "mosek" = mosek_solver)
                                  private$solver <- function(qp){
                                    sol <- solver(qp)
                                    return(matrix(sol, private$n1, private$n2))
                                  }
                                  private$tol <- tol
                                  private$niter <- as.numeric(niter)
                                  
                                  private$f_user <- function(X1, X2,
                                                             Y1, Y2,
                                                             cost, G,
                                                             param,
                                                             b, d,
                                                             lambda) {
                                    sum((X2 - sweep(crossprod(G, X1), MARGIN = 1, FUN = "/", STAT = private$b,
                                                    check.margin = FALSE))^2 * private$b / (private$p)) + 
                                      lambda[1] * sum(cost * G) +
                                      lambda[2] * sum(param^2)/private$p^2 +
                                      lambda[3] * 0.5 * sum(G^2)
                                    
                                  }
                                  private$df_user <-  function(X1, X2,
                                                               Y1, Y2,
                                                               cost, G,
                                                               param,
                                                               b, d,
                                                               lambda) {
                                    2 * tcrossprod(X1, X2 - sweep(crossprod(G, X1), MARGIN = 1, FUN = "/", STAT = private$b,
                                                                  check.margin = FALSE))  * private$b / (private$p) + 
                                      lambda[1] * cost +
                                    lambda[3] * G
                                  }
                                  private$solve_param_user <- function(...) {
                                  }
                                  
                                  private$X1 <- X1
                                  private$X2 <- X2
                                  
                                  
                                  private$Y1 <- Y1
                                  private$Y2 <- Y2
                                  
                                  private$p  <- ncol(X1)
                                  private$n1 <- nrow(X1)
                                  private$n2 <- nrow(X2)
                                  
                                  private$param <- rep(0, private$p)
                                  
                                  private$lambda <- c(cost = lambda_c, coefficients = lambda_b,
                                                      gamma = lambda_g)
                                  
                                  private$run_warm_start <- FALSE
                                  
                                }
                              )
                              
                                
)

nojm_nooutcome <- function(method, penalty) {
  out <- list()
  pen.fun <- switch(penalty,
                    "none" = function(G){return(0.0)},
                    "L2" = function(G) {return(0.5 * sum(G^2))},
                    "variance" = function(G) {return(0.5 * sum(rowSums(G)^2))},
                    "entropy" = function(G) {
                      posG <- G[G > 0]
                      return(sum(posG * log(posG)))}
  )
  dpen.fun <- switch(penalty,
                     "none" = function(G){return(0.0)},
                     "L2" = function(G) {return(G)},
                     "variance" = function(G) {return(matrix(rowSums(G), nrow(G), ncol(G)))},
                     "entropy" = function(G) {
                       l_G <- log(G)
                       l_G[is.infinite(l_G)] <- -1e6
                       return(matrix(1, nrow(G), ncol(G)) + 
                                l_G                        
                       )
                     })
  if (method == "Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           b, d,
                           lambda) {
      sum(cost * G) +
        lambda[3] * pen.fun(G)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             b, d,
                             lambda) {
        cost +
        lambda[3] * dpen.fun(G)
    }
    out$solve_param_user <- function() {
      
    }
  } else if (method == "Constrained Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           b, d,
                           lambda
                           ) {
      lambda[3] * pen.fun(G)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             b, d,
                             lambda) {
        lambda[3] * dpen.fun(G)
    }
    out$solve_param_user <- function() {
      
    }
  }
  
  return(out)
}


jm_nooutcome <- function(method, penalty) {
  out <- list()
  pen.fun <- switch(penalty,
                    "none" = function(G){return(0.0)},
                    "L2" = function(G) {return(0.5 * sum(G^2))},
                    "variance" = function(G) {return(0.5 * sum(rowSums(G)^2))},
                    "entropy" = function(G) {
                      posG <- G[G > 0]
                      return(sum(posG * log(posG)))}
  )
  dpen.fun <- switch(penalty,
                     "none" = function(G){return(0.0)},
                     "L2" = function(G) {return(G)},
                     "variance" = function(G) {return(matrix(rowSums(G), nrow(G), ncol(G)))},
                     "entropy" = function(G) {
                       l_G <- log(G)
                       l_G[is.infinite(l_G)] <- -1e6
                       return(matrix(1, nrow(G), ncol(G)) + 
                                                       l_G                        
                     )
                     })
  if (method == "Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           b, d,
                           lambda) {
      sum((X2 - crossprod(G, X1) / b)^2 * b / d) + 
        lambda[1] * sum(cost * G) +
        lambda[3] * pen.fun(G)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             b, d,
                             lambda) {
      2 * tcrossprod(X1, X2 - crossprod(G, X1) / b)  * b / d + 
        lambda[1] * cost +
        lambda[3] * dpen.fun(G)
    }
    out$solve_param_user <- function() {
      
    }
  } else if (method == "Constrained Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           b, d,
                           lambda
    ) {
      sum((X2 - crossprod(G, X1) / b)^2 * b / d) + 
      lambda[3] * pen.fun(G)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             b, d,
                             lambda) {
      2 * tcrossprod(X1, X2 - crossprod(G, X1) / b)  * b / d + 
      lambda[3] * dpen.fun(G)
    }
    out$solve_param_user <- function() {
      
    }
  }
  
  return(out)
}

jm_outcome_linear <- function(method) {
  out <- list()
  if (method == "Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           lambda,
                           b, d) {
      sum((X2 - crossprod(G, X1) / b)^2) * b / (d) + 
      sum((X2 %*% param[-1] + param[1] - crossprod(G, Y1) )^2) +
        lambda[1] * sum(cost * G) +
        lambda[2] * sum(param^2)/d^2
      lambda[3] * 0.5 * sum(G^2)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             lambda,
                             b, d) {
      2 * tcrossprod(X1, X2 - crossprod(G, X1) / b)  * b / (d) + 
      2 * tcrossprod(Y1, X2 %*% param[-1] + param[1] - crossprod(G, Y1) / b)  * b / (d) +
        lambda[1] * cost +
        lambda[3] * G
    }
    
  } else if (method == "Constrained Wasserstein") {
    out$f_user <- function(X1, X2,
                           Y1, Y2,
                           cost, G,
                           param,
                           lambda,
                           b, d) {
      sum((X2 - crossprod(G, X1) / b)^2) * b / (d) + 
        sum((X2 %*% param[-1] + param[1] - crossprod(G, Y1) )^2) +
        lambda[2] * sum(param^2) / d^2
        lambda[3] * sum(G^2)
      
    }
    out$df_user <-  function(X1, X2,
                             Y1, Y2,
                             cost, G,
                             param,
                             lambda,
                             b, d) {
      2 * tcrossprod(X1, X2 - crossprod(G, X1) / b)  * b / (d) + 
        2 * tcrossprod(Y1, X2 %*% param[-1] + param[1] - crossprod(G, Y1) / b)  * b / (d) +
        lambda[3] * G
    }
  }
  out$solve_param_user <- function(X1, X2,
                                   Y1, Y2,
                                   cost, G,
                                   param,
                                   lambda,
                                   b, d) {
    inv_mat(crossprod(X2) + diag(lambda[2] / (d^2), d, d) ) %*% 
      crossprod(X2, crossprod(G, Y1)) / b
  }
  
  return(out)
}

#based on scipy implementaitaon
scalar_search_armijo <- function(phi, phi0, derphi0, x, dx, c1=1e-4, alpha0=1, amin=0) {
  #   Minimize over alpha, the function ``phi(alpha)``.
  #   Uses the interpolation algorithm (Armijo backtracking) as suggested by
  #   Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57
  #   alpha > 0 is assumed to be a descent direction.
  #   Returns
  #   -------
  #   alpha
  #   phi1
  
  phi_a0 = phi(x + dx * alpha0)
  if (phi_a0 <= phi0 + c1 * alpha0 * derphi0) {
    return(list(alpha = alpha0, ph1 = phi_a0))
  }
  
  # Otherwise, compute the minimizer of a quadratic interpolant:
  alpha1 = -(derphi0) * alpha0^2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
  phi_a1 = phi(x + dx * alpha1)
  if ((phi_a1 <= (phi0 + c1 * alpha1 * derphi0) ) ) { #&& (alpha1 >= 0) ) {
    return(list(alpha = alpha1, phi1 = phi_a1))
  }
  if (alpha1 < 0) alpha1 <- alpha0 - 0.01  #avoids the negative step size
  # Otherwise, loop with cubic interpolation until we find an alpha which
  # satisfies the first Wolfe condition (since we are backtracking, we will
  # assume that the value of alpha is not too small and satisfies the second
  # condition.
  a <- b <- alpha2 <- phi_a2 <- NULL
  while (alpha1 > amin) {      # we are assuming alpha>0 is a descent direction
    factor = alpha0^2 * alpha1^2 * (alpha1 - alpha0)
    a = alpha0^2 * (phi_a1 - phi0 - derphi0 * alpha1) - 
            alpha1^2 * (phi_a0 - phi0 - derphi0 * alpha0)
    a = a / factor
    b = -alpha0^3 * (phi_a1 - phi0 - derphi0 * alpha1) + 
            alpha1^3 * (phi_a0 - phi0 - derphi0 * alpha0)
    b = b / factor
    
    alpha2 = (-b + sqrt(abs(b^2 - 3 * a * derphi0))) / (3.0 * a)
    phi_a2 = phi(x + dx * alpha2)
    if ((phi_a2 <= (phi0 + c1 * alpha2 * derphi0)) ) {# && alpha2 >= 0) {
      return(list(alpha = alpha2, phi1 = phi_a2))
    }
    
    if ((alpha1 - alpha2) > alpha1 / 2.0 || (1 - alpha2/alpha1) < 0.96) {
      alpha2 = alpha1 / 2.0
      alpha0 = alpha1
      alpha1 = alpha2
      phi_a0 = phi_a1
      phi_a1 = phi_a2
    }
  }
  
  # ## brute force
  # phi_a3 <- NULL
  # for (alpha3 in seq(alpha0, amin, by = -0.01)) {
  #   phi_a3 = phi(x + dx * alpha3)
  #   if (phi_a3 <= (phi0 + c1 * alpha3 * derphi0) ) {
  #     return(list(alpha = alpha3, phi1 = phi_a3))
  #   }
  # }
    
  # Failed to find a suitable step length
  return(list(alpha = NULL, phi1 = phi_a1))
}
  

