# conditional gradient descent
cg <- function(optimizer, verbose = TRUE) {
  stopifnot(inherits(optimizer, "cgOptimizer")) #needs to be R6 method
  
  optimizer$solve_G()
  if (verbose) pb <- txtProgressBar(min = 0, max = floor(optimizer$get_niter()/10), style = 3)
  for (i in 1:optimizer$get_niter()) {
    optimizer$solve_param()
    if (optimizer$converged() && i > 1) break
    optimizer$solve_S()
    optimizer$step()
    if (verbose && i %% 10 == 0) setTxtProgressBar(pb, i/10)
  }
  if (verbose) close(pb)
  # optimizer$solve_param()
  return(invisible(optimizer))
}

#parent class
{
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
                               step = function() {
                                 private$cur_iter <- private$cur_iter + 1L
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
}

# divergence parent class
{
  wassDiv <- R6::R6Class("wassDiv", 
                            inherit = cgOptimizer,
                            public = list(
                              converged = function() {
                                private$f_val <- sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
                                # print(private$f_val)
                                f_val_diff <- abs(private$f_val - private$f_val_old)
                                # if ((i %% 100) == 0) message(round(sum(private$a * y[z==0]), digits = 3), ", ", appendLF = FALSE)
                                # message(round(private$a * y[z==0]))
                                # message(private$f_val, ", ", appendLF = FALSE)
                                # diff_param <- (private$param - private$param_old)
                                # paramnorm <-  sum(diff_param^2)/sum(private$param_old^2)
                                # if(is.nan(paramnorm) & sum(diff_param^2) == 0 ) paramnorm <- 0
                                # diff_G <- (private$G - private$G_old)
                                # Gnorm <- sum(diff_G^2)/sum(private$G_old^2)
                                # if(is.nan(Gnorm) & sum(diff_G^2) == 0 ) Gnorm <- 0
                                nan.check <- isTRUE(any(is.nan(private$a)) || is.nan(private$f_val))
                                
                                converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol)  ||
                                  isTRUE(sum(abs(private$a_old - private$a)) < private$tol) ||
                                  nan.check
                                # message(sum(abs(private$a_old - private$a)), ", ", appendLF = FALSE)
                                # ||
                                #   isTRUE(f_val_diff < private$tol) 
                                # || (isTRUE( sqrt(paramnorm) < private$tol) &&
                                #   isTRUE( sqrt(Gnorm) < private$tol) )
                                # private$G_old <- private$G
                                # private$param_old <- private$param
                                if (nan.check) warning("NaN found in parameters! Try reducing the step size.")
                                private$f_val_old <- private$f_val
                                private$a_old <- private$a
                                return(converged)
                              },
                              f = function() {
                                return(
                                  sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
                                )
                              },
                              df = function() {
                                return(private$f_pot)
                      
                              },
                              get_G = function() {
                                return(private$a)
                              },
                              get_param = function() {
                                return(list(f = private$f_pot,
                                            g = private$g_pot))
                              },
                              get_niter = function() {
                                return(private$niter)
                              },
                              get_weight = function() {
                                NULL
                              },
                              return_cw = function() {
                                out <- list(w0 = private$a,
                                            w1 = private$b,
                                            gamma = self$get_weight(),
                                            estimand = "ATT",
                                            method = "Wasserstein",
                                            args = list(
                                              constraint = list(joint = NULL,
                                                                penalty = private$lambda,
                                                                margins = NULL),
                                              penalty = "entropy",
                                              power = private$power,
                                              metric = private$metric,
                                              niter = private$niter,
                                              cur_iter = private$cur_iter,
                                              add.margins = private$add.margins,
                                              joint.mapping = FALSE,
                                              add.divergence = TRUE,
                                              conditional.gradient = TRUE))
                                class(out) <- "causalWeights"
                                return(out)
                              },
                              solve_G = function() {
                                self$solve_param()
                                private$f_val_old <- self$f()
                                private$a_old <- private$a
                              },
                              solve_param = function() {
                                
                                sol <- sinkhorn_geom(x = private$X1, y = private$X2, 
                                                     a = private$a, 
                                                     b = private$b, power = private$p, 
                                                     blur = private$lambda, reach = private$sinkhorn_args$reach, 
                                                     diameter = private$sinkhorn_args$diameter,
                                                     scaling = private$sinkhorn_args$scaling, 
                                                     truncate = private$sinkhorn_args$truncate,
                                                     metric = "Lp", kernel = private$sinkhorn_args$kernel,
                                                     cluster_scale=private$sinkhorn_args$cluster_scale, 
                                                     debias=private$sinkhorn_args$debias, 
                                                     verbose=private$sinkhorn_args$verbose, 
                                                     backend=private$sinkhorn_args$backend)
                                private$f_pot <- sol$f
                                private$g_pot <- sol$g
                                
                              },
                              solve_S = function() {
                                if (private$search != "mirror-nocgd" && 
                                    private$search != "mirror-accelerated") {
                                  private$qp <- qp_dual(f = private$f_pot,
                                                        g = private$g_pot,
                                                        b = private$b)
                                  private$S <- private$solver(private$qp)$sol
                                  
                                  # private$S <- as.numeric(private$f_pot == min(private$f_pot))
                                  
                                  # grad <- self$df()
                                  # private$S <- renormalize(as.numeric(grad == min(grad)))
                                }
                                
                              },
                              step = function() {
                                
                                private$cur_iter <- private$cur_iter + 1L
                                f_val <- self$f()
                                df_val <- self$df()
                                
                                old_a <- private$a
                                
                                if (private$search == "armijo") {
                                  deltaG <- private$S - old_a
                                  derphi0 <- sum(deltaG * df_val)
                                  f <- function(x, dx, alpha, ...) {
                                    
                                    # xnew <- simplex_proj(x)
                                    # proj <- x * exp( dx * alpha)
                                    # private$a <- proj / sum(proj)
                                    private$a <- x + dx * alpha
                                    self$solve_param()
                                    
                                    
                                    loss <- self$f() #self$f()
                                    if (loss < 0) return(f_val)
                                    return(loss)
                                  }
                                  
                                  if (-derphi0 >  private$tol) {
                                    search_res <-  scalar_search_armijo(phi = f,
                                                                        phi0 = f_val,
                                                                        derphi0 = derphi0,
                                                                        x = old_a,
                                                                        dx = c(deltaG),
                                                                        c1 = 1e-4, alpha0 = 0.99, 
                                                                        amin = 0)
                                    if ( !is.null(search_res$alpha)) {
                                      search_res$alpha = min(1, max(search_res$alpha,0))
                                      private$a <- old_a + search_res$alpha * deltaG
                                    }
                                    
                                  } else {
                                    private$a <- old_a
                                  }
                                } else if (private$search == "fixed") {
                                  deltaG <- private$S - old_a
                                  derphi0 <- sum(deltaG * df_val)
                                  if (-derphi0 >  private$tol) {
                                    alpha <- min(max(private$stepsize * -derphi0 / sum(deltaG^2), 0), 1)
                                    private$a <- old_a + alpha * deltaG
                                  }
                                } else if (private$search == "mirror") {
                                  
                                  deltaG <- exp(private$S) - exp(old_a)
                                  derphi0 <- sum(deltaG * df_val)
                                  if (-derphi0 >  private$tol) {
                                    # eta <- private$stepsize / sqrt(i)
                                    alpha <- min(private$stepsize * -derphi0 / sum(deltaG^2))
                                    prop <- log(exp(old_a) + deltaG * alpha)
                                    private$a <- prop/sum(prop)
                                  }
                                } else if (private$search == "mirror-nocgd") {
                                  # epsilon <- private$stepsize / sqrt(private$cur_iter)
                                  # epsilon <- private$stepsize * sqrt(2 * private$maxcost / private$cur_iter)
                                  # if (f_val > private$tol) {
                                  #   M_k = sum(df_val^2)
                                  #   h_k <- epsilon/M_k
                                  #   prop <- old_a * exp(-h_k * df_val)
                                  # } else {
                                  #   dg  = -df_val
                                  #   M_k = sum(dg^2)
                                  #   h_k <- epsilon/M_k
                                  #   prop <- old_a * exp(-h_k * dg)
                                  # }
                                  eta <- private$stepsize / sqrt(private$cur_iter)
                                  prop <- old_a * exp(-eta * df_val)
                                  private$a <- prop/sum(prop)
                                } else if (private$search == "mirror-accelerated") {
                                  
                                  step_ <- (private$cur_iter + 1)/2
                                  
                                  
                                  private$a_tilde <- private$a_tilde * exp(- private$stepsize * step_ * df_val)
                                  private$a_tilde <- private$a_tilde/sum(private$a_tilde)
                                  
                                  private$a_hat <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                  
                                  step_ <- step_ + 0.5
                                  
                                  private$a <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                  # update parameters next
                                  
                                }
                                
                                
                              },
                              initialize = function(X1, X2, 
                                                    cost,
                                                    qp_solver = c("mosek", "gurobi", 
                                                                  "cplex"),
                                                    lambda = 100,
                                                    add.margins = FALSE,
                                                    metric = dist.metrics(),
                                                    power = 2,
                                                    niter = 1000,
                                                    tol = 1e-7,
                                                    search = "mirror-accelerated",
                                                    stepsize = 1e-1,
                                                    sample_weight = NULL,
                                                    reach = NULL,
                                                    diameter = NULL,
                                                    scaling = 0.5, truncate = 5,
                                                    kernel = NULL,
                                                    cluster_scale = NULL, 
                                                    debias = TRUE, 
                                                    verbose = FALSE, backend='auto',
                                                    ...
                              ) {
                                metric    <- match.arg(metric, dist.metrics())
                                # private$add.margins <- isTRUE(add.margins)
                                # private$add.mapping <- isTRUE(add.mapping)
                                private$penalty <- "entropy"
                                if (missing(niter) || length(niter) == 0) niter <- 1000
                                if (missing(tol) || length(tol) == 0) tol <= 1e-7 
                                private$search <- match.arg(search, c("armijo","mirror",
                                                                      "fixed", "mirror-nocgd",
                                                                      "mirror-accelerated"
                                ))
                                
                                private$cur_iter <- 0
                                
                                private$n1 <- nrow(X1)
                                private$n2 <- nrow(X2)
                                
                                x <- rbind(X1,X2)
                                z <- c(rep(0, private$n1), rep(1, private$n2))
                                
                                sw  <- get_sample_weight(sample_weight, z)
                                private$a <- private$a_hat <- private$a_tilde <- sw$a
                                private$b <- sw$b
                                
                                
                                if (missing(cost) || length(cost) == 0) {
                                  private$cost <- cost_fun(x = x, z = z, power = power, metric = metric,
                                                           estimand = "ATT")^power
                                  # if (add.margins) {
                                  #   costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
                                  #                                                       estimand = "ATT")),
                                  #               list(private$cost^(1/power)))
                                  # } else {
                                  # costqp <- private$cost^(1/power)
                                  # }
                                } else {
                                  # if (add.margins) {
                                  #   costqp <- cost
                                  #   private$cost   <- cost[[length(cost)]]^(power)
                                  # } else {
                                  private$cost <- cost^(power)
                                  # }
                                }
                                private$maxcost <- max(private$cost)
                                if (is.null(diameter)) diameter <- private$maxcost
                                private$stepsize <- stepsize * sqrt(2 * private$maxcost)
                                # if (add.mapping) private$cost <- private$cost / max(private$cost)
                                
                                # private$penalty_list <- list(margins = qp_constraint$margins)
                                
                                # private$qp <- qp_lin_comb(private$a)
                                
                                # stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
                                # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
                                # private$cost_idx <- grep("cost", names(private$qp$obj$L))
                                private$solver <- switch(qp_solver,
                                                         "cplex" = cplex_solver,
                                                         "gurobi" = gurobi_solver,
                                                         "mosek" = mosek_solver)
                                private$tol <- tol
                                private$niter <- as.numeric(niter)
                                
                                if (metric == "mahalanobis") {
                                  U <- inv_sqrt_mat(cov(x), symmetric = TRUE)
                                  
                                  update <- (x - matrix(colMeans(x), nrow = nrow(x),
                                                        ncol = ncol(x), byrow = TRUE)) %*% U
                                  
                                  X1 <- update[1:nrow(X1),,drop = FALSE]
                                  X2 <- update[-(1:nrow(X1)),,drop = FALSE]
                                  
                                } else if (metric == "sdLp") {
                                  update <- scale(x)
                                  X1 <- update[1:nrow(X1),,drop = FALSE]
                                  X2 <- update[-(1:nrow(X1)),,drop = FALSE]
                                }
                                
                                private$metric <- metric
                                
                                private$p <- power
                                private$X1 <- X1
                                private$X2 <- X2
                                
                                private$d  <- ncol(X1)
                                private$n1 <- nrow(X1)
                                private$n2 <- nrow(X2)
                                
                                private$lambda <- lambda
                                
                                private$sinkhorn_args <- list(blur = private$lambda, reach = reach, 
                                                              diameter = diameter,
                                                              scaling = scaling, 
                                                              truncate = truncate,
                                                              metric = "Lp", kernel = kernel,
                                                              cluster_scale = cluster_scale, 
                                                              debias = debias, 
                                                              verbose = verbose, 
                                                              backend = backend)
                                
                                
                                
                              }
                            ),
                            private = list(
                              "a" = "numeric",
                              "a_hat" = "numeric",
                              "a_old" = "numeric",
                              "a_tilde" = "numeric",
                              # "add.mapping" = "logical",
                              # "add.margins" = "logical",
                              "b" = "numeric",
                              "cost" = "numeric",
                              "cost_idx" = "numeric",
                              "cur_iter" = "integer",
                              "d" = "numeric",
                              "f_pot" = "numeric",
                              "f_val" = "numeric",
                              "f_val_old" = "numeric",
                              "G" = "numeric",
                              "G_old" = "numeric",
                              "g_pot" = "numeric",
                              "lambda" = "numeric",
                              "maxcost" = "numeric",
                              "metric" = "character",
                              "n1" = "numeric",
                              "n2" = "numeric",
                              "niter" = "numeric",
                              "p" = "numeric",
                              "penalty" = "character",
                              "S" = "numeric",
                              "search" = "character",
                              "stepsize" = "numeric",
                              "sinkhorn_args" = "list",
                              "solver" = "function",
                              "tol" = "numeric",
                              "X1" = "matrix",
                              "X2" = "matrix"
                            )
                            
                            
  )
}

#entropy based divergence
{
  wassDivEnt <- R6::R6Class("wassDivEnt", 
                           inherit = wassDiv,
                             public = list(
                               converged = function() {
                                 private$f_val <- sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
                                 # print(private$f_val)
                                 f_val_diff <- abs(private$f_val - private$f_val_old)
                                 # if ((i %% 100) == 0) message(round(sum(private$a * y[z==0]), digits = 3), ", ", appendLF = FALSE)
                                 # message(round(private$a * y[z==0]))
                                 # message(private$f_val, ", ", appendLF = FALSE)
                                 # diff_param <- (private$param - private$param_old)
                                 # paramnorm <-  sum(diff_param^2)/sum(private$param_old^2)
                                 # if(is.nan(paramnorm) & sum(diff_param^2) == 0 ) paramnorm <- 0
                                 # diff_G <- (private$G - private$G_old)
                                 # Gnorm <- sum(diff_G^2)/sum(private$G_old^2)
                                 # if(is.nan(Gnorm) & sum(diff_G^2) == 0 ) Gnorm <- 0
                                 nan.check <- isTRUE(any(is.nan(private$a)) || is.nan(private$f_val))
                                 
                                 converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol)  ||
                                   isTRUE(sum(abs(private$a_old - private$a)) < private$tol) ||
                                   nan.check
                                 # message(sum(abs(private$a_old - private$a)), ", ", appendLF = FALSE)
                                 # ||
                                 #   isTRUE(f_val_diff < private$tol) 
                                 # || (isTRUE( sqrt(paramnorm) < private$tol) &&
                                 #   isTRUE( sqrt(Gnorm) < private$tol) )
                                 # private$G_old <- private$G
                                 # private$param_old <- private$param
                                 if (nan.check) warning("NaN found in parameters! Try reducing the step size.")
                                 private$f_val_old <- private$f_val
                                 private$a_old <- private$a
                                 return(converged)
                               },
                               # f = function() {
                               #   return(
                               #     sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
                               #   )
                               # },
                               # df = function() {
                               #   return(private$f_pot)
                               #   
                               #   sol1 <- sinkhorn_geom(x = private$X1, y = private$X2, 
                               #                        a = private$a, 
                               #                        b = private$b, power = private$p, 
                               #                        blur = private$lambda, reach = private$sinkhorn_args$reach, 
                               #                        diameter = private$sinkhorn_args$diameter,
                               #                        scaling = private$sinkhorn_args$scaling, 
                               #                        truncate = private$sinkhorn_args$truncate,
                               #                        metric = "Lp", kernel = private$sinkhorn_args$kernel,
                               #                        cluster_scale=private$sinkhorn_args$cluster_scale, 
                               #                        debias=FALSE, 
                               #                        verbose=private$sinkhorn_args$verbose, 
                               #                        backend=private$sinkhorn_args$backend)
                               #   sol2 <- sinkhorn_geom(x = private$X1, y = private$X1, 
                               #                         a = private$a, 
                               #                         b = private$a, power = private$p, 
                               #                         blur = private$lambda, reach = private$sinkhorn_args$reach, 
                               #                         diameter = private$sinkhorn_args$diameter,
                               #                         scaling = private$sinkhorn_args$scaling, 
                               #                         truncate = private$sinkhorn_args$truncate,
                               #                         metric = "Lp", kernel = private$sinkhorn_args$kernel,
                               #                         cluster_scale=private$sinkhorn_args$cluster_scale, 
                               #                         debias=FALSE, 
                               #                         verbose=private$sinkhorn_args$verbose, 
                               #                         backend=private$sinkhorn_args$backend)
                               #   fmat <- matrix(sol1$f, private$n1, private$n2)
                               #   gmat <- matrix(sol1$g, private$n1, private$n2, byrow = TRUE)
                               #   eta <- (fmat + gmat - private$cost)/private$lambda
                               #   pi_raw <- exp(eta) 
                               #   
                               #   cost2 <- cost_calc_lp(private$X1,private$X1,private$p)^private$p
                               #   fmat <- matrix(sol2$f, private$n1, private$n1)
                               #   gmat <- matrix(sol2$g, private$n1, private$n1, byrow = TRUE)
                               #   eta <- (fmat + gmat - cost2)/private$lambda
                               #   pi_raw2 <- exp(eta) 
                               #   
                               #   sol1$f - sol2$f - private$lambda * c((pi_raw - 1) %*% private$b) +
                               #     private$lambda * c((pi_raw2 - 1) %*% private$a)
                               #   
                               #   
                               # },
                               # get_G = function() {
                               #   return(private$a)
                               # },
                               # get_param = function() {
                               #   return(list(f = private$f_pot,
                               #               g = private$g_pot))
                               # },
                               get_niter = function() {
                                 return(private$niter)
                               },
                               get_weight = function() {
                                 
                                 fmat <- matrix(private$f_pot, private$n1, private$n2)
                                 gmat <- matrix(private$g_pot, private$n1, private$n2, byrow = TRUE)
                                 eta <- (fmat + gmat - private$cost)/private$lambda
                                 pi_raw <- exp(eta) * matrix(private$a, private$n1, private$n2) * matrix(private$b, private$n1, private$n2, byrow = TRUE)
                                 pi <- round_pi(pi_raw, private$a, private$b)
                                 return( pi )
                               },
                               return_cw = function() {
                                 out <- list(w0 = private$a,
                                             w1 = private$b,
                                             gamma = self$get_weight(),
                                             estimand = "ATT",
                                             method = "Wasserstein",
                                             args = list(
                                               constraint = list(joint = NULL,
                                                                 penalty = private$lambda,
                                                                 margins = NULL),
                                               penalty = "entropy",
                                               power = private$power,
                                               metric = private$metric,
                                               niter = private$niter,
                                               cur_iter = private$cur_iter,
                                               add.margins = private$add.margins,
                                               joint.mapping = FALSE,
                                               add.divergence = TRUE,
                                               conditional.gradient = TRUE))
                                 class(out) <- "causalWeights"
                                 return(out)
                               },
                               solve_G = function() {
                                 self$solve_param()
                                 private$f_val_old <- self$f()
                                 private$a_old <- private$a
                               },
                               solve_param = function() {
                                 
                                 sol <- sinkhorn_geom(x = private$X1, y = private$X2, 
                                                      a = private$a, 
                                                      b = private$b, power = private$p, 
                                                      blur = private$lambda, reach = private$sinkhorn_args$reach, 
                                                      diameter = private$sinkhorn_args$diameter,
                                                      scaling = private$sinkhorn_args$scaling, 
                                                      truncate = private$sinkhorn_args$truncate,
                                                      metric = "Lp", kernel = private$sinkhorn_args$kernel,
                                                      cluster_scale=private$sinkhorn_args$cluster_scale, 
                                                      debias=private$sinkhorn_args$debias, 
                                                      verbose=private$sinkhorn_args$verbose, 
                                                      backend=private$sinkhorn_args$backend)
                                 private$f_pot <- sol$f
                                 private$g_pot <- sol$g
                                 
                               },
                               solve_S = function() {
                                 if (private$search != "mirror-nocgd" && 
                                     private$search != "mirror-accelerated") {
                                   private$qp <- qp_dual(f = private$f_pot,
                                                         g = private$g_pot,
                                                         b = private$b)
                                   private$S <- private$solver(private$qp)$sol
                                   
                                   # private$S <- as.numeric(private$f_pot == min(private$f_pot))
                                   
                                   # grad <- self$df()
                                   # private$S <- renormalize(as.numeric(grad == min(grad)))
                                 }
                                 
                               },
                               step = function() {
                                 
                                 private$cur_iter <- private$cur_iter + 1L
                                 f_val <- self$f()
                                 df_val <- self$df()
                                 
                                 old_a <- private$a
                                 
                                 if (private$search == "armijo") {
                                   deltaG <- private$S - old_a
                                   derphi0 <- sum(deltaG * df_val)
                                   f <- function(x, dx, alpha, ...) {
                                     
                                     # xnew <- simplex_proj(x)
                                     # proj <- x * exp( dx * alpha)
                                     # private$a <- proj / sum(proj)
                                     private$a <- x + dx * alpha
                                     self$solve_param()
                                     
                                     
                                     loss <- self$f() #self$f()
                                     if (loss < 0) return(f_val)
                                     return(loss)
                                   }
                                   
                                   if (-derphi0 >  private$tol) {
                                     search_res <-  scalar_search_armijo(phi = f,
                                                                         phi0 = f_val,
                                                                         derphi0 = derphi0,
                                                                         x = old_a,
                                                                         dx = c(deltaG),
                                                                         c1 = 1e-4, alpha0 = 0.99, 
                                                                         amin = 0)
                                     if ( !is.null(search_res$alpha)) {
                                       search_res$alpha = min(1, max(search_res$alpha,0))
                                       private$a <- old_a + search_res$alpha * deltaG
                                     }
                                     
                                   } else {
                                     private$a <- old_a
                                   }
                                 } else if (private$search == "fixed") {
                                   deltaG <- private$S - old_a
                                   derphi0 <- sum(deltaG * df_val)
                                   if (-derphi0 >  private$tol) {
                                     alpha <- min(max(private$stepsize * -derphi0 / sum(deltaG^2), 0), 1)
                                     private$a <- old_a + alpha * deltaG
                                   }
                                 } else if (private$search == "mirror") {
                                   
                                   deltaG <- exp(private$S) - exp(old_a)
                                   derphi0 <- sum(deltaG * df_val)
                                   if (-derphi0 >  private$tol) {
                                     # eta <- private$stepsize / sqrt(i)
                                     alpha <- min(private$stepsize * -derphi0 / sum(deltaG^2))
                                     prop <- log(exp(old_a) + deltaG * alpha)
                                     private$a <- prop/sum(prop)
                                   }
                                 } else if (private$search == "mirror-nocgd") {
                                   # epsilon <- private$stepsize / sqrt(private$cur_iter)
                                   # epsilon <- private$stepsize * sqrt(2 * private$maxcost / private$cur_iter)
                                   # if (f_val > private$tol) {
                                   #   M_k = sum(df_val^2)
                                   #   h_k <- epsilon/M_k
                                   #   prop <- old_a * exp(-h_k * df_val)
                                   # } else {
                                   #   dg  = -df_val
                                   #   M_k = sum(dg^2)
                                   #   h_k <- epsilon/M_k
                                   #   prop <- old_a * exp(-h_k * dg)
                                   # }
                                   eta <- private$stepsize / sqrt(private$cur_iter)
                                   prop <- old_a * exp(-eta * df_val)
                                   private$a <- prop/sum(prop)
                                 } else if (private$search == "mirror-accelerated") {
                                   
                                   step_ <- (private$cur_iter + 1)/2
                                   
                                   
                                   private$a_tilde <- private$a_tilde * exp(- private$stepsize * step_ * df_val)
                                   private$a_tilde <- private$a_tilde/sum(private$a_tilde)
                                   
                                   private$a_hat <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                   
                                   step_ <- step_ + 0.5
                                   
                                   private$a <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                   # update parameters next
                                   
                                 }
                                 

                               },
                             initialize = function(X1, X2, 
                                                   cost,
                                                   qp_solver = c("mosek", "gurobi", 
                                                                 "cplex"),
                                                   lambda = 100,
                                                   add.margins = FALSE,
                                                   metric = dist.metrics(),
                                                   power = 2,
                                                   niter = 1000,
                                                   tol = 1e-7,
                                                   search = "mirror-accelerated",
                                                   stepsize = 1e-1,
                                                   sample_weight = NULL,
                                                   reach = NULL,
                                                   diameter = NULL,
                                                   scaling = 0.5, truncate = 5,
                                                   kernel = NULL,
                                                   cluster_scale = NULL, 
                                                   debias = TRUE, 
                                                   verbose = FALSE, backend='auto',
                                                   ...
                             ) {
                               metric    <- match.arg(metric, dist.metrics())
                               # private$add.margins <- isTRUE(add.margins)
                               # private$add.mapping <- isTRUE(add.mapping)
                               private$penalty <- "entropy"
                               if (missing(niter) || length(niter) == 0) niter <- 1000
                               if (missing(tol) || length(tol) == 0) tol <= 1e-7 
                               private$search <- match.arg(search, c("armijo","mirror",
                                                                     "fixed", "mirror-nocgd",
                                                                     "mirror-accelerated"
                                                                     ))
                               
                               private$cur_iter <- 0
                               
                               private$n1 <- nrow(X1)
                               private$n2 <- nrow(X2)
                               
                               x <- rbind(X1,X2)
                               z <- c(rep(0, private$n1), rep(1, private$n2))
                               
                               sw  <- get_sample_weight(sample_weight, z)
                               private$a <- private$a_hat <- private$a_tilde <- sw$a
                               private$b <- sw$b
                               
                               
                               if (missing(cost) || length(cost) == 0) {
                                 private$cost <- cost_fun(x = x, z = z, power = power, metric = metric,
                                                          estimand = "ATT")^power
                                 # if (add.margins) {
                                 #   costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
                                 #                                                       estimand = "ATT")),
                                 #               list(private$cost^(1/power)))
                                 # } else {
                                 # costqp <- private$cost^(1/power)
                                 # }
                               } else {
                                 # if (add.margins) {
                                 #   costqp <- cost
                                 #   private$cost   <- cost[[length(cost)]]^(power)
                                 # } else {
                                 private$cost <- cost^(power)
                                 # }
                               }
                               private$maxcost <- max(private$cost)
                               if (is.null(diameter)) diameter <- private$maxcost
                               private$stepsize <- stepsize * sqrt(2 * private$maxcost)
                               # if (add.mapping) private$cost <- private$cost / max(private$cost)
                               
                               # private$penalty_list <- list(margins = qp_constraint$margins)
                               
                               # private$qp <- qp_lin_comb(private$a)
                               
                               # stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
                               # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
                               # private$cost_idx <- grep("cost", names(private$qp$obj$L))
                               private$solver <- switch(qp_solver,
                                                "cplex" = cplex_solver,
                                                "gurobi" = gurobi_solver,
                                                "mosek" = mosek_solver)
                               private$tol <- tol
                               private$niter <- as.numeric(niter)
                               
                               if (metric == "mahalanobis") {
                                 U <- inv_sqrt_mat(cov(x), symmetric = TRUE)
                                 
                                 update <- (x - matrix(colMeans(x), nrow = nrow(x),
                                                           ncol = ncol(x), byrow = TRUE)) %*% U
                                 
                                 X1 <- update[1:nrow(X1),,drop = FALSE]
                                 X2 <- update[-(1:nrow(X1)),,drop = FALSE]
                                 
                               } else if (metric == "sdLp") {
                                 update <- scale(x)
                                 X1 <- update[1:nrow(X1),,drop = FALSE]
                                 X2 <- update[-(1:nrow(X1)),,drop = FALSE]
                               }
                               
                               private$metric <- metric
                               
                               private$p <- power
                               private$X1 <- X1
                               private$X2 <- X2
                               
                               private$d  <- ncol(X1)
                               private$n1 <- nrow(X1)
                               private$n2 <- nrow(X2)
                               
                               private$lambda <- lambda
                               
                               private$sinkhorn_args <- list(blur = private$lambda, reach = reach, 
                                                             diameter = diameter,
                                                             scaling = scaling, 
                                                             truncate = truncate,
                                                             metric = "Lp", kernel = kernel,
                                                             cluster_scale = cluster_scale, 
                                                             debias = debias, 
                                                             verbose = verbose, 
                                                             backend = backend)
                               
                               
                               
                             }
                           ),
                           private = list(
                             "a" = "numeric",
                             "a_hat" = "numeric",
                             "a_old" = "numeric",
                             "a_tilde" = "numeric",
                             # "add.mapping" = "logical",
                             # "add.margins" = "logical",
                             "b" = "numeric",
                             "cost" = "numeric",
                             "cost_idx" = "numeric",
                             "cur_iter" = "integer",
                             "d" = "numeric",
                             "f_pot" = "numeric",
                             "f_val" = "numeric",
                             "f_val_old" = "numeric",
                             "G" = "numeric",
                             "G_old" = "numeric",
                             "g_pot" = "numeric",
                             "lambda" = "numeric",
                             "maxcost" = "numeric",
                              "metric" = "character",
                             "n1" = "numeric",
                             "n2" = "numeric",
                             "niter" = "numeric",
                             "p" = "numeric",
                             "penalty" = "character",
                             "S" = "numeric",
                             "search" = "character",
                             "stepsize" = "numeric",
                             "sinkhorn_args" = "list",
                             "solver" = "function",
                             "tol" = "numeric",
                             "X1" = "matrix",
                             "X2" = "matrix"
                           )
                           
                           
  )
}

#L2 based divergence
{
  wassDivL2 <- R6::R6Class("wassDivL2", 
                            inherit = cgOptimizer,
                            public = list(
                              converged = function() {
                                private$f_val <- self$f()
                                f_val_diff <- abs(private$f_val - private$f_val_old)
                                message(private$f_val, ", ", appendLF = FALSE)
                                # diff_param <- (private$param - private$param_old)
                                # paramnorm <-  sum(diff_param^2)/sum(private$param_old^2)
                                # if(is.nan(paramnorm) & sum(diff_param^2) == 0 ) paramnorm <- 0
                                # diff_G <- (private$G - private$G_old)
                                # Gnorm <- sum(diff_G^2)/sum(private$G_old^2)
                                # if(is.nan(Gnorm) & sum(diff_G^2) == 0 ) Gnorm <- 0
                                
                                converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol)  ||
                                  isTRUE(sum(abs(private$a_old - private$a)) < private$tol)
                                # message(sum(abs(private$a_old - private$a)), ", ", appendLF = FALSE)
                                # ||
                                #   isTRUE(f_val_diff < private$tol) 
                                # || (isTRUE( sqrt(paramnorm) < private$tol) &&
                                #   isTRUE( sqrt(Gnorm) < private$tol) )
                                # private$G_old <- private$G
                                # private$param_old <- private$param
                                private$f_val_old <- private$f_val
                                private$a_old <- private$a
                                return(converged)
                              },
                              f = function() {
                                return(
                                  sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
                                )
                              },
                              df = function() {
                                return(private$f_pot)
                              },
                              get_G = function() {
                                return(private$a)
                              },
                              get_param = function() {
                                return(list(f = private$f_pot,
                                            g = private$g_pot))
                              },
                              get_niter = function() {
                                return(private$niter)
                              },
                              get_weight = function() {
                                
                                fmat <- matrix(private$f_pot, private$n1, private$n2)
                                gmat <- matrix(private$g_pot, private$n1, private$n2, byrow = TRUE)
                                eta <- (fmat + gmat - private$cost)
                                return((eta * (eta > 0))/private$lambda)
                              },
                              return_cw = function() {
                                out <- list(w0 = private$a,
                                            w1 = private$b,
                                            gamma = self$get_weight(),
                                            estimand = "ATT",
                                            method = "Wasserstein",
                                            args = list(
                                              constraint = list(joint = NULL,
                                                                penalty = private$lambda,
                                                                margins = NULL),
                                              penalty = "L2",
                                              power = private$power,
                                              metric = private$metric,
                                              niter = private$niter,
                                              cur_iter = private$cur_iter,
                                              add.margins = private$add.margins,
                                              joint.mapping = FALSE,
                                              add.divergence = TRUE,
                                              conditional.gradient = TRUE))
                                class(out) <- "causalWeights"
                                return(out)
                              },
                              solve_G = function() {
                                self$solve_param()
                                private$f_val_old <- self$f()
                                private$a_old <- private$a
                              },
                              solve_param = function() {
                                
                                fit_xy <- lbfgsb3c::lbfgsb3(par = private$optimizer_xy$init(),
                                                            fn = private$optimizer_xy$obj,
                                                            gr = private$optimizer_xy$grad,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control = private$control
                                )
                                
                                fit_xx <- lbfgsb3c::lbfgsb3(par = private$optimizer_xx$init(),
                                                            fn = private$optimizer_xx$obj,
                                                            gr = private$optimizer_xx$grad,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control = private$control
                                )
                                
                                if (length(private$q_pot) == 0) {
                                  fit_yy <- lbfgsb3c::lbfgsb3(par = private$optimizer_yy$init(),
                                                              fn = private$optimizer_yy$obj,
                                                              gr = private$optimizer_yy$grad,
                                                              lower = -Inf,
                                                              upper = Inf,
                                                              control = private$control
                                  )
                                  private$q_pot <- c(private$optimizer_yy$get_f(fit_yy$par))
                                }
                                
                                private$f_pot <- c(private$optimizer_xy$get_f(fit_xy$par) - private$optimizer_xx$get_f(fit_xx$par))
                                private$g_pot <- c(private$optimizer_xy$get_g(fit_xy$par) - private$q_pot)
                              },
                              solve_S = function() {
                                private$qp <- qp_dual(f = private$f_pot, 
                                                      g = private$g_pot, 
                                                      b = private$b)
                                private$S <- private$solver(private$qp)$sol
                              },
                              initialize = function(X1, X2, 
                                                    cost = NULL,
                                                    qp_solver = c("mosek", "gurobi", 
                                                                  "cplex"),
                                                    lambda = 100,
                                                    add.margins = FALSE,
                                                    metric = dist.metrics(),
                                                    power = 2,
                                                    niter = 1000,
                                                    tol = 1e-7,
                                                    sample_weight = NULL,
                                                    ...
                              ) {
                                metric    <- match.arg(metric, dist.metrics())
                                # private$add.margins <- isTRUE(add.margins)
                                # private$add.mapping <- isTRUE(add.mapping)
                                private$penalty <- "L2"
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
                                  cost_joint <- cost_fun(x = x, z = z, power = power, metric = metric,
                                                           estimand = "ATT")
                                  cost_a <- cost_fun(x = rbind(X1,X1), z = c(rep(1, nrow(X1)),
                                                                             rep(0, nrow(X1))),
                                                     power = power, metric = metric,
                                                      estimand = "ATT")
                                  cost_b <- cost_fun(x = rbind(X2,X2), z = c(rep(1, nrow(X2)),
                                                                             rep(0, nrow(X2))),
                                                     power = power, metric = metric,
                                                     estimand = "ATT")
                                  private$cost <- cost_joint^power
                                  # if (add.margins) {
                                  #   costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
                                  #                                                       estimand = "ATT")),
                                  #               list(private$cost^(1/power)))
                                  # } else {
                                  # costqp <- private$cost^(1/power)
                                  # }
                                } else {
                                  # if (add.margins) {
                                  #   costqp <- cost
                                  #   private$cost   <- cost[[length(cost)]]^(power)
                                  # } else {
                                  private$cost <- cost$joint^(power)
                                  cost_joint <- cost$joint
                                  cost_a <- cost$a
                                  cost_b <- cost$b
                                  # }
                                }
                                private$q_pot <- numeric(0)
                                # if (add.mapping) private$cost <- private$cost / max(private$cost)
                                
                                # private$penalty_list <- list(margins = qp_constraint$margins)
                                
                                # private$qp <- qp_lin_comb(private$a)
                                
                                # stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
                                # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
                                # private$cost_idx <- grep("cost", names(private$qp$obj$L))
                                private$solver <- switch(qp_solver,
                                                         "cplex" = cplex_solver,
                                                         "gurobi" = gurobi_solver,
                                                         "mosek" = mosek_solver)
                                private$tol <- tol
                                private$niter <- as.numeric(niter)
                                private$cur_iter <- 0
                                
                                private$metric <- metric
                                
                                private$p <- power
                                private$X1 <- X1
                                private$X2 <- X2
                                
                                private$d  <- ncol(X1)
                                private$n1 <- nrow(X1)
                                private$n2 <- nrow(X2)
                                
                                private$lambda <- lambda
                                
                                private$optimizer_xy <- otDualL2$new(lambda = private$lambda,
                                                                     cost = cost_joint,
                                                                     p = private$p,
                                                                     a = private$a,
                                                                     b = private$b)
                                
                                private$optimizer_xx <- otDualL2_self$new(lambda = private$lambda,
                                                                          cost = cost_a,
                                                                          p = private$p,
                                                                          a = private$a,
                                                                          b = private$a)
                                
                                private$optimizer_yy <- otDualL2_self$new(lambda = private$lambda,
                                                                          cost = cost_b,
                                                                          p = private$p,
                                                                          a = private$b,
                                                                          b = private$b)
                                
                                
                              }
                            ),
                            private = list(
                              "a" = "numeric",
                              "a_hat" = "numeric",
                              "a_old" = "numeric",
                              "a_tilde" = "numeric",
                              "b" = "numeric",
                              "cost" = "numeric",
                              "cost_idx" = "numeric",
                              "cur_iter" = "integer",
                              "d" = "numeric",
                              "f_pot" = "numeric",
                              "f_val" = "numeric",
                              "f_val_old" = "numeric",
                              "G" = "numeric",
                              "G_old" = "numeric",
                              "g_pot" = "numeric",
                              "lambda" = "numeric",
                              "metric" = "character",
                              "n1" = "numeric",
                              "n2" = "numeric",
                              "niter" = "numeric",
                              "optimizer_xx" = "R6",
                              "optimizer_xy" = "R6",
                              "optimizer_yy" = "R6",
                              "p" = "numeric",
                              "q_pot" = "numeric",
                              "penalty" = "character",
                              "S" = "numeric",
                              "sinkhorn_args" = "list",
                              "solver" = "function",
                              "tol" = "numeric",
                              "X1" = "matrix",
                              "X2" = "matrix"
                            )
                            
                            
  )
}
# {
#   wassDivL2 <- R6::R6Class("wassDivL2", 
#                            inherit = cgOptimizer,
#                            public = list(
#                              public = list(
#                                converged = function() {
#                                  private$f_val <- self$f()
#                                  f_val_diff <- abs(private$f_val - private$f_val_old)
#                                  # diff_param <- (private$param - private$param_old)
#                                  # paramnorm <-  sum(diff_param^2)/sum(private$param_old^2)
#                                  # if(is.nan(paramnorm) & sum(diff_param^2) == 0 ) paramnorm <- 0
#                                  # diff_G <- (private$G - private$G_old)
#                                  # Gnorm <- sum(diff_G^2)/sum(private$G_old^2)
#                                  # if(is.nan(Gnorm) & sum(diff_G^2) == 0 ) Gnorm <- 0
#                                  
#                                  converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol) 
#                                  # ||
#                                  #   isTRUE(f_val_diff < private$tol) 
#                                  # || (isTRUE( sqrt(paramnorm) < private$tol) &&
#                                  #   isTRUE( sqrt(Gnorm) < private$tol) )
#                                  # private$G_old <- private$G
#                                  # private$param_old <- private$param
#                                  private$f_val_old <- private$f_val
#                                  return(converged)
#                                },
#                                f = function() {
#                                  return(
#                                    sum(private$f_pot * private$a) + sum(private$g_pot * private$b)
#                                  )
#                                },
#                                df = function() {
#                                  return(private$f_pot)
#                                },
#                                get_G = function() {
#                                  return(private$G)
#                                },
#                                get_param = function() {
#                                  return(list(f = private$f_pot,
#                                              g = private$g_pot))
#                                },
#                                get_niter = function() {
#                                  return(private$niter)
#                                },
#                                
#                                return_cw = function() {
#                                  out <- list(w0 = private$a,
#                                              w1 = private$b,
#                                              gamma = self$get_weight(),
#                                              estimand = "ATT",
#                                              method = "Wasserstein",
#                                              args = list(
#                                                constraint = list(joint = NULL,
#                                                                  penalty = private$lambda,
#                                                                  margins = NULL),
#                                                penalty = "L2",
#                                                power = private$power,
#                                                metric = private$metric,
#                                                niter = private$niter,
#                                                cur_iter = private$cur_iter,
#                                                add.margins = private$add.margins,
#                                                joint.mapping = FALSE,
#                                                add.divergence = TRUE,
#                                                conditional.gradient = TRUE))
#                                  class(out) <- "causalWeights"
#                                  return(out)
#                                },
#                                solve_G = function() {
#                                  private$solve_param()
#                                  private$f_val_old <- self$f()
#                                  
#                                },
#                                solve_param = function() {
#                                  
#                                  fit_xy <- lbfgsb3c::lbfgsb3(par = private$optimizer_xy$init(),
#                                                           fn = private$optimizer_xy$obj,
#                                                           gr = private$optimizer_xy$grad,
#                                                           lower = -Inf,
#                                                           upper = Inf,
#                                                           control = private$control
#                                  )
#                                  
#                                  fit_xx <- lbfgsb3c::lbfgsb3(par = private$optimizer_xx$init(),
#                                                              fn = private$optimizer_xx$obj,
#                                                              gr = private$optimizer_xx$grad,
#                                                              lower = -Inf,
#                                                              upper = Inf,
#                                                              control = private$control
#                                  )
#                                  
#                                  if (is.null(private$q_pot)) {
#                                    fit_yy <- lbfgsb3c::lbfgsb3(par = private$optimizer_yy$init(),
#                                                                fn = private$optimizer_yy$obj,
#                                                                gr = private$optimizer_yy$grad,
#                                                                lower = -Inf,
#                                                                upper = Inf,
#                                                                control = private$control
#                                    )
#                                    private$q_pot <- c(private$optimizer_yy$get_f(fit_yy$par))
#                                  }
#                                  
#                                  private$f_pot <- c(private$optimizer_xy$get_f(fit_xy$par) - private$optimizer_xx$get_f(fit_xx$par))
#                                  private$g_pot <- c(private$optimizer_xy$get_g(fit_xy$par) - private$q_pot)
#                                  
#                                },
#                                solve_S = function() {
#                                  private$qp <- qp_dual(f = private$f_pot, 
#                                                        g = private$g_pot, 
#                                                        b = private$b)
#                                  private$S <- private$solver(private$qp)
#                                },
#                                step = function() {
#                                  private$cur_iter <- private$cur_iter + 1L
#                                  olda <- private$a
#                                  deltaG <- private$S - olda
#                                  f_val <- self$f()
#                                  df_val <- self$df()
#                                  
#                                  f <- function(x, ...) {
#                                    
#                                    xnew <- simplex_proj(x)
#                                    private$a <- xnew
#                                    private$optimizer_xy$update_a(xnew)
#                                    private$optimizer_xx$update_a(xnew)
#                                    private$solve_param()
#                                    
#                                    loss <- private$f()
#                                    if(loss < 0) return(f_val)
#                                    return(loss)
#                                  }
#                                  # search_res <- rje::armijo(
#                                  #   fun = f,
#                                  #   x = private$G,
#                                  #   dx = deltaG,
#                                  #   beta = 3,
#                                  #   sigma = 1e-10,
#                                  #   grad = df_val,
#                                  #   maximise = FALSE,
#                                  #   searchup = TRUE,
#                                  #   adj.start = 1
#                                  # )
#                                  # private$G <- private$G + search_res$adj * deltaG
#                                  derphi0 <- sum(deltaG * df_val)
#                                  if (-derphi0 >  private$tol) {
#                                    search_res <-  scalar_search_armijo(phi = f,
#                                                                        phi0 = f_val,
#                                                                        derphi0 = derphi0,
#                                                                        x = olda,
#                                                                        dx = c(deltaG),
#                                                                        c1 = 1e-4, alpha0 = .99, amin = 0)
#                                    if ( !is.null(search_res$alpha)) {
#                                      search_res$alpha = min(1, max(search_res$alpha,0))
#                                      private$G <- olda + search_res$alpha * deltaG
#                                    }
#                                    
#                                  } else {
#                                    private$G <- olda
#                                  }
#                                  private$optimizer_xy$update_a(private$G)
#                                  private$optimizer_xx$update_a(private$G)
#                                  private$a <- private$G
#                                  
#                                }
#                              ),
#                              initialize = function(X1, X2, 
#                                                    cost,
#                                                    qp_solver = c("mosek", "gurobi", 
#                                                                  "cplex"),
#                                                    lambda = 100,
#                                                    add.margins = FALSE,
#                                                    metric = dist.metrics(),
#                                                    power = 2,
#                                                    niter = 1000,
#                                                    tol = 1e-7,
#                                                    sample_weight = NULL,
#                                                    ...
#                              ) {
#                                metric    <- match.arg(metric, dist.metrics())
#                                # private$add.margins <- isTRUE(add.margins)
#                                # private$add.mapping <- isTRUE(add.mapping)
#                                private$penalty <- "L2"
#                                if (missing(niter) || length(niter) == 0) niter <- 1000
#                                if (missing(tol) || length(tol) == 0) tol <= 1e-7 
#                                
#                                
#                                private$n1 <- nrow(X1)
#                                private$n2 <- nrow(X2)
#                                
#                                x <- rbind(X1,X2)
#                                z <- c(rep(0, private$n1), rep(1, private$n2))
#                                
#                                sw  <- get_sample_weight(sample_weight, z)
#                                private$a <- sw$a
#                                private$b <- sw$b
#                                
#                                if (missing(cost) || length(cost) == 0) {
#                                  private$cost <- cost_fun(x = x, z = z, power = power, metric = metric,
#                                                           estimand = "ATT")^power
#                                  # if (add.margins) {
#                                  #   costqp <- c(lapply(1:ncol(X1), function(d) cost_fun(x = x[,d,drop = FALSE], z = z, power = power, metric = metric,
#                                  #                                                       estimand = "ATT")),
#                                  #               list(private$cost^(1/power)))
#                                  # } else {
#                                    # costqp <- private$cost^(1/power)
#                                  # }
#                                } else {
#                                  # if (add.margins) {
#                                  #   costqp <- cost
#                                  #   private$cost   <- cost[[length(cost)]]^(power)
#                                  # } else {
#                                    private$cost <- cost^(power)
#                                  # }
#                                }
#                                private$metric <- metric
#                                # if (add.mapping) private$cost <- private$cost / max(private$cost)
#                                
#                                # private$penalty_list <- list(margins = qp_constraint$margins)
#                                
#                                # private$qp <- qp_lin_comb(private$a)
#                                
#                                # stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
#                                # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
#                                # private$cost_idx <- grep("cost", names(private$qp$obj$L))
#                                solver <- switch(qp_solver,
#                                                 "cplex" = cplex_solver,
#                                                 "gurobi" = gurobi_solver,
#                                                 "mosek" = mosek_solver)
#                                private$tol <- tol
#                                private$niter <- as.numeric(niter)
#                                
#                                private$p <- power
#                                private$X1 <- X1
#                                private$X2 <- X2
#                                
#                                private$d  <- ncol(X1)
#                                private$n1 <- nrow(X1)
#                                private$n2 <- nrow(X2)
#                                
#                                private$lambda <- c(cost = lambda$joint, coefficients = 0.0,
#                                                    gamma = lambda$penalty)
#                                
#                                private$optimizer_xx <- otDualL2_self$new(lambda = private$lambda,
#                                                                          cost = private$cost,
#                                                                          p = private$p,
#                                                                          a = private$a,
#                                                                          b = private$a)
#                                private$optimizer_xy <- otDualL2$new(lambda = private$lambda,
#                                                                          cost = private$cost,
#                                                                          p = private$p,
#                                                                          a = private$a,
#                                                                          b = private$b)
#                                private$optimizer_yy <- otDualL2_self$new(lambda = private$lambda,
#                                                                          cost = (private$cost)^(1/power),
#                                                                          p = private$p,
#                                                                          a = private$b,
#                                                                          b = private$b)
#                                
#                              }
#                            ),
#                            private = list(
#                              "a" = "numeric",
#                              # "add.mapping" = "logical",
#                              # "add.margins" = "logical",
#                              "b" = "numeric",
#                              "cost" = "numeric",
#                              "cost_idx" = "numeric",
#                              "cur_iter" = "integer",
#                              "d" = "numeric",
#                              "f_pot" = "numeric",
#                              "f_val" = "numeric",
#                              "f_val_old" = "numeric",
#                              "G" = "numeric",
#                              "G_old" = "numeric",
#                              "g_pot" = "numeric",
#                              "lambda" = "numeric",
#                              "n1" = "numeric",
#                              "n2" = "numeric",
#                              "niter" = "numeric",
#                              "optimizer_xx" = "R6",
#                              "optimizer_xy" = "R6",
#                              "optimizer_yy" = "R6",
#                              "p" = "numeric",
#                              "penalty" = "character",
#                              "S" = "numeric",
#                              "solver" = "function",
#                              "tol" = "numeric",
#                              "X1" = "matrix",
#                              "X2" = "matrix"
#                            )
#                            
#                            
#   )
# }


# calculate the wass. distance using CGD
fullWassCGD <- R6::R6Class("wassCGD",
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
                                          
                                          private$qp <- qp_lin_comb(private$a)
                                          
                                          stopifnot(all(c("obj", "LC","bounds","nvar") %in% names(private$qp)))
                                          # names(private$qp$obj$L) <- c(rep("cost",length(private$cost)), rep("pen", length(private$qp$obj$L) - length(cost)))
                                          private$cost_idx <- grep("cost", names(private$qp$obj$L))
                                          solver <- switch(qp_solver,
                                                           "cplex" = cplex_solver,
                                                           "gurobi" = gurobi_solver,
                                                           "mosek" = mosek_solver)
                                          private$solver <- function(qp){
                                            sol <- solver(qp)
                                            return(matrix(sol$sol, private$n1, private$n2))
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
        return(matrix(sol$sol, private$n1, private$n2))
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
      private$cur_iter <- 0
    }
  )
                                              
                                              
)

# calculate joint weight and mapping
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
                                  private$cur_iter <- 0
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
                                    return(matrix(sol$sol, private$n1, private$n2))
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

# helper functions, no joint mapping, no outcome
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

# joint mapping, without the outcome
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

# joint mapping with linear outcome
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
  
  phi_a0 = phi(x, dx, alpha0)
  if (phi_a0 <= phi0 + c1 * alpha0 * derphi0) {
    return(list(alpha = alpha0, ph1 = phi_a0))
  }
  
  # Otherwise, compute the minimizer of a quadratic interpolant:
  alpha1 = -(derphi0) * alpha0^2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
  phi_a1 = phi(x, dx, alpha1)
  if ((phi_a1 <= (phi0 + c1 * alpha1 * derphi0) )  && (alpha1 >= 0) ) {
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
    if(factor == 0) break
    a = alpha0^2 * (phi_a1 - phi0 - derphi0 * alpha1) - 
            alpha1^2 * (phi_a0 - phi0 - derphi0 * alpha0)
    a = a / factor
    b = -alpha0^3 * (phi_a1 - phi0 - derphi0 * alpha1) + 
            alpha1^3 * (phi_a0 - phi0 - derphi0 * alpha0)
    b = b / factor
    
    alpha2 = (-b + sqrt(abs(b^2 - 3 * a * derphi0))) / (3.0 * a)
    phi_a2 = phi(x,  dx, alpha2)
    if ((phi_a2 <= (phi0 + c1 * alpha2 * derphi0)) && alpha2 >= 0) {
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
  
wass_div_opt <- function(X0, X1, lambda, a = NULL, b = NULL, niter = 1e3, step = 1e-1, 
                         tol = 1e-5, ...) {
  
  convergence <- function(old, new, tol) {
    diff <- abs(old - new) 
    vals <- diff /old
    vals[diff == 0] <- 0
    all(vals < tol)
  }
  
  n <- nrow(X0)
  m <- nrow(X1)
  
  if(is.null(a)) a <- rep(1/n, n)
  if(is.null(b)) b <- rep(1/m, m)
  # b_pot <- sinkhorn_geom(x = X1, y =X1,
  #                        a = b, b = b,
  #                        ...)
  # 
  # q_ <- b_pot$g
  
  a_hat <- a_tilde <- a_ <-  a
  a_old_ <- rep(0,length(a))
  
  norm_const <- norm_const_old <- 1
  
  pot_update <- step_ <- f <- NULL
  f_old <- Inf
  
  for (i in 1:niter) {
    
    step_ <- (i + 1)/2
    
    a_ <- (1 - 1/step_) * a_hat + 1/step_ * a_tilde
    
    a_old_ <- a_
    
    pot_update <- sinkhorn_geom(x = X0, y =X1,
                                a = a_, b = b, blur = lambda,
                                ...)
    
    f <- sum(a_ * pot_update$f) + sum(b * pot_update$g)
    
    if(convergence(f, f_old, tol)) break
    f_old <- f
    
    a_tilde <- a_tilde * exp(- step * step_ * pot_update$f)
    a_tilde <- a_tilde/sum(a_tilde)
    
    a_hat <- (1 - 1/step_) * a_hat + 1/step_ * a_tilde
    
    
  }
  
  return(a_)
}
