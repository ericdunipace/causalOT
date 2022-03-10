# conditional gradient descent
cg <- function(optimizer, verbose = TRUE) {
  stopifnot(inherits(optimizer, "cgOptimizer")) #needs to be R6 method
  
  optimizer$solve_G()
  # if (verbose) pb <- txtProgressBar(min = 0, max = floor(optimizer$get_niter()/10), style = 3)
  if (verbose) pb <- txtProgressBar(min = 0, max = optimizer$get_niter(), style = 3)
  for (i in 1:optimizer$get_niter()) {
    optimizer$solve_param()
    if (optimizer$converged() && i > 1) {
      if(verbose) message("\nConverged")
      break
    }
    # cat("iter:", i, ", objective: ", optimizer$f(),"\n")
    optimizer$solve_S()
    optimizer$step()
    # if (verbose && i %% 10 == 0) setTxtProgressBar(pb, i/10)
    if (verbose) setTxtProgressBar(pb, i)
  }
  # optimizer$solve_S()
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
                              calc_gamma = function() {
                                # lambdas <- 10^(-3:log10(private$lambda))
                                # weights <- lapply(lambdas, function(blur) {
                                #   optprob <- sinkhorn_geom(x = private$X1, y = private$X2,
                                #                            a = private$a,
                                #                            b = private$b, power = private$p,
                                #                            blur = blur, reach = private$sinkhorn_args$reach,
                                #                            diameter = private$sinkhorn_args$diameter,
                                #                            scaling = private$sinkhorn_args$scaling,
                                #                            truncate = private$sinkhorn_args$truncate,
                                #                            metric = "Lp", kernel = private$sinkhorn_args$kernel,
                                #                            cluster_scale=private$sinkhorn_args$cluster_scale,
                                #                            debias=FALSE,
                                #                            verbose=private$sinkhorn_args$verbose,
                                #                            backend=private$sinkhorn_args$backend);
                                #   raw_pi <- private$dual_to_primal(optprob$f, optprob$g, blur);
                                #   weights <- list(w0 = private$a,
                                #                   w1 = private$b,
                                #                   gamma = round_pi(raw_pi/sum(raw_pi), private$a, private$b),
                                #                   method = "Wasserstein", estimand = "ATT",
                                #                   args = list(penalty = blur,
                                #                               div.pen = private$lambda));
                                #   class(weights) <- "causalWeights";
                                #   return(weights)
                                # })
                                
                                # med.cost <- median(private$cost)
                                # cost_p  <- (private$cost)^(1/private$p)
                                # lambdas <- 10^(-1:log10(private$lambda))
                                # weights <- lapply(lambdas, function(blur) {
                                #   tplan <- approxOT::transport_plan_given_C(mass_x = private$a,
                                #                                             mass_y = private$b,
                                #                                             p = private$p,
                                #                                             cost = cost_p,
                                #                                             method = "sinkhorn",
                                #                                             epsilon = blur / med.cost,
                                #                                             niter = 1e4)
                                #   gamma <- matrix(0, private$n1, private$n2)
                                #   gamma[dist_2d_to_1d(tplan$from, tplan$to, private$n1, private$n2)] <- tplan$mass
                                #   weights <- list(w0 = private$a,
                                #                   w1 = private$b,
                                #                   gamma = gamma,
                                #                   method = "Wasserstein", estimand = "ATT",
                                #                   args = list(penalty = blur,
                                #                               div.pen = private$lambda));
                                #   class(weights) <- "causalWeights";
                                #   return(weights)
                                # })
                                # 
                                # n.boot <- 1e2
                                # K <- 10
                                # R <- 10
                                # eval.method <- match.arg(eval.method)
                                # wass.method <- "sinkhorn"
                                # wass.iter <- 1e3
                                # estimand <- "ATT"
                                # method <- "Wasserstein"
                                # solver <- "LBFGS"
                                # 
                                # weight <- eval_weights(weights, 
                                #                       args = list(x0 = private$X1,
                                #                                                x1 = private$X2,
                                #                                                x = rbind(private$X1,
                                #                                                          private$X2),
                                #                                                z = c(rep(0,private$n1),
                                #                                                      rep(1, private$n2)),
                                #                                                grid = lambdas,  
                                #                                                n.boot = n.boot,
                                #                                                K = K, 
                                #                                                R = R,
                                #                                                eval.method = eval.method,
                                #                                                wass.method = wass.method,
                                #                                                wass.iter = wass.iter,
                                #                                                sample_weight = list(a = private$a,
                                #                                                                     b = private$b,
                                #                                                                     total = renormalize(c(private$a * private$n1,
                                #                                                                                         private$b * private$n2))),
                                #                                                estimand = "ATT", 
                                #                                                method = "Wasserstein", 
                                #                                                solver = private$search, 
                                #                                                metric = private$metric,
                                #                                                p = private$p, 
                                #                                                cost = private$cost, 
                                #                                                add.joint = TRUE,
                                #                                                add.margins = FALSE, 
                                #                                                joint.mapping = FALSE,
                                #                                                neg.weights = FALSE,
                                #                                                cgd = FALSE, 
                                #                                                verbose = isTRUE(verbose)))
                                # private$gamma <- weight$gamma
                                # private$gamma.lambda <- weight$args$penalty
                                
                                med.cost <- median(private$cost)
                                cost_p  <- (private$cost)^(1/private$p)
                                min.lambda <- max(private$cost)/abs(.Machine$double.min.exp) + 0.01
                                tplan <- approxOT::transport_plan_given_C(mass_x = private$a,
                                                                          mass_y = private$b,
                                                                          p = private$p,
                                                                          cost = cost_p,
                                                                          method = "sinkhorn",
                                                                          epsilon = min.lambda / med.cost,
                                                                          niter = 1e2)
                                mass <- tplan$mass
                                mass[mass < 0] <- 0
                                private$gamma <- matrix(0, private$n1, private$n2)
                                private$gamma[dist_2d_to_1d(tplan$from, tplan$to, private$n1, private$n2)] <- renormalize(mass)
                                
                                # blur <- 1e-3
                                # optprob <- sinkhorn_geom(x = private$X1, y = private$X2,
                                #             a = private$a,
                                #             b = private$b, power = private$p,
                                #             blur = blur, reach = private$sinkhorn_args$reach,
                                #             diameter = private$sinkhorn_args$diameter,
                                #             scaling = private$sinkhorn_args$scaling,
                                #             truncate = private$sinkhorn_args$truncate,
                                #             metric = "Lp", kernel = private$sinkhorn_args$kernel,
                                #             cluster_scale = private$sinkhorn_args$cluster_scale,
                                #             debias=FALSE,
                                #             verbose=private$sinkhorn_args$verbose,
                                #             backend=private$sinkhorn_args$backend)
                                # raw_pi <- private$dual_to_primal(optprob$f, optprob$g, blur)
                                # Total <- sum(raw_pi)
                                # if (is.character(private$cost)) {
                                #   private$cost <- cost_calc_lp(private$X1, private$X2,
                                #                                p = private$p, direction = "rowwise")^private$p
                                # }
                                # private$gamma <- round_pi(optprob$f + blur * log(private$a), 
                                #                           optprob$g + blur * log(private$b), 
                                #                           private$cost, blur, 
                                #                           private$a, private$b)
                                # private$gamma <- round_pi(raw_pi/Total, 
                                #                           private$a, private$b)
                                
                              },
                              calc_scm_bary_proj = function() {

                                op <- list()
                                # op$obj <-

                              },
                              converged = function() {
                                private$f_val <- self$f()
                                f_val_diff <- abs(private$f_val - private$f_val_old)
                                
                                # nan.check <- isTRUE(any(is.nan(private$a)) || is.nan(private$f_val))
                                
                                nan.check <- isTRUE(is.nan(private$f_val))
                                
                                # converged <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol)  ||
                                #   isTRUE(sum(abs(private$a_old - private$a)) < private$tol) ||
                                #   nan.check
                                
                                relsame <- isTRUE(f_val_diff / abs(private$f_val_old) < private$tol)  || nan.check
                                private$converged.count <- if(relsame) {
                                  private$converged.count + as.integer(relsame)
                                } else {
                                  0L
                                }
                                  
                                if(private$search == "LBFGS") {
                                  converged <- isTRUE(private$converged.count >= 2)
                                } else {
                                  converged <- isTRUE(private$converged.count >= 1)
                                }
                                
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
                              get_a = function() {
                                return(private$a)
                              },
                              get_niter = function() {
                                return(private$niter)
                              },
                              get_weight = function() {
                                
                                # fmat <- matrix(c(private$f_pot$numpy()), private$n1, private$n2)
                                # gmat <- matrix(c(private$g_pot$numpy()), private$n1, private$n2, byrow = TRUE)
                                # eta <- (fmat + gmat - private$cost)/private$lambda
                                # pi_raw <- exp(eta) * matrix(private$a, private$n1, private$n2) * matrix(private$b, private$n1, private$n2, byrow = TRUE)
                                # pi <- round_pi(pi_raw, private$a, private$b)
                                # return( pi )
                                if(!is.matrix(private$gamma)) self$calc_gamma()
                                return(private$gamma)
                              },
                              return_cw = function(...) {
                                # estimand <- match.arg(estimand, c("ATT","ATC","ATE"))
                                param <- self$get_param()
                                out <- list(w0 = private$a,
                                            w1 = private$b,
                                            gamma = self$get_weight(),
                                            args = list(
                                              dual = list(f = param$f,
                                                          g = param$g),
                                              solver = private$prog_solver,
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
                                              conditional.gradient = if(private$search %in% c("armijo")){TRUE}else{FALSE},
                                              search = private$search),
                                            estimand = "ATT",
                                            method = "Wasserstein")
                                class(out) <- "causalWeights"
                                return(out)
                              },
                              solve_G = function() {
                                # self$solve_param()
                                private$f_val_old <- private$sinkhorn_args$diameter
                                private$a_old <- private$a
                                if (private$search == "armijo") {
                                  private$op <- private$op_update(f = rep(1, private$n1),
                                                                  g = private$g_pot,
                                                                  a = private$a,
                                                                  b = private$b,
                                                                  op = private$op)
                                  private$a <- private$solver(private$op)$sol
                                }
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
                                if (private$search == "armijo" ) {
                                  private$op <- private$op_update(f =  self$get_param()$f,
                                                                  g = private$g_pot,
                                                                  a = private$a,
                                                                  b = private$b,
                                                                  op = private$op)
                                  private$S <- private$solver(private$op)$sol
                                  # private$S <- as.numeric(private$f_pot == min(private$f_pot))
                                  
                                  # grad <- self$df()
                                  # private$S <- renormalize(as.numeric(grad == min(grad)))
                                } else if (private$search == "LBFGS") {
                                  private$scheduler$step(private$f_val)
                                }
                                
                              },
                              step = function() {
                                
                                private$cur_iter <- private$cur_iter + 1L
                                
                                
                                
                                old_a <- private$a
                                
                                if (private$search == "armijo") {
                                  df_val <- c(self$df())
                                  f_val <- self$f()
                                  
                                  deltaG <- private$S - old_a
                                  derphi0 <- sum(deltaG * df_val)
                                  f <- function(x, dx, alpha, ...) {
                                    
                                    # xnew <- simplex_proj(x)
                                    # proj <- x * exp( dx * alpha)
                                    # private$a <- proj / sum(proj)
                                    private$a <- x + dx * alpha
                                    private$a[private$a < 0] <- 0
                                    private$a <- renormalize(private$a)
                                    l_a <- log(private$a)
                                    l_a[is.infinite(l_a)] <- (-.Machine$double.xmax)
                                    private$pydat$l_at$data <- private$torch$DoubleTensor(l_a)$contiguous()$to(private$device)
                                    private$pydat$at$data <- private$torch$softmax(private$pydat$l_at$detach(), 0L)$to(private$device)
                                    self$solve_param()
                                    
                                    
                                    loss <- self$f() #self$f()
                                    # if (loss < 0) return(f_val)
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
                                      private$a[private$a < 0] <- 0
                                      private$a <- renormalize(private$a)
                                      if(private$python_running) {
                                        l_a <- log(private$a)
                                        l_a[is.infinite(l_a)] <- (-.Machine$double.xmax)
                                        private$pydat$l_at$data <- private$torch$DoubleTensor(l_a)$contiguous()$to(private$device)
                                        private$pydat$at$data <- private$torch$softmax(private$pydat$l_at$detach(), 0L)$to(private$device)
                                      }
                                    }
                                    
                                  } else {
                                    private$a <- old_a
                                  }
                                } else if (private$search == "mirror") {
                                  df_val <- c(self$df())
                                  eta <- private$stepsize / sqrt(private$cur_iter)
                                  prop <- old_a * exp(-eta * df_val)
                                  private$a <- prop/sum(prop)
                                  
                                  if (private$python_running) {
                                    l_a <- log(private$a)
                                    l_a[is.infinite(l_a)] <- (-.Machine$double.xmax)
                                    private$pydat$l_at$data <- private$torch$DoubleTensor(l_a)$contiguous()
                                    private$pydat$at$data <- private$torch$softmax(private$pydat$l_at$detach(), 0L)
                                  }
                                  
                                } else if (private$search == "mirror-accelerated") {
                                  df_val <- c(self$df())
                                  step_ <- (private$cur_iter + 1)/2
                                  
                                  
                                  private$a_tilde <- private$a_tilde * exp(- private$stepsize * step_ * df_val)
                                  private$a_tilde <- private$a_tilde/sum(private$a_tilde)
                                  
                                  private$a_hat <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                  
                                  step_ <- step_ + 0.5
                                  
                                  private$a <- (1 - 1/step_) * private$a_hat + 1/step_ * private$a_tilde
                                  
                                  if (private$python_running) {
                                    l_a <- log(private$a)
                                    l_a[is.infinite(l_a)] <- (-.Machine$double.xmax)
                                    private$pydat$l_at$data <- private$torch$DoubleTensor(l_a)$to(private$device)
                                    private$pydat$at$data <- private$torch$softmax(private$pydat$l_at$detach(), 0L)$to(private$device)
                                  }
                                  
                                  # update parameters next
                                  
                                } else if (private$search == "LBFGS" || private$search == "pgd") {
                                  private$optimizer$zero_grad()
                                  private$optimizer$step(private$closure)
                                  
                                  private$pydat$at <- private$torch$softmax(private$pydat$l_at$detach(), 0L)
                                  private$a <- c(private$pydat$at$cpu()$numpy())
                                  
                                  if(private$search == "pgd") {
                                    private$op <- private$op_update(f =  self$get_param()$f,
                                                                    g = private$g_pot,
                                                                    a = private$a,
                                                                    b = private$b,
                                                                    op = private$op)
                                    private$a <- private$solver(private$op)$sol
                                    l_a <- log(private$a)
                                    l_a[is.infinite(l_a)] <- (-.Machine$double.xmax)
                                    private$pydat$l_at$data <- private$torch$DoubleTensor(l_a)$contiguous()$to(private$device)
                                    private$pydat$at$data <- private$torch$softmax(private$pydat$l_at$detach(), 0L)$to(private$device)
                                    
                                  }
                                }
                                
                                
                              },
                              initialize = function(X1, X2, 
                                                    cost,
                                                    prog_solver = supported.solvers(),
                                                    lambda = 100,
                                                    add.margins = FALSE,
                                                    metric = dist.metrics(),
                                                    power = 2,
                                                    niter = 1000,
                                                    tol = 1e-7,
                                                    search = c("LBFGS",
                                                               "pgd",
                                                               "armijo",
                                                               "mirror",
                                                               "mirror-accelerated"),
                                                    stepsize = 1e-1,
                                                    sample_weight = NULL,
                                                    reach = NULL,
                                                    diameter = NULL,
                                                    scaling = 0.5, truncate = 5,
                                                    kernel = NULL,
                                                    cluster_scale = NULL, 
                                                    debias = TRUE, 
                                                    verbose = FALSE, backend='auto',
                                                    balance.function.formula = NULL,
                                                    balance.function.delta = NULL,
                                                    ...
                              ) {
                                metric    <- match.arg(metric, dist.metrics())
                               
                                private$penalty <- "entropy"
                                if (missing(niter) || length(niter) == 0) niter <- 1000
                                if (missing(tol) || length(tol) == 0) tol <= 1e-7 
                                private$search <- match.arg(search)
                                private$python_running <- FALSE
                                
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
                                
                                
                                if(!is.null(balance.function.formula) && !is.na(balance.function.formula)) {
                                  private$prog_solver <- match.arg(prog_solver)
                                  private$solver <- switch(private$prog_solver,
                                                           # "cplex" = cplex_solver,
                                                           # "gurobi" = gurobi_solver,
                                                           "mosek" = mosek_solver,
                                                           "osqp" = osqp_solver)
                                  if(private$search != "pgd" && private$search != "armijo") private$search <- "armijo"
                                  
                                  # get balance functions
                                  if (is.null(balance.function.delta)) balance.function.delta <- 0.05
                                  
                                  form <- form_all_squares(balance.function.formula, colnames(private$X2))
                                  
                                  form.temp <- as.character(form[length(form)])
                                  form <- as.formula(paste0("~ 0 +", form.temp))
                                  
                                  BC <- list(source = model.matrix(formula(form), data.frame(private$X1)),
                                             target = model.matrix(formula(form), data.frame(private$X2)),
                                             K = balance.function.delta)
                                  
                                  if ( all(BC$source[,1] == 1)) BC$source <- BC$source[,-1]
                                  if ( all(BC$target[,1] == 1)) BC$target <- BC$target[,-1]
                                  
                                  private$op <- switch(private$search,
                                                       "pgd" = qp_proj(f = rep(1,private$n1), 
                                                                       g = rep(1, private$n2), 
                                                                       a = private$a, 
                                                                       b = private$b, 
                                                                       BC = BC),
                                                       "armijo" = lp_min_constraint(f = rep(1,private$n1), 
                                                                                    g = rep(1, private$n2), 
                                                                                    a = private$a, 
                                                                                    b = private$b, 
                                                                                    BC = BC)
                                  )
                                  private$op_update <- switch(private$search,
                                                              "pgd" = qp_proj_update,
                                                              "armijo" = lp_min_constraint_update)
                                } else {
                                  if(private$search == "pgd" || private$search == "armijo") {
                                    warning("Can't do projected gradient descent (pgd) or armijo without balance functions! Setting optimizer to LBFGS")
                                    private$search <- "LBFGS"
                                  }
                                }
                                
                                
                                private$lambda <- lambda
        
                                private$sinkhorn_args <- list(blur = private$lambda, reach = reach, 
                                                              diameter = diameter,
                                                              scaling = scaling, 
                                                              truncate = truncate,
                                                              metric = "Lp", kernel = kernel,
                                                              cluster_scale = cluster_scale, 
                                                              debias = TRUE, # debias, 
                                                              verbose = verbose, 
                                                              backend = backend)
                                # sets up python function
                                
                                if(private$search == "LBFGS" || private$search == "pgd") {
                                  private$python_running <- TRUE
                                  private$np <- reticulate::import("numpy", convert = TRUE)
                                  private$torch <- reticulate::import("torch", convert = TRUE)
                                  private$geomloss <- reticulate::import("geomloss", convert = TRUE)
                                  use_cuda <- private$torch$cuda$is_available()
                                  private$device <- private$torch$device(if(use_cuda){"cuda"} else {"cpu"})
                                }
                                private$converged.count <- 0L
                                private$specific_initialize()
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
                           "closure" = "function",
                           "converged.count" = "integer",
                           "cost" = "numeric",
                           "cost_idx" = "numeric",
                           "cur_iter" = "integer",
                           "d" = "numeric",
                           "device" = "python.builtin.module",
                           "dual_to_primal" = function(f, g, blur) {
                             if(is.character(private$cost)) {
                               private$cost <- cost_calc_lp(private$X1, private$X2,
                                                                                         p = private$p, direction = "rowwise")^private$p
                             }
                             return(tcrossprod(private$a, private$b) * exp((matrix(f, private$n1, private$n2) +
                                                                      matrix(g, private$n1, private$n2, byrow = TRUE) -
                                                                      private$cost)/blur))
                             
                           },
                           "f_pot" = "numeric",
                           "f_val" = "numeric",
                           "f_val_old" = "numeric",
                           "G" = "numeric",
                           "G_old" = "numeric",
                           "g_pot" = "numeric",
                           "gamma" = "matrix",
                           "geomloss" = "python.builtin.module",
                           "lambda" = "numeric",
                           "maxcost" = "numeric",
                           "metric" = "character",
                           "otModel" = "python.builtin.object",
                           "optimizer" = "python.builtin.object",
                           "n1" = "numeric",
                           "n2" = "numeric",
                           "niter" = "numeric",
                           "np" = "python.builtin.module",
                           "op" = "list",
                           "op_update" = "character",
                           "p" = "numeric",
                           "prog_solver" = "character",
                           "pydat" = "list",
                           "penalty" = "character",
                           "python_running" ="logical",
                           "S" = "numeric",
                           "search" = "character",
                           "scheduler" = "python.builtin.object",
                           "specific_initialize" = function(){},
                           "stepsize" = "numeric",
                           "sinkhorn_args" = "list",
                           "solver" = "function",
                           "tol" = "numeric",
                           "torch" = "python.builtin.module",
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
                               f = function() {
                                 return(
                                   c(private$torch$add(private$torch$dot(private$f_pot, private$pydat$at$detach()), 
                                               private$torch$dot(private$g_pot, private$pydat$bt))$cpu()$numpy())
                                 )
                               },
                               df = function() {
                                 return(private$f_pot$cpu()$numpy())
                               },
                               get_param = function() {
                                 return(
                                   list(f = as.numeric(private$f_pot$cpu()$numpy()),
                                       g = as.numeric(private$g_pot$cpu()$numpy()))
                                 )
                               },
                               solve_param = function() {
                                 
                                 sol <- private$otModel$forward(private$pydat$at$detach(), 
                                                                private$pydat$xt, 
                                                                private$pydat$bt, 
                                                                private$pydat$yt)
                                 private$f_pot <- sol[[1]]$squeeze()
                                 private$g_pot <- sol[[2]]$squeeze()
                                 # sol <- sinkhorn_geom(x = private$X1, y = private$X2, 
                                 #                      a = private$a, 
                                 #                      b = private$b, power = private$p, 
                                 #                      blur = private$lambda, reach = private$sinkhorn_args$reach, 
                                 #                      diameter = private$sinkhorn_args$diameter,
                                 #                      scaling = private$sinkhorn_args$scaling, 
                                 #                      truncate = private$sinkhorn_args$truncate,
                                 #                      metric = "Lp", kernel = private$sinkhorn_args$kernel,
                                 #                      cluster_scale=private$sinkhorn_args$cluster_scale, 
                                 #                      debias=private$sinkhorn_args$debias, 
                                 #                      verbose=private$sinkhorn_args$verbose, 
                                 #                      backend=private$sinkhorn_args$backend)
                                 # 
                                 # private$f_pot <- sol$f
                                 # private$g_pot <- sol$g
                                 
                                 
                                 
                               }
                               
                           ),
                           private = list(
                             "specific_initialize" = function() {
                               private$python_running <- TRUE
                               private$np <- reticulate::import("numpy", convert = TRUE)
                               private$torch <- reticulate::import("torch", convert = TRUE)
                               private$geomloss <- reticulate::import("geomloss", convert = TRUE)
                               
                               use_cuda <- private$torch$cuda$is_available()
                               dtype <- if(use_cuda){private$torch$cuda$DoubleTensor} else {private$torch$DoubleTensor}
                               private$device <- private$torch$device(if(use_cuda){"cuda"} else {"cpu"})
                               
                               if (private$sinkhorn_args$backend == "tensorized" || (private$n1 <= 5000 && private$n2 <= 5000 && private$sinkhorn_args$backend != "multiscale" && private$sinkhorn_args$backend != "online")) {
                                 if (private$p == 2) {
                                   cost <- private$geomloss$utils$squared_distances
                                 } else if (private$p == 1) {
                                   reticulate::source_python(file = lp_python_path)
                                   cost <- l1_loss
                                 } else {
                                   # reticulate::source_python(file = lp_python_path)
                                   # cost <- lp_loss
                                   cost <- paste0("Sum(Pow(X-Y,", private$p,"))")
                                   private$sinkhorn_args$backend <- "online"
                                   private$p <- 1L
                                 }
                               } else if (private$n1 > 5000 || private$n2 > 5000 || private$sinkhorn_args$backend == "multiscale" || private$sinkhorn_args$backend == "online") {
                                 if(private$sinkhorn_args$backend == "tensorized") private$sinkhorn_args$backend <- "online"
                                 if (private$p == 2) {
                                   cost <- "SqDist(X,Y)"
                                 } else if (private$p == 1) {
                                   cost <- "Sum(Abs(X - Y))"
                                 } else {
                                   cost <- paste0("Sum(Pow(X-Y,", private$p,"))")
                                   private$p <- 1L
                                 }
                                 # pykeops <- reticulate::import("pykeops", convert = TRUE)
                                 # pykeops$clean_pykeops()
                               } else {
                                 cost <- NULL
                               }
                               
                               # if(search != "mirror" && search != "mirror-accelerated") {
                               
                               
                               private$pydat <- list()
                               private$pydat$xt <- dtype(private$np$array(private$X1))$contiguous()
                               private$pydat$yt <- dtype(private$np$array(private$X2))$contiguous()
                               private$pydat$at <- dtype(private$a)$contiguous()
                               private$pydat$l_at <- private$torch$autograd$Variable(dtype(log(private$a))$contiguous(), requires_grad = TRUE)
                               private$pydat$bt <- dtype(private$b)$contiguous()
                               
                               # private$pydat$l_at <- private$pydat$l_at$to(device)
                               # private$pydat$at <- private$pydat$at$to(device)
                               
                               # sets up python function
                               
                               private$otModel <- private$geomloss$SamplesLoss("sinkhorn", p = private$p, 
                                                                               blur = private$sinkhorn_args$blur, 
                                                                               reach = private$sinkhorn_args$reach,
                                                                               diameter = private$sinkhorn_args$diameter, 
                                                                               scaling = private$sinkhorn_args$scaling, 
                                                                               cost = cost, kernel = private$sinkhorn_args$kernel,
                                                                               cluster_scale = private$sinkhorn_args$cluster_scale,
                                                                               debias = private$sinkhorn_args$debias,
                                                                               potentials = TRUE,
                                                                               verbose = private$sinkhorn_args$verbose,
                                                                               backend = private$sinkhorn_args$backend)
                               if(private$search == "LBFGS" || private$search == "pgd") {
                                 private$optimizer <- private$torch$optim$LBFGS(params = list(private$pydat$l_at),
                                                                                lr = private$stepsize
                                                                                , line_search_fn = "strong_wolfe"
                                 )
                                 private$scheduler <- private$torch$optim$lr_scheduler$ReduceLROnPlateau(optimizer = private$optimizer, mode = "min", patience = 0L)
                                 
                                 private$closure <- function() { #needed for LBFGS search
                                   private$optimizer$zero_grad()
                                   private$pydat$at <- private$torch$softmax(private$pydat$l_at, 0L)
                                   pot <- private$otModel$forward(private$pydat$at$detach(),
                                                                  private$pydat$xt,
                                                                  private$pydat$bt,
                                                                  private$pydat$yt)
                                   loss <- private$torch$add(private$torch$dot(pot[[1]]$squeeze(), private$pydat$at), 
                                                             private$torch$dot(pot[[2]]$squeeze(), private$pydat$bt))
                                   loss$backward()
                                   return(loss)
                                 }
                               }
                             }
                           )
                           
                           
  )
}

#L2 based divergence
{
  {
    wassDivL2 <- R6::R6Class("wassDivL2", 
                              inherit = wassDiv,
                              public = list(
                                f = function() {
                                  return(
                                    c(private$f_val)
                                  )
                                },
                                get_param = function() {
                                  return(
                                    list(f = private$f_pot,
                                         g = private$g_pot)
                                  )
                                },
                                solve_param = function() {
                                  if (private$cur_iter == 0) {
                                    sol <- private$otModel$forward()
                                  } else {
                                    sol <- private$otModel$update_a(private$a)
                                  }
                                  
                                  private$f_pot <- sol$f
                                  private$g_pot <- sol$g
                                  private$f_val <- sol$loss
                                }
                                
                              ),
                              private = list(
                                "specific_initialize" = function() {
                                  
                                  private$otModel <- otDualL2$new(x = private$X1, y = private$X2,
                                                                  a = private$a,
                                                                  b = private$b,
                                                                  p = private$p,
                                                                  lambda = private$lambda,
                                                                  solver = private$prog_solver,
                                                                  debias = TRUE,
                                                                  cost = if(!is.null(private$cost) && !is.character(private$cost)) {
                                                                    private$cost^(1/private$p)
                                                                  } else {
                                                                    NULL
                                                                  },
                                                                  control = list(maxit = private$niter,
                                                                                 lmm = 20L))
                                  if(private$search == "LBFGS" || private$search == "pgd") {
                                    # sets up python function
                                    private$python_running <- TRUE
                                    private$np <- reticulate::import("numpy", convert = TRUE)
                                    private$torch <- reticulate::import("torch", convert = TRUE)
                                    private$geomloss <- reticulate::import("geomloss", convert = TRUE)
                                    
                                    private$pydat <- list()
                                    private$pydat$at <- private$torch$DoubleTensor(private$a)$contiguous()
                                    private$pydat$bt <- private$torch$DoubleTensor(private$b)$contiguous()
                                    private$pydat$l_at <- private$torch$autograd$Variable(private$torch$DoubleTensor(log(private$a))$contiguous(), requires_grad = TRUE)
                                    
                                    private$optimizer <- private$torch$optim$LBFGS(params = list(private$pydat$l_at),
                                                                                   lr = private$stepsize
                                                                                   , line_search_fn = "strong_wolfe"
                                    )
                                    private$closure <- function() { #needed for LBFGS search
                                      private$optimizer$zero_grad()
                                      private$pydat$at <- private$torch$softmax(private$pydat$l_at, 0L)
                                      pot <- private$otModel$update_a(c(private$pydat$at$detach()$numpy()))
                                      ft <- private$torch$DoubleTensor(pot$f)
                                      gt <- private$torch$DoubleTensor(pot$g)
                                      # potentials_loss <- private$torch$add(private$torch$dot(ft, private$pydat$at), 
                                      #                                                        private$torch$dot(gt, private$pydat$bt))
                                      # loss <- private$torch$sub(private$torch$DoubleTensor(list(pot$loss)), potentials_loss$detach())
                                      # loss <- private$torch$add(loss, potentials_loss)
                                      loss <- private$torch$add(private$torch$dot(ft, private$pydat$at), 
                                                                           private$torch$dot(gt, private$pydat$bt))
                                      loss$backward()
                                      return(loss)
                                    }
                                  }
                                }
                              )
                              
                              
    )
  }
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
                                                              qp_solver = supported.solvers(),
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
                                                           # "cplex" = cplex_solver,
                                                           # "gurobi" = gurobi_solver,
                                                           "mosek" = mosek_solver,
                                                           "osqp" = osqp_solver)
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
                          qp_solver = supported.solvers(),
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
                       # "cplex" = cplex_solver,
                       # "gurobi" = gurobi_solver,
                       "mosek" = mosek_solver,
                       "osqp" = osqp_solver)
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
                                                      qp_solver = supported.solvers(),
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
                                                # "cplex" = cplex_solver,
                                                # "gurobi" = gurobi_solver,
                                                "mosek" = mosek_solver,
                                                "osqp" = osqp_solver)
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
  alpha1 = -(derphi0) * alpha0 ^ 2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
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
