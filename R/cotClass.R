# setup simpler L2 and ent distances with osqp and LBFGS, respectively
# for ent dist, consider using torch so have GPU access with keops...

# cot classes

setOldClass(c("gridSearchClass","R6"))
setOldClass("balanceFunction")
setOldClass(c("OT","R6"))
setOldClass("torch_tensor")
COT <- R6::R6Class("COT",
  # inherit = gridSearchClass,
  public = list(
   a = "numeric",
   b = "numeric",
   n = "numeric",
   m = "numeric",
   ot = "OT",
   penalty.fun = "character",
   bf = "balanceFunction",
   runbf = "logical",
   param = "torch_tensor",
   optim = "torch_optimizer",
   scheduler = "lr_scheduler",
   niter = "integer",
   tol = "numeric",
   evalBoot = function(w_tilde, b_tilde) {
     # browser()
     torch::with_no_grad({
       self$ot$a <- w_tilde
       self$ot$b <- b_tilde
       loss <- energy_dist(self$ot)$item()
       self$ot$b <- self$b
     }) 
     return(
       loss
     )
   },
   initialize = function(source, target,
                         a = NULL, b = NULL,
                         options = list(NULL)) {
     # browser()
     if (missing(source) || missing(target)) stop("Must provide source and target")
     if (!inherits(options, "cotOptions")) options <- do.call(cotOptions, options)

     if (is.null(options$lambda)) {
       lambda <- 10
     } else {
       lambda <- options$lambda
     }
     self$ot <- OT$new(x = source, y = target, a = a, b = b,
                   penalty = lambda, cost = options$cost,
                   p = options$p, debias = options$debias,
                   tensorized = options$cost.online,
                   diameter=NULL)
     use_cuda <- torch::cuda_is_available() && torch::cuda_device_count()>1
     if (use_cuda) {
       device <-  torch::torch_device("cuda")
       if (!self$ot$tensorized) {
         rkeops::compile4float64()
         rkeops::compile4gpu()
         rkeops::use_gpu()
       }
     } else {
       device <-  torch::torch_device("cpu")
       if (!self$ot$tensorized) {
         rkeops::compile4float64()
         rkeops::use_cpu()
       }
     }
     self$a <- as.numeric(self$ot$a$to(device = torch::torch_device("cpu")))
     self$b <- as.numeric(self$ot$b$to(device = torch::torch_device("cpu")))
     self$n <- as.numeric(self$ot$n)
     self$m <- as.numeric(self$ot$m)
     
     if (!is.null(options$balance.formula)) {
       if (is.null(colnames(source))) colnames(source) <- paste0("X", 1:ncol(source))
       if (is.null(colnames(target))) colnames(target) <- colnames(source)
       tbf <- terms(formula(options$balance.formula))
       attr(tbf, "intercept") <- 0
       source.bf <- model.matrix(tbf, data.frame(source))
       target.bf <- model.matrix(tbf, data.frame(target))
       self$bf <- balanceFunction$new(source = source.bf, 
                                      target = target.bf,
                                      a = as.numeric(self$ot$a),
                                      b = as.numeric(self$ot$b),
                                      delta = options$delta)
     } else {
       self$bf <- list(NULL)
     }
     
     self$runbf <- FALSE
     if(inherits(self$bf, "balanceFunction")) self$runbf <- TRUE
     
     self$param  <- torch::torch_zeros(self$ot$n - 1,
                                        dtype = torch::torch_double(),
                                       device = device,
                                        requires_grad = TRUE)
     # self$param  <- torch::torch_zeros(self$ot$n,
     #                                   dtype = torch::torch_double(),
     #                                   requires_grad = TRUE)
     if ( inherits(options$torch.optimizer, "torch_optimizer_generator") ) {
       private$torch_optim <- list(
         fun = options$torch.optimizer,
         args = options$solver.options
       )
       self$optim <- do.call(private$torch_optim$fun,
                             c(list(params = list(self$param)), private$torch_optim$args))
       
       if (inherits(self$optim, "optim_lbfgs")) {
         cb <- causalOT:::.cubic_interpolate
         tmpfun <- get(".cubic_interpolate", envir = asNamespace("torch"))
         environment(cb) <- environment(tmpfun)
         attributes(cb) <- attributes(tmpfun)
         assignInNamespace(".cubic_interpolate", cb, "torch" )
         
         ln_srch <- private$torch_optim$args$line_search_fn
         no_ls <- (is.null(ln_srch) || is.na(ln_srch) || ln_srch != "strong_wolfe")
         if(no_ls && lambda != 0 ) {
           warning(" Torch's LBFGS doesn't work well without 'strong_wolfe' line search on this problem. Specify it with line_search_fn = 'strong_wolfe' in the options argument.")
         }
       }
     } else {
       private$torch_optim <- list(
         fun = NULL,
         args = NULL
       )
     }

     if(inherits(options$torch.scheduler, "lr_scheduler") &&
        inherits(self$optim, "torch_optimizer") ) {
       self$scheduler <- options$torch.scheduler(c(list(optimizer = self$optim),
                                                   options$scheduler.options)
       )
     } else {
       self$scheduler <- list(step = function(value){invisible()})
     }
     self$penalty.fun <- options$penalty.function
     
     self$niter = options$niter
     self$tol = options$tol
     # private$sinkhorn_cot = torch::jit_trace(sinkhorn_cot, self$ot, self$niter, self$tol)

     # run setup for various methods
     if (self$penalty.fun == "L2") {
       private$L2_setup(options)
     } else if (self$penalty.fun == "entropy" && !self$ot$debias && self$ot$tensorized) {
       private$solver.options <- lbfgs3c_control(options$solver.options)
     } else if (self$penalty.fun == "entropy" && 
                ( self$ot$debias || !self$ot$tensorized) ) {
       private$forward <- private$cot_forward
     } #else {
     #   stop("Condition not found in COT setup. Please report this bug")
     # }
     
     # fix bincount function
     tb <- causalOT:::cot_torch_bincount
     tmpfun <- get("torch_bincount", envir = asNamespace("torch"))
     environment(tb) <- environment(tmpfun)
     attributes(tb) <- attributes(tmpfun)
     assignInNamespace("torch_bincount", tb, "torch" )
     
     return(invisible(self))
   }
  ),
  active = list(
    penalty = function(value) { # sets penalty values, and returns them if wanted
      if (missing(value)) {
        if(self$runbf) {
          return(c(lambda = self$ot$penalty,
                      delta = self$bf$delta))
        } else {
          return(c(lambda = self$ot$penalty ))
        }
      } else {
        if (is.list(value)) {
          if (!all(names(value) %in% c("lambda", "delta")) ) {
            stop("If penalties are provided as a list, must be with names in c('lamba', 'delta')")
          }
          lambda <- value$lambda
          delta <- value$delta
        } else {
          if (any(!is.null(names(value))) ) {
            if (!all(names(value) %in% c("lambda", "delta")) ) {
              stop("If penalties are provided as a named vector, must be with names in c('lamba', 'delta')")
            }
            lambda <- value["lambda"]
            delta <- value["delta"]
            names(lambda) <- NULL
            names(delta) <- NULL
          } else {
            lambda <- value[1]
            delta <- NULL
            if(length(value) > 1) warning("For unnamed vectors, only the first number is used to set the OT penalty, lambda. Other values are ignored.")
          }
          
        }
        if(!is.na(lambda) && !is.null(lambda) ) self$ot$penalty <- lambda
        if(!is.na(delta) && !is.null(delta) && self$runbf) self$bf$delta <- delta
        
        if (inherits(private$torch_optim$fun, "torch_optimizer_generator")) {
          self$optim <- do.call(private$torch_optim$fun,
                                c(list(params = self$param),
                                  private$torch_optim$args) )
        }
      }
    },
    weight = function(value) { # sets weights from the outside!
      # browser()
      if(missing(value)) {
        return(private$transform())
      } else {
        n_param <- length(self$param) + 1
        # n_param <- length(self$param)
        n_value  <- length(value)
        if(n_param != n_value) {
          stop("Length of weight for warm start is not equal to length of weight parameter.")
        }
        # browser()
        torch::autograd_set_grad_mode(enabled = FALSE)
        
        temp_log <- private$inv_transform(value)
        self$param$copy_(as.numeric(temp_log))
        
        torch::autograd_set_grad_mode(enabled = TRUE)
        
        # if (inherits(private$torch_optim$fun, "torch_optimizer_generator")) {
        #   self$optim <- do.call(private$torch_optim$fun,
        #                         c(list(params = self$param),
        #                           private$torch_optim$args) )
        # }
        return(invisible(self))
      }

    }
  ),
  private = list(
    delta.idx = "numeric",
    forward = "function",
    solver = "R6",
    solver.options = "list",
    transform = function() {
      # torch::with_no_grad({
      #   m <- self$param$mean()
      #   self$param$sub_(m)
      # }) 
      full_param <- torch::torch_cat(
        list(torch::torch_zeros(1, dtype = torch::torch_double()),
             self$param)
      )
      mirror_softmax(full_param)
      # mirror_softmax(self$param)
    },
    torch_optim = "function",
    inv_transform = function(value) {
      min_neg <- round(log(.Machine$double.xmin)) - 50.0
      logs <- torch::torch_log(value)
      logs[logs< min_neg] <- min_neg
      logs <- logs[2:length(logs)] - logs[1]
      # logs <- logs - logs[1]
      return(logs)
    }
  )
)
COT$set("public", "solve",
  function(penalty, w) {
    
    # get previous weights for warm start
    # if (length(w) > 0) {
    #   w <- w[[1]]
    #   n_param <- length(self$param)
    #   if (length(w) == n_param) {
    #     self$weight <- w
    #   } else {
    #     warning("Warm start not used. Proposed starting parameter of different length than params")
    #   }
    # }
    
    #update current penalty parameters
    if(missing(penalty) || any(is.na(penalty)) ) {
      stop("Penalty term must be given and be a named list or named vector")
    } else {
      self$penalty <- penalty
    }
    
    # estimate weights
    private$optimize_weights()

    return(as.numeric(self$weight))
  }
)
COT$set("private", "optimize_weights",
  function() {
    lambda <- self$ot$penalty
    
    if(lambda == 0 && !self$runbf) { # all methods for 0 with BF
      private$NNM_opt()
    } else if (self$ot$debias || (self$penalty.fun == "entropy" && !self$ot$tensorized)) { #sinkhorn divergence
      if (is.infinite( lambda ) ) {
        private$forward <- private$energy_dist_forward
      } else {
        private$forward <- private$cot_forward
      }
      private$torch_optimizer_loop()
    } else if (is.infinite(lambda) && !self$ot$debias) { #L2 and ent ot distance, infinite
      self$weight <- rep(1.0 / self$ot$n, self$ot$n)
    } else if (!self$ot$debias && self$ot$tensorized) {
      if(self$penalty.fun == "L2" || (lambda == 0 && self$runbf)) {
        if(!inherits(private$solver,  "osqp_model") ) private$L2_setup(list(NULL))
        private$L2_optimizer()
      } else if (self$penalty.fun == "entropy") {
        private$ent_optimizer()
      } else {
        stop("condition not considered in COT$optimize_weights. Please report this bug!")
      }
      
    } else {
      stop("condition not considered in COT$optimize_weights. Please report this bug!")
    }
    return(invisible(self))
  }
)

COT$set("private", "NNM_opt", 
  function() {
    
    C_xy <- self$ot$C_xy
    if (!self$ot$tensorized) { 
      x = as.matrix(C_xy$data$x)
      y = as.matrix(C_xy$data$y)
      d = ncol(x)
      
      # browser()
      use_cuda <- torch::cuda_is_available() && torch::cuda_device_count()>1
      rkeops::compile4float64()
      if (use_cuda) {
        rkeops::compile4gpu()
        rkeops::use_gpu()
      }
      argmin_op <- rkeops::keops_kernel(
        formula = paste0("ArgMin_Reduction(", C_xy$fun, ", 1)"),
        args = c(
          paste0("X = Vi(",d,")"),
          paste0("Y = Vj(",d,")"))
      )
      mins = torch::torch_tensor(c(argmin_op(list(x,y))) + 1, 
                                 dtype = torch::torch_int64())
        
    } else {
      mins = C_xy$data$argmin(1)
    }
    w_nnm = torch::torch_bincount(self = mins, weights = self$ot$b, minlength = self$ot$n)
    self$weight <- w_nnm
  }        
)

COT$set("private", "torch_optimizer_loop",
function() {
  loss_0 <- loss <- self$ot$diameter * 1000.
  tol <- self$tol
  
  #closure function for lbfgs
  if ( inherits(self$optim, "optim_lbfgs") ) {
    closure <- function() {
      # browser()
      self$optim$zero_grad()
      loss <- private$forward()
      loss$backward()
      return(loss)
    }
    for (i in 1:self$niter) {
      # browser()
      # self$optim$zero_grad()
      # loss <- private$forward()
      # loss$backward()
      loss <- self$optim$step(closure)
      if (converged(loss$item(), loss_0, tol)) break
      self$scheduler$step(loss)
      loss_0 <- loss$item()
    }
  } else {
    for (i in 1:self$niter) {
      self$optim$zero_grad()
      loss <- private$forward()
      loss$backward()
      self$optim$step()
      if (converged(loss$item(), loss_0, tol)) break
      self$scheduler$step(loss)
      loss_0 <- loss$item()
    }
  }
}
)

COT$set("private", "L2_optimizer",
  function() {
    # browser()
    P_new <- rep(self$ot$penalty, self$ot$n * self$ot$m)
    n <- self$ot$n
    m <- self$ot$m
    
    private$solver$Update(Px = P_new)
    Pi <- osqp_R6_solve(private$solver, 
                        self$bf$delta, 
                        private$delta.idx, 
                        rep(1/(n*m),n * m))
    self$weight <- rowSums(matrix(Pi, self$ot$n, self$ot$m))
  }
)

COT$set("private", "ent_optimizer",
  function() {
    runbf <- self$runbf
    n   <- self$ot$n
    m   <- self$ot$m
    lambda <- as.double(self$ot$penalty)
    C_xy <- as.matrix(self$ot$C_xy$data)
    
    if (runbf) {
      source <- self$bf$source
      target <- self$bf$target_vector
      delta <- as.double(self$bf$delta)
      k <- ncol(source)
      bounds <- rbind(cbind(rep(-Inf,m ), rep(Inf, m)),
      cbind(rep(0, 2*k), rep(Inf, 2*k)))
    } else {
      bounds <- cbind(rep(-Inf,m ), rep(Inf, m ))
      source <- matrix(0,0,0)
      target <- as.double(0.)
      delta <- as.double(0.)
      k <- 0
    }
    
    init <- c(rep(0, m),  rep(lambda, 2 * k))
    # browser()
    
    par <- lbfgs3c_R6_solve(init = init, options = private$solver.options,
                                    bounds = bounds,
                                    objective = cotEntropy_obj_,
                                    gradient = cotEntropy_grad_,
                                    source_ = source,
                                    target_ = target,
                                    cost_ = C_xy,
                                    b_ = as.double(self$ot$b),
                                    delta = delta,
                                    lambda = lambda
                                    )
    # get OT matrix
    g   <- par[1:m]
    eta <- -scale(C_xy/ lambda, center = g/lambda, scale = FALSE)
    if ( runbf ) { # if there are balance functions
      k <-  ncol(source)
      beta_unconst <- par[(m + 1):(m + 2 * k)]
      beta <- beta_unconst[1:k] - beta_unconst[-c(1:k)]
      eta <- eta - c(source %*% beta)/lambda
    }

    #set weights
    self$weight <- rowSums(exp(eta - logSumExp(c(eta))))
  }
)

COT$set("public", "gridInit",
function(grid, length) {

  runbf <- self$runbf
  delta <- grid$delta
  lambda <- grid$lambda
  if (is.numeric(lambda)) {
    lambda <- sort(lambda, decreasing = TRUE)
    if(!all(lambda >= 0)) stop("Penalty terms for COT must be numbers >= 0 or Inf.")
    delta <- sort(grid$delta, decreasing = TRUE)
    if(runbf && !all(delta >= 0)) stop("Balancing function constraints for COT must be numbers >= 0.")

  } else if (missing(grid) || is.na(lambda) || is.null(lambda)) {
    diam <- self$ot$diameter
    lambda <- c(Inf,
                exp(seq(log(1e2) + log(diam), log(1e-4) + log(diam), length.out = length)),
                0)

    if(runbf) {
      delta <- self$bf$gridInit(NULL, length)
    }
  } else {
    stop("Grid values weren't provided but I can't initialize any. Report this bug!")
  }

  if(runbf) { # get approximately good value for delta
    w <- vector("list", length(delta))
    for(d in seq_along(delta)) {
      w[[d]] <- self$bf$solve(delta[d], w[[d-1]])
    }
    metric <- bootStrap_(w,
                         as.integer(100),
                         self$bf)
    delta_sel <- delta[which.min(metric)]

    penalty_list <- lapply(lambda, function(l) c(lambda = l, delta = delta_sel))
  } else {
    penalty_list <- lambda
  }

  return(penalty_list)
}
)

#forward functions
COT$set("private", "energy_dist_forward",
function() {
  self$ot$a <- self$weight
  return(energy_dist(self$ot))
})

COT$set("private", "cot_forward",
function() {
  self$ot$a <- self$weight
  # browser()
  self$ot$sinkhorn_cot(niter = self$niter, tol = self$tol)
  loss <- sinkhorn_dist(self$ot)

  return(loss)
})

# COT options function
#' Options available for the COT method
#'
#' @param lambda The penalty parameter for the entropy penalized optimal transport. Default is NULL. Can be a single number or a set of numbers to try.
#' @param delta The bound for balancing functions if they are being used. Only available for biased entropy penalized optimal transport. Can be a single number or a set of numbers to try.
#' @param penalty.function The penalty to use for the optimal transport distances. Should be one of "entropy" or "L2".
#' @param debias Should debiased optimal transport be used? TRUE or FALSE. Only available for `penalty.function`="entropy". 
#' @param p The power of the cost function to use for the cost.
#' @param cost A function to calculate the pairwise costs. Should take arguments `x1`, `x2`, and `p`. Default is NULL.
#' @param cost.online Should an online cost algorithm be used? One of "auto", "online", or "tensorized". "tensorized" is the offline option.
#' @param balance.formula Formula for the balancing functions.
#' @param grid.length The number of penalty parameters to explore in a grid search if none are provided in arguments `lambda` or `delta`.
#' @param torch.optimizer The torch optimizer to use for methods using debiased entropy penalized optimal transport.
#' @param torch.scheduler The scheduler for the optimizer
#' @param niter The number of iterations to run the solver
#' @param nboot The number of iterations for the bootstrap to select final penalty parameters.
#' @param tol The tolerance for convergence
#' @param ... Arguments passed to the solvers. See details
#'
#' @return A list of class `cotOptions` with the following slots
#' \itemize{
#' \item `lambda`The penalty parameter for the optimal transport distance
#' \item `delta`The constraint for the balancing functions
#' \item `penalty.function`The penalty function to use
#' \item `debias`TRUE or FALSE if debiased optimal transport distances are used
#' \item `balance.formula`
#' \item `grid.length` The number of parameters to check in a grid search of best parameters
#' \item `p` The power of the cost function
#' \item `cost.online` Whether online costs are used
#' \item `cost` The user supplied cost function if supplied.
#' \item `torch.optimizer` The `torch` optimizer used for Sinkhorn Divergences
#' \item `torch.scheduler` The scheduler for the `torch` optimizer
#' \item `solver.options` The arguments to be passeed to the `torch.optimizer`
#' \item `scheduler.options` The arguments to be passeed to the `torch.scheduler`
#' \item `niter` The number of iterations to run the solver
#' \item `nboot` The number of bootstrap samples
#' \item `tol` The tolerance for convergence.
#' }
#' @export
#' 
#' @details 
#' # Solvers and distances
#' The function is setup to direct the COT optimizer to run three basic methods: debiased entropy penalized optimal transport (Sinkhorn Divergences), entropy penalized optimal transport (Sinkhorn Distances), or L2 penalized optimal transport. Each of these options utilizes a specific solver.
#' 
#' ## Sinkhorn Distances
#' The optimal transport problem solved is \eqn{min_w OT_\lambda(w,b) } where \deqn{OT_\lambda(w,b) = \sum_{ij} C(x_i, x_j) P_{ij} + \lambda \sum_{ij} P_{ij}\log(P_{ij}),} such that the rows of the matrix \eqn{P_{ij}} sum to \eqn{w} and the columns sum to \eqn{b}. In this case \eqn{C(,)} is the cost between units i and j. 
#' 
#' The solver for this distance is provided by [lbfgsb3c()][lbfgsb3c::lbfgsb3c()] and so arguments in dots are passed to the function. We can also supply balancing functions via the formula argument in `balance.formula` (see below).
#' 
#' ## L2 OT Distances
#' In this case, the optimal transport problem is again \eqn{min_w OT_\lambda(w,b) } where \deqn{OT_\lambda(w,b) = \sum_{ij} C(x_i, x_j) P_{ij} + \frac{\lambda}{2} \sum_{ij} P_{ij}^2,}
#' 
#' The solver for this distance is provided by [osqp()][osqp::osqp()] and so arguments in dots are passed to the solver via that packages [osqpSettings()][osqp::osqpSettings()] function. We can also supply balancing functions via the formula argument in `balance.formula` (see below).
#' 
#' ## Sinkhorn Divergences
#' The Sinkhorn Divergence solves \deqn{min_w OT_\lambda(w,b) - 0.5 OT_\lambda(w,w) - 0.5 * OT_\lambda(b,b).} The solver for this function uses the `torch` package in `R` and by default will use the LBFGS solver. Your desired `torch` optimizer can be passed via `torch.optimizer` with a scheduler passed via `torch.scheduler`. GPU support is available as detailed in the `torch` package. Additional arguments in `...` are passed as extra arguments to the `torch` optimizer and schedulers as appropriate.
#' 
#' # Function balancing
#' There may be certain functions of the covariates that we wish to balance within some tolerance, \eqn{\delta}. For these functions \eqn{B}, we will desire
#' \deqn{\frac{\sum_{i: Z_i = 0} w_i B(x_i) - \sum_{j: Z_j = 1} B(x_j)/n_1}{\sigma} \leq \delta}, where in this case we are targeting balance with the treatment group for the ATT. $\sigma$ is the pooled standard deviation prior to balancing. This is currently only available for \eqn{L_2} penalized distances and Sinkhorn Distances.
#' 
#' # Cost functions
#' The cost function specifies pairwise distances. If argument `cost` is NULL, the function will default to using \eqn{L_p^p} distances with a default \eqn{p = 2} supplied by the argument `p`.So for `p = 2`, the cost between units \eqn{x_i} and \eqn{x_j} will be \deqn{C(x_i, x_j) = \frac{1}{2} \| x_i - x_j \|_2^2.}
#' If `cost` is provided, it should be a function that takes arguments `x1`, `x2`, and `p`.
#' 
#'
#' @examples
#' opts <- cotOptions(lambda = 1e3, torch.optimizer = torch::optim_lbfgs)
#' opts <- cotOptions(lambda = NULL)
#' opts <- cotOptions(lambda = seq(0.1, 100, length.out = 7))
cotOptions <- function(lambda = NULL,
                       delta = NULL,
                       penalty.function = "entropy",
                       debias = TRUE,
                       p = 2.0,
                       cost = NULL,
                       cost.online = "auto",
                       balance.formula = NULL,
                       grid.length = 7L,
                       torch.optimizer = torch::optim_lbfgs,
                       torch.scheduler = torch::lr_reduce_on_plateau,
                       niter = 1e3,
                       nboot = 100L,
                       tol = 1e-7, ...) { # dots are the solver and scheduler args
  mc <- match.call()
  used.args <- as.list(mc)[-1]
  # browser()
  gsOpts <- gridSearchOptions(nboot = nboot, grid.length = grid.length)
  
  nboot <- gsOpts$nboot
  grid.length <- gsOpts$grid.length
  
  output <- list()
  if(arg_not_used(lambda)) {
    output["lambda"] <- list(NULL)
  } else {
    if(any(lambda < 0)) stop("lambda must be >= 0")
    output$lambda <- sort(lambda, decreasing = TRUE)
  }

  if(arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  if ( arg_not_used(penalty.function) ) {
    output$penalty.function <- "entropy"
  } else {
    output$penalty.function <- match.arg( penalty.function, c("entropy", "L2") )
  }
  
  if ( arg_not_used(debias) ) {
    if(output$penalty.function == "entropy") {
      output$debias  <- TRUE
    } else if (output$penalty.function == "L2") {
      output$debias <- FALSE
    } else {
      stop("cotOption error: penalty function not found!")
    }
  } else {
    output$debias  <- isTRUE( debias )
    if (is.null(used.args$debias) && output$penalty.function == "L2") output$debias <- FALSE
  }
  if( arg_not_used(balance.formula) ) {
    output["balance.formula"] <- list(NULL)
  } else {
    balance.formula <- as.character(balance.formula)
    bf_split <- strsplit(balance.formula, "~")
    output$balance.formula <- paste0("~ 0 +", bf_split[[1]][2])
  }
  if ( arg_not_used(grid.length) ) {
    output$grid.length <- 7L
  } else {
    output$grid.length <- as.integer(grid.length)
    if(grid.length <= 0) stop("grid.length must be greater than 0")
  }
  if (!is.null(output$lambda) && !is.null(output$delta)) {
    output["grid.length"] <- list(NULL)
  }
  
  if ( arg_not_used(p) ) {
    output[["p"]] <- 2.0
  } else {
    p <- as.numeric(p)
    stopifnot(p >= 1.0)
    output[["p"]] <- p
  }
  
  if ( arg_not_used(cost.online) ) {
    output$cost.online <- "auto"
  } else {
    output$cost.online <- match.arg(cost.online,
                                             c("auto",
                                               "tensorized",
                                               "online"))
  }
  if ( arg_not_used(cost) ) {
    output["cost"] <- list(NULL)
  } else {
    output["cost"]  <- cost
  }
  if ( (missing(torch.optimizer) || is.null(torch.optimizer)) && !output$debias) {
    output["torch.optimizer"] <- list(NULL)
    output["solver.options"] <- list(NULL)
  } else {
    if (!inherits(torch.optimizer, "torch_optimizer_generator")) {
      stop("torch.optimizer must be a torch_optimizer_generator function or NULL")
    }
    output$torch.optimizer <- torch.optimizer
    output$solver.options <- list(...)[...names() %in% formalArgs(torch.optimizer)]
  }
  if (missing(torch.scheduler) || is.null(torch.scheduler) || is.null(torch.optimizer)) {
    output["torch.scheduler"] <- list(NULL)
    output["scheduler.options"] <- list(NULL)
  } else {
    if (!inherits(torch.scheduler, "lr_scheduler")) {
      stop("torch.scheduler must be a lr_scheduler function or NULL")
    }
    output$torch.scheduler <- torch.scheduler
    output$scheduler.options <- list(...)[...names() %in% formalArgs(torch.scheduler)]
  }

  if (arg_not_used (niter)) {
    output$niter <- 1000L
  } else {
    output$niter <- as.integer(niter)
  }

  if (arg_not_used (nboot)) {
    output$nboot <- 100L
  } else {
    output$nboot <- as.integer(nboot)
  }

  if (arg_not_used(tol)) {
    output$tol <- 1e-7
  } else {
    output$tol <- as.double(tol)
  }
  
  # check combinations
  if(output$debias && is.null(output$torch.optimizer)) {
    stop("Must supply a 'torch_optimizer_generator' object in options argument 'torch.optimizer' when option debias is TRUE")
  }
  
  if(output$debias && output$penalty.function == "L2") {
    warning("No debias options with L2 penalty. Defaulting to debias=FALSE")
    output$debias <- FALSE
  }
  
  if (output$penalty.function == "L2") {
    output$solver.options <- list(...)[...names() %in% formalArgs(osqp::osqpSettings)]
  }
  if (!output$debias && output$penalty.function == "entropy" && output$cost.online != "online") {
    output$solver.options <- lbfgs3c_control(...)
  }

  class(output) <- "cotOptions"
  return(output)
}

COT$set("private", "L2_setup", # TODO how to setup this function... don't need to add online support all at once
 function(options) {
   # browser()
   solver_args <- options$solver.options
   osqp_options <- do.call(osqp::osqpSettings, solver_args)
   n <- self$ot$n
   m <- self$ot$m
   nvar <- n * m
   P <- Matrix::Diagonal(nvar, x = self$ot$penalty)
   q <- c(as.numeric(self$ot$C_xy$data))

   A_bounds <- Matrix::Diagonal(nvar, 1)
   l_bounds <- rep(0, nvar)
   u_bounds <- rep(Inf, nvar)

   A_b   <- vec_to_col_constraints(n,m)
   l_b   <- u_b <- as.numeric(self$ot$b)

   A_sum <- Matrix::sparseMatrix(i = rep(1,nvar), j = 1:nvar, x = 1)
   l_sum <- u_sum <- 1

   if (self$runbf) {
     A_bf <- Matrix::crossprod(self$bf$A, 
                               vec_to_row_constraints(n,m))
     u_bf <- rep(self$bf$delta, nrow(A_bf))
     l_bf <- -u_bf
     private$delta.idx <- 1:length(l_bf)
     
   } else {
     A_bf <- l_bf <- u_bf <- NULL
   }

   A <- rbind(A_bf, A_bounds, A_sum, A_b)
   l <- c(l_bf, l_bounds, l_sum, l_b)
   u <- c(u_bf, u_bounds, u_sum, u_b)

   
   private$solver <- osqp::osqp(P = P, q = q, A = A, l = l, u = u, pars = osqp_options)
   # self$solve <- function(penalty, w) {
   # 
   #   if (length(w) > 0) {
   #     w <- w[[1]]
   #     n_param <- length(self$param)
   #     if( length(w) == n_param ) {
   #       self$weight <- w
   #       model$WarmStart(x = w)
   #     } else {
   #       warning("Warm start not used. Proposed starting parameter of different length than params")
   #     }
   #   }
   # 
   #   if (self$runbf) {
   #     lambda <- as.numeric(penalty[["lambda"]])
   #     delta  <- as.numeric(penalty[["delta"]])
   #   } else {
   #     lambda <- as.numeric(penalty)
   #     delta  <- NULL
   #   }
   #   private$solver$Update
   #   self$weight <- osqp_R6_solve(private$solver, delta, 
   #                                private$delta.idx, self$weight)
   #   return(self$weight)
   # }
 }

)
# COT$set("private", "entropy_setup", # TODO how to setup this function... don't need to add online support all at once
#   function(options) {
#     
#     browser()
#     private$solver <- list(objective = causalOT:::cotEntropy_obj_,
#                            gradient = causalOT:::cotEntropy_grad_)
#     
#     
#   }
# 
# )

#### SCM ####
SCM <- R6::R6Class("SCM",
 public = list(
   source = "matrix",
   target = "matrix",
   a = "vector",
   b = "vector",
   A = "matrix",
   bf = "R6",
   gridInit = function(...) {
     return(c(NA_real_))
   },
   solve = function(penalty = NULL, w = NULL) {

     n_t <- nrow(self$target)
     n_s <- nrow(self$source)
     d   <- ncol(self$source)
     
     w_list <- vector("list", n_t)
     res <- vector("list", 1)
     addl_0 <- rep(0,d)
     w_init <- c(rep(1/n_s,n_s), addl_0)
    
     # update BF penalty function
     if(inherits(self$bf, "balanceFunction") && !is.null(penalty)) {
       if(is.list(penalty)) penalty <- penalty[[1]]
       if(is.numeric(penalty) && penalty >= 0) {
         self$bf$delta <- as.numeric(penalty)
       } else {
         warning("provided penalty for balanceFunction method is not a number > 0. Not used!")
       }
     }
     
     #update target
     t_idx <- private$target_idx
     l <- private$solver$GetData(element = "l")
     u <- private$solver$GetData(element = "u")

     for (j in 1:n_t) {
       l[t_idx] <- u[t_idx] <- c(self$target[j,])
       private$solver$Update(l = l, u = u)
       res[[1]] <- osqp_R6_solve(private$solver, self$bf$delta, 
                            self$bf$delta_idx, w_init, normalize = FALSE)
       res[[1]][res[[1]] < 0] <- 0
       w_list[[j]] <- renormalize(res[[1]][1:n_s]) * self$b[j]
     }
     return(Reduce("+", w_list))
   },
   initialize = function(source, target,
                         a, b,
                         method, options) {
     # browser()
     if(!inherits(options, "scmOptions")) {
       if(!is.list(options)) stop("options must be a list or output of scmOptions function")
       options <- do.call(scmOptions, options)
     }
     
     self$source <- as.matrix(source)
     self$target <- as.matrix(target)
     self$a <- private$check_weights(a, self$source)
     self$b <- private$check_weights(b, self$target)

     n <- nrow(self$source)
     d <- ncol(self$source)
     if (!is.null(options$balance.formula)) {
       if (is.null(colnames(source))) colnames(source) <- paste0("X", 1:ncol(source))
       if (is.null(colnames(target))) colnames(target) <- colnames(source)
       tbf <- terms(formula(options$balance.formula))
       attr(tbf, "intercept") <- 0
       source.bf <- model.matrix(tbf, data.frame(source))
       target.bf <- model.matrix(tbf, data.frame(target))
       self$bf <- balanceFunction$new(source = source.bf, 
                                      target = target.bf,
                                      a = as.numeric(self$ot$a),
                                      b = as.numeric(self$ot$b),
                                      delta = options$delta)
       # self$runbf <- TRUE
       k <- ncol(target.bf)
     } else {
       # self$runbf <- FALSE
       self$bf <- list(NULL)
       k <- 0
     }

     nvars <- n + d

     # q <- c(self$source %*% c(self$target[1,]), rep(0, d))
     q <- NULL
     P <- Matrix::sparseMatrix(i = n + 1:d, j = n + 1:d, x = 1)

     A_quad <- cbind(Matrix::Matrix(t(self$source)),
                     Matrix::Diagonal(d, x = 1)
                )
     u_quad <- l_quad <- c(self$target[1,])

     A_bounds <- rbind(cbind(Matrix::Diagonal(n, x = 1),
                       Matrix::Matrix(matrix(0, n, d))
     ),
     Matrix::sparseMatrix(i = rep(1,n), j = 1:n, x = 1, dims = c(1, nvars)))
     l_bounds <- c(rep(0,   n), 1)
     u_bounds <- c(rep(Inf, n), 1)

     if(k > 0) {
       A_bf <- cbind(
         Matrix::Matrix(
           t(self$bf$A)
         ),
         Matrix::Matrix(matrix(0, k, d)))
       l_bf <- rep(-self$bf$delta, k)
       u_bf <- rep(self$bf$delta, k)
       self$bf$delta_idx <- 1:length(l_bf)
     } else {
       A_bf <- l_bf <- u_bf <- NULL
     }

     l <- c(l_bf, l_bounds, l_quad)
     u <- c(u_bf, u_bounds, u_quad)
     A <- rbind(A_bf, A_bounds, A_quad)
     
     private$target_idx <- (k + n + 2):length(l)

     private$solver <- osqp::osqp(P = P, q = q,
                               A = A, l = l, u = u,
                               pars = options$solver.options)

   }
 ),
 private = list(
   check_weights = function(a, x) {
     if(missing(a) || is.null(a) || all(is.na(a))) {
       a <- rep(1.0/nrow(x), nrow(x))
     } else {
       a <- renormalize(a)
     }
     return(a)
   },
   target_idx = "numeric",
   solver = "R6"
 )
)


# @param lambda Penalty parameter on the weights. Not currently used but here because a plan is to add it.
# @param delta The constraint parameter for the balancing functions. Not currently used.
# @param grid.length The number of penalty parameters to try. Not currently used.
# @param nboot The number of bootstrap samples. Not currently used.
# @param balance.formula The formula that denotes the covariate functions to balance. Not currently used.

#' Options for the SCM Method
#'
#' @param ... Arguments passed to the [osqpSettings()][osqp::osqpSettings()] function which solves the problem.
#'
#' @return A list with arguments to pass to [osqpSettings()][osqp::osqpSettings()]
#' @export
#' 
#' @details Options for the solver used in the optimization of the Synthetic Control Method of Abadie and Gardeazabal (2003).
#'
#' @examples
#' opts <- scmOptions()
scmOptions <- function(
    # lambda = NULL,
    #                    delta = NULL,
    #                    grid.length = 7L,
    #                    nboot = 1000L,
                       # balance.formula = NULL,
                       ...) { # dots are the osqp args
  mc <- match.call()
  used.args <- as.list(mc)[-1]
  # browser()
  
  grid.length <- 7L
  nboot <- 1000L
  delta <- NULL
  balance.formula <- NULL
  gsOpts <- gridSearchOptions(nboot = nboot, grid.length = grid.length)
  
  nboot <- gsOpts$nboot
  grid.length <- gsOpts$grid.length
  
  output <- list()
  # lambda currently not used but may consider in future with following code
  # if(arg_not_used(lambda)) {
  #   output["lambda"] <- list(NULL)
  # } else {
  #   if(any(lambda < 0)) stop("lambda must be >= 0")
  #   output$lambda <- sort(lambda, decreasing = TRUE)
  # }
  output$lambda <- NULL
  
  if(arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  
  # also only L2 penalty at this time but may consider in future
  # if ( arg_not_used(penalty) ) {
  #   output$penalty <- "entropy"
  # } else {
  #   output$penalty <- match.arg( penalty, c("entropy", "L2") )
  # }
  
  if( arg_not_used(balance.formula) ) {
    output["balance.formula"] <- list(NULL)
  } else {
    balance.formula <- as.character(balance.formula)
    bf_split <- strsplit(balance.formula, "~")
    output$balance.formula <- paste0("~ 0 +", bf_split[[1]][2])
  }
  
  # only for delta at this time
  if ( arg_not_used(grid.length) ) {
    output$grid.length <- 7L
  } else {
    output$grid.length <- as.integer(grid.length)
    if(grid.length <= 0) stop("grid.length must be greater than 0")
  }
  
  if (!is.null(nboot)) {
    output$nboot <- as.integer(nboot)
  } else {
    output$nboot <- 1000L
  }
  
  # if (!is.null(output$lambda) && !is.null(output$delta)) {
  #   output["grid.length"] <- list(NULL)
  # }
  if ( !is.null(output$delta)) {
    output["grid.length"] <- list(NULL)
  }
  
  # if (output$penalty == "L2") {
    output$solver.options <- list(...)[...names() %in% formalArgs(osqp::osqpSettings)]
  # }
  
  class(output) <- "scmOptions"
  return(output)
}

#### general functions ####

balanceDistributions <- function(source, target,
                                 a, b,
                                 method, options) {
  if(method == "SCM") {
    if (!is.list(options)) options <- list(options)
    
    if(!inherits(options, "scmOptions")) options <- do.call(scmOptions, options)
    
    prob <- SCM$new(source = source, target = target,
                    a = a, b = b,
                    options = options)
  } else if ( method %in% cot_methods() ) {
    if (!is.list(options)) options <- list(options)
    
    if(!inherits(options, "cotOptions")) options <- do.call(cotOptions, options)
    options["lambda"] <- list(switch(method,
                             "EnergyBW" = Inf,
                             "NNM" = 0,
                             options$lambda))
    
    prob <- COT$new(source = source, target = target,
                    a = a, b= b,
                    options = options)
  } else {
    stop("Method not found in function balanceDistributions.")
  }
    
  return(list(problem = prob,
              options = options))

}
# 
# setMethod("grid_solve", signature(object = "balanceDistributions",
#                                   penalty = "numeric",
#                                   w = "numeric"),
# function(object, penalty, w) {
# 
# }
# )
