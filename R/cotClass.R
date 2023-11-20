
# setOldClass(c("gridSearchClass","R6"))
# setOldClass("balanceFunction")
# setOldClass(c("OT","R6"))
# setOldClass("torch_tensor")

#' @include balanceFunctions.R
#' @include gridSearch.R
#' @include OT.R

#### NEW COT CLASS ####
#' COT base class
COT <- R6::R6Class(
  classname = "COT",
  public = {list(
    niter = "integer",
    tol = "numeric",
    lambda_bootstrap = "numeric",
    nboot = "integer",
    gridInit = function(...) {
      private$optimizer$penalty$lambda
    },
    grid_search = function(...) {
      
      private$optimizer$choose_hyperparameters(n_boot_lambda = self$nboot, 
                                               n_boot_delta = self$nboot,
                                               lambda_bootstrap = self$lambda_bootstrap)
      
      lambda <- private$optimizer$selected_lambda
      if (length( private$optimizer$selected_delta ) == 1 && is.numeric(private$optimizer$selected_delta) ) {
        delta <- private$optimizer$selected_delta
      } else if (length( private$optimizer$selected_delta ) == 1 && length(private$optimizer$selected_delta[[1]]) == 1 && is.list(private$optimizer$selected_delta)) {
        delta <- private$optimizer$selected_delta[[1]]
      } else {
        delta <- NULL
      }
      # private$weights_metrics
      # browser()
      return(
        list(
          weight = self$weights,
          penalty = c(lambda = lambda,
                      delta  = delta),
          metric  = private$optimizer$info()$hyperparam.metrics,
          penalty.grid = private$optimizer$penalty
        )
      )
    },
    solve = function(...) {
      # niter = 100L, tol = 1e-5, optimizer = c("torch", "frank-wolfe"),
      # torch_optim = torch::optim_lbfgs,
      # torch_scheduler = torch::lr_reduce_on_plateau,
      # torch_args = NULL,
      # osqp_args = NULL,
      # quick.balance.function = TRUE
      private$optimizer$solve(niter = self$niter, 
                              tol = self$tol,
                              torch_optim = private$torch_optim,
                              torch_scheduler = private$torch_sched,
                              torch_args = private$torch_optim_args,
                              osqp_args = private$osqp_args,
                              quick.balance.function = private$quick.balance.function)
      
    },
    initialize = function(source, target,
                          a = NULL, b = NULL,
                          options = list(NULL)) {
      # browser()
      if (missing(source) || missing(target)) stop("Must provide source and target")
      if (!inherits(options, "cotOptions")) options <- do.call(cotOptions, options)
      
      self$niter = options$niter
      self$tol = options$tol
      self$nboot = options$nboot
      self$lambda_bootstrap = options$lambda.bootstrap
      
      #setup bf
      if (!is.null(options$balance.formula)) {
        if (is.null(colnames(source))) colnames(source) <- paste0("X", 1:ncol(source))
        if (is.null(colnames(target))) colnames(target) <- colnames(source)
        dfs <- data.frame(source)
        tbf <- terms(formula(options$balance.formula), data = dfs)
        attr(tbf, "intercept") <- 0
        source.bf <- model.matrix(tbf, dfs)
        target.bf <- model.matrix(tbf, data.frame(target))
        
        target.sd <- matrixStats::colSds(target.bf)
        target.bf <- scale(target.bf, center = FALSE, scale = target.sd)
        
        balance.functions <- scale(source.bf, center = FALSE, scale = target.sd)
        target.values <- colMeans(target.bf)
        
      } else {
        balance.functions <- NULL
        target.values <- NA_real_
      }
      
      private$source <- Measure(x = source, weights = a, adapt = "weights",
                                balance.functions = balance.functions,
                                target.values = target.values,
                                dtype = options$dtype,
                                device = options$device)
      private$target <- Measure(x = target, weights = b, adapt = "none",
                                dtype = options$dtype,
                                device = options$device)
      
      
      # select dual or primal optimization
      opt <- options$opt.direction
      stopifnot("Option `opt.direction` should be one of `dual` or `primal`." = opt %in% c("primal","dual"))
      
      runbf <- !(all(is.na(private$source$balance_target)))
      n_lambda <- length(options$lambda)
      
      which.opt.fun <- if (n_lambda == 1 && options$lambda == 0 && !runbf ) {
        "nnm"
      } else if (n_lambda == 1 && options$lambda == 0 && runbf ) {
        "dual"
      } else if (n_lambda == 1 && is.infinite(options$lambda) ) {
        "primal"
      } else if ( opt == "dual" ) {
        "dual"
      } else if ( opt == "primal" ) {
        "primal"
      } else {
        stop("You found a bug! Please report this issue")
      }
      
      opt_fun <- switch(which.opt.fun,
                        "primal" = OTProblem_,
                        "dual" = cotDualTrain,
                        "nnm" = NNM)
      
      private$optimizer <- opt_fun$new(private$source, private$target) 
      # the order is important, first measure must be the one to be adapted
      
      # 
      private$quick.balance.function <- isTRUE(options$quick.balance.function)
      
      
      
      private$torch_sched <- options$torch.scheduler
      private$torch_optim <- options$torch.optimizer
      private$torch_optim_args <- c(options$solver.options, options$scheduler.options)
      private$osqp_args <- options$osqp.options
      
      names_args <- names(options)
      optim_args_names <- names(formals(private$optimizer$setup_arguments))
      opt_args <- options[match(optim_args_names,
                                names_args, nomatch = 0L)]
      do.call(private$optimizer$setup_arguments, opt_args)
      
      return(invisible(self))
    }
  )},
  active = {list(
    b = function(value) {
      as.numeric(private$target$weights$to(device = "cpu"))
    },
    weights = function(value) {
      if(missing(value)) {
        as.numeric(private$source$weights$to(device = "cpu"))
      } else {
        private$source$weights <- value
      }
    },
    penalty = function(value) {
      if(missing(value)) {
        return(private$optimizer$penalty)
      } else {
        stopifnot(is.list(value))
        stopifnot(all(names(value) %in% c("lambda", "delta")))
        stopifnot(all(c("lambda", "delta") %in% names(value)))
        private$optimizer$penalty <- value
      }
      
    }
  )},
  private = {list(
    optimizer = "OTProblem",
    osqp_args = "list",
    quick.balance.function = "logical",
    source = "Measure",
    target = "Measure",
    torch_optim = "torch_optim",
    torch_sched = "torch_sched",
    torch_optim_args = "list"
  )}
)

#### options function to make things hopefully easier for users ####

# COT options function
#' Options available for the COT method
#'
#' @param lambda The penalty parameter for the entropy penalized optimal transport. Default is NULL. Can be a single number or a set of numbers to try.
#' @param delta The bound for balancing functions if they are being used. Only available for biased entropy penalized optimal transport. Can be a single number or a set of numbers to try.
#' @param opt.direction Should the optimizer solve the primal or dual problems. Should be one of "dual" or "primal" with a default of "dual" since it is typically faster.
#' @param debias Should debiased optimal transport be used? TRUE or FALSE.
#' @param p The power of the cost function to use for the cost.
#' @param cost.function A function to calculate the pairwise costs. Should take arguments `x1`, `x2`, and `p`. Default is NULL.
#' @param cost.online Should an online cost algorithm be used? One of "auto", "online", or "tensorized". "tensorized" is the offline option.
#' @param diameter The diameter of the covariate space, if known. Default is NULL.
#' @param balance.formula Formula for the balancing functions.
#' @param quick.balance.function TRUE or FALSE denoting whether balance function constraints should be selected via a linear program (TRUE) or just checked for feasibility (FALSE). Default is TRUE. 
#' @param grid.length The number of penalty parameters to explore in a grid search if none are provided in arguments `lambda` or `delta`.
#' @param torch.optimizer The torch optimizer to use for methods using debiased entropy penalized optimal transport. If `debiased` is FALSE or `opt.direction` is "primal", will default to [torch::optim_lbfgs()]. Otherwise [torch::optim_rmsprop()] is used.
#' @param torch.scheduler The scheduler for the optimizer. Defaults to [torch::lr_multiplicative()].
#' @param niter The number of iterations to run the solver
#' @param nboot The number of iterations for the bootstrap to select the final penalty parameters.
#' @param lambda.bootstrap The penalty parameter to use for the bootstrap hyperparameter selection of lambda.
#' @param tol The tolerance for convergence
#' @param device An object of class `torch_device` denoting which device the data will be located on. Default is NULL which will try to use a gpu if available.
#' @param dtype An object of class `torch_dtype` that determines data type of the data, i.e. double, float, integer. Default is NULL which will try to select for you.
#' @param ... Arguments passed to the solvers. See details
#'
#' @return A list of class `cotOptions` with the following slots
#' \itemize{
#' \item `lambda`The penalty parameter for the optimal transport distance
#' \item `delta`The constraint for the balancing functions
#' \item `opt.direction` Whether to solve the primal or dual optimization problems
#' \item `debias`TRUE or FALSE if debiased optimal transport distances are used
#' \item `balance.formula` The formula giving how to generate the balancing functions.
#' \item `quick.balance.function` TRUE or FALSE whether quick balance functions will be run.
#' \item `grid.length` The number of parameters to check in a grid search of best parameters
#' \item `p` The power of the cost function
#' \item `cost.online` Whether online costs are used
#' \item `cost.function` The user supplied cost function if supplied.
#' \item `diameter` The diameter of the covariate space.
#' \item `torch.optimizer` The `torch` optimizer used for Sinkhorn Divergences
#' \item `torch.scheduler` The scheduler for the `torch` optimizer
#' \item `solver.options` The arguments to be passeed to the `torch.optimizer`
#' \item `scheduler.options` The arguments to be passeed to the `torch.scheduler`
#' \item `osqp.options` Arguments passed to the `osqp` function if quick balance functions are used.
#' \item `niter` The number of iterations to run the solver
#' \item `nboot` The number of bootstrap samples
#' \item `lambda.bootstrap` The penalty parameter to use for the bootstrap hyperparameter selection.
#' \item `tol` The tolerance for convergence.
#' \item `device` An object of class `torch_device`.
#' \item `dtype` An object of class `torch_dtype`.
#' }
#' @export
#' 
#' @details 
#' # Solvers and distances
#' The function is setup to direct the COT optimizer to run two basic methods: debiased entropy penalized optimal transport (Sinkhorn Divergences) or entropy penalized optimal transport (Sinkhorn Distances). 
#' 
#' ## Sinkhorn Distances
#' The optimal transport problem solved is \eqn{min_w OT_\lambda(w,b) } where \deqn{OT_\lambda(w,b) = \sum_{ij} C(x_i, x_j) P_{ij} + \lambda \sum_{ij} P_{ij}\log(P_{ij}),} such that the rows of the matrix \eqn{P_{ij}} sum to \eqn{w} and the columns sum to \eqn{b}. In this case \eqn{C(,)} is the cost between units i and j. 
#' 
#' ## Sinkhorn Divergences
#' The Sinkhorn Divergence solves \deqn{min_w OT_\lambda(w,b) - 0.5 OT_\lambda(w,w) - 0.5 * OT_\lambda(b,b).} The solver for this function uses the `torch` package in `R` and by default will use the `optim_rmsprop` solver. Your desired `torch` optimizer can be passed via `torch.optimizer` with a scheduler passed via `torch.scheduler`. GPU support is available as detailed in the `torch` package. Additional arguments in `...` are passed as extra arguments to the `torch` optimizer and schedulers as appropriate.
#' 
#' # Function balancing
#' There may be certain functions of the covariates that we wish to balance within some tolerance, \eqn{\delta}. For these functions \eqn{B}, we will desire
#' \deqn{\frac{\sum_{i: Z_i = 0} w_i B(x_i) - \sum_{j: Z_j = 1} B(x_j)/n_1}{\sigma} \leq \delta}, where in this case we are targeting balance with the treatment group for the ATT. \eqn{\sigma} is the pooled standard deviation prior to balancing.
#' 
#' # Cost functions
#' The cost function specifies pairwise distances. If argument `cost.function` is NULL, the function will default to using \eqn{L_p^p} distances with a default \eqn{p = 2} supplied by the argument `p`. So for `p = 2`, the cost between units \eqn{x_i} and \eqn{x_j} will be \deqn{C(x_i, x_j) = \frac{1}{2} \| x_i - x_j \|_2^2.}
#' If `cost.function` is provided, it should be a function that takes arguments `x1`, `x2`, and `p`: `function(x1, x2, p){...}`.
#' 
#'
#' @examples
#' if ( torch::torch_is_installed()) {
#' opts1 <- cotOptions(lambda = 1e3, torch.optimizer = torch::optim_rmsprop)
#' opts2 <- cotOptions(lambda = NULL)
#' opts3 <- cotOptions(lambda = seq(0.1, 100, length.out = 7))
#' }
cotOptions <- function(lambda = NULL,
                       delta = NULL,
                       opt.direction = c("dual", "primal"),
                       debias = TRUE,
                       p = 2.0,
                       cost.function = NULL,
                       cost.online = "auto",
                       diameter = NULL,
                       balance.formula = NULL,
                       quick.balance.function = TRUE,
                       grid.length = 7L,
                       torch.optimizer = torch::optim_rmsprop,
                       torch.scheduler = torch::lr_multiplicative,
                       niter = 2e3,
                       nboot = 100L,
                       lambda.bootstrap = 0.05,
                       tol = 1e-4, 
                       device = NULL,
                       dtype = NULL,
                       ...) { # dots are the solver and scheduler args
  # setup_arguments function arguments
  # lambda, delta, 
  # grid.length = 7L,
  # cost.function = NULL, 
  # p = 2,
  # cost.online = "auto",
  # debias = TRUE,
  # diameter = NULL, niter = 1000L,
  # tol = 1e-3
  
  mc <- match.call()
  used.args <- as.list(mc)[-1]
  # browser()
  gsOpts <- gridSearchOptions(nboot = nboot, grid.length = grid.length)
  
  nboot <- gsOpts$nboot
  grid.length <- gsOpts$grid.length
  
  output <- list()
  if(missing(lambda) || arg_not_used(lambda)) {
    output["lambda"] <- list(NULL)
  } else {
    if(any(lambda < 0)) stop("lambda must be >= 0")
    output$lambda <- sort(lambda, decreasing = TRUE)
  }

  if(missing(delta) || arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  
  if(missing(opt.direction) || arg_not_used(opt.direction)) {
    output[["opt.direction"]] <- "dual"
  } else {
    output[["opt.direction"]] <- match.arg(opt.direction, c("dual", "primal"))
  }
  
  if (missing(debias) || arg_not_used(debias) ) {
      output$debias  <- TRUE
  } else {
    output$debias  <- isTRUE( debias )
  }
  if ((! output$debias) && output$opt.direction == "dual") output$opt.direction <- "primal"
  
  if(missing(balance.formula) || arg_not_used(balance.formula) ) {
    output["balance.formula"] <- list(NULL)
  } else {
    balance.formula <- as.character(balance.formula)
    bf_split <- strsplit(balance.formula, "~")
    output$balance.formula <- paste0("~ 0 +", bf_split[[1]][2])
  }
  
  if(missing(quick.balance.function) ||  arg_not_used(quick.balance.function) ) {
    output[["quick.balance.function"]] <- TRUE
  } else {
    output[["quick.balance.function"]] <- isTRUE(quick.balance.function)
  }
  if (missing(grid.length) ||  arg_not_used(grid.length) ) {
    output$grid.length <- 7L
  } else {
    output$grid.length <- as.integer(grid.length)
    if (grid.length <= 0) stop("grid.length must be greater than 0")
  }
  if (!is.null(output$lambda) && !is.null(output$delta)) {
    output["grid.length"] <- list(NULL)
  }
  
  if (missing(p) || arg_not_used(p) ) {
    output[["p"]] <- 2.0
  } else {
    stopifnot(is.numeric(p)) 
    stopifnot(p >= 1.0)
    output[["p"]] <- p
  }
  
  if ( missing(cost.online) || arg_not_used(cost.online) ) {
    output$cost.online <- "auto"
  } else {
    output$cost.online <- match.arg(cost.online,
                                             c("auto",
                                               "tensorized",
                                               "online"))
  }
  if (missing(cost.function) || arg_not_used(cost.function) ) {
    output["cost.function"] <- list(NULL)
  } else {
    output["cost.function"]  <- cost.function
  }
  if (missing(diameter) ||  arg_not_used(diameter) ) {
    output["diameter"] <- list(NULL)
  } else {
    output["diameter"]  <- as.numeric(diameter)
  }
  
  if ( is.null(torch.optimizer) ) {
    output["torch.optimizer"] <- list(NULL)
    output["solver.options"] <- list(NULL)
  } else if (missing(torch.optimizer) ) {
    output$torch.optimizer <- switch(output$opt.direction,
                                     "primal" = torch::optim_lbfgs,
                                     "dual" = torch::optim_rmsprop)
    output$solver.options <- list(...)[...names() %in% methods::formalArgs(output$torch.optimizer)]
    
    if(output$opt.direction == "dual" && is.null(output$solver.options$momentum) ) {
      output$solver.options$momentum <- 0.9
    }
  } else {
    if (!inherits(torch.optimizer, "torch_optimizer_generator")) {
      stop("torch.optimizer must be a torch_optimizer_generator function or NULL")
    }
    output$torch.optimizer <- torch.optimizer
    output$solver.options <- list(...)[...names() %in% methods::formalArgs(output$torch.optimizer)]
  }
  
  if ( missing(torch.scheduler) && missing(torch.optimizer) ) {
    if (inherits(output$torch.optimizer, "optim_lbfgs") ) {
      output$torch.scheduler <- torch::lr_reduce_on_plateau
    } else if ( inherits(output$torch.optimizer, "optim_rmsprop")) {
      output$torch.scheduler <- torch::lr_multiplicative
      output$scheduler.options <- list(lr_lambda = function(epoch) {0.99})
    }
    
  } else if (missing(torch.scheduler) || is.null(torch.scheduler) || is.null(torch.optimizer)) {
    output["torch.scheduler"] <- list(NULL)
    output["scheduler.options"] <- list(NULL)
  } else {
    if (!inherits(torch.scheduler, "lr_scheduler")) {
      stop("torch.scheduler must be a lr_scheduler function or NULL")
    }
    output$torch.scheduler <- torch.scheduler
    output$scheduler.options <- list(...)[...names() %in% methods::formalArgs(torch.scheduler)]
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
  
  if (missing(lambda.bootstrap) || arg_not_used (lambda.bootstrap)) {
    output$lambda.bootstrap <- 0.05
  } else {
    output$lambda.bootstrap <- as.numeric(lambda.bootstrap)
    stopifnot("lambda.bootstrap must be > 0"= lambda.bootstrap>0)
  }

  if (arg_not_used(tol)) {
    output$tol <- 1e-7
  } else {
    output$tol <- as.double(tol)
  }
  
  
  output$device <- cuda_device_check(device)
  output$dtype  <- cuda_dtype_check(dtype, output$device)
  
  # check combinations
  # if(output$debias && is.null(output$torch.optimizer)) {
  #   stop("Must supply a 'torch_optimizer_generator' object in options argument 'torch.optimizer' when option debias is TRUE")
  # }
  
  output$osqp.options <- list(...)[...names() %in% methods::formalArgs(osqp::osqpSettings)]
  
  # if(output$debias && output$penalty.function == "L2") {
  #   warning("No debias options with L2 penalty. Defaulting to debias=FALSE")
  #   output$debias <- FALSE
  # }
  
  # if (output$penalty.function == "L2") {
  #   output$solver.options <- list(...)[...names() %in% formalArgs(osqp::osqpSettings)]
  # }
  # if (!output$debias && output$penalty.function == "entropy" && output$cost.online != "online") {
  #   output$solver.options <- lbfgs3c_control(...)
  # }
  
  if(inherits(output$torch.optimizer, "optim_lbfgs") && output$opt.direction == "dual") {
    warning("LBFGS doesn't work with the dual optimization function. Switching to optim_rmsprop.")
    output$torch.optimizer <- torch::optim_rmsprop
    output$solver.options <- list(...)[...names() %in% methods::formalArgs(output$torch.optimizer)]
  }
  
  

  class(output) <- "cotOptions"
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
    penalty <- prob$penalty
    options["lambda"] <- list(penalty$lambda)
    options["delta"] <- list(penalty$delta)
   
  } else {
    stop("Method not found in function balanceDistributions.")
  }
    
  return(list(problem = prob,
              options = options))

}



