#' causalWeights class
#'
#' @slot w0 A slot with the weights for the control group. 
#' @slot w1 The weights for the treated group. 
#' @slot gamma The trasportation matrix. If estimand is "ATE", will be a list with the transportation plan for each treatment group to balance towards the overall treatment. 
#' @slot estimand A character denoting the estimand targeted by the weights. One of "ATT","ATC", or "ATE". 
#' @slot method A character denoting the method used to estimate the weights. 
#' @slot args The other arguments used to construct the weights. 
#' 
#' @docType class
#' 
#' @export
setClass("causalWeights", slots = c(w0 = "numeric", w1 = "numeric", 
                                    gamma = "matrix",estimand = "character",
                                    method = "character", args = "list"))

#' sampleWeights class
#'
#' @slot a The sample weights for the fist group
#' @slot b The sample weights for the second group
#' @slot total The sample weights for the overall sample
#'
#' @docType class
#' 
#'
#' @export
setClass("sampleWeights", slots = c(a = "numeric", b = "numeric", total = "numeric"))

#' Estimate causal weights
#'
#' @param data Either a matrix, a data.frame, or a DataSim class. Arguments "balance.covariates" and "treatment.indicator" must be provided in the `...` arguments if data is of class data.frame or matrix.
#' @param constraint The constraints or penalties for the weights. See details.
#' @param estimand The estimand of interest. One of "ATT","ATC", or "ATE".
#' @param method The method to estimate the causal weights. Must be one of the methods returned by [supported.methods()][supported.methods()].
#' @param formula The formula for creating the design matrix used in various methods. See details.
#' @param transport.matrix Should the method calculate the transportation matrix if not done as a part of the method (TRUE/FALSE)? Default is FALSE.
#' @param grid.search Should hyperparameters be selected by a grid search? Only available for "SBW" and "COT"/"Wasserstein" methods.
#' @param ... Many additional arguments are possible depending on the chosen method. See details for more information. Arguments "balance.covariates" and "treatment.indicator" must be provided if data is of class data.frame or matrix.
#' 
#' @details 
#' We detail some of the particulars of the function arguments below.
#' 
#' ## `data`
#' The following classes are recognized by the `data` variable.
#'  
#' ### DataSim class
#' The DataSim class is provided by this package for simulations. You can pass a DataSim object (once data has been simulated) to this function
#' and it will be recognized and handled appropriately.
#' 
#' ### data.frame or matrix
#' If the `data` argument is of class `data.frame` or `matrix`, then additional arguments are necessary to pass in the dots (`...`).
#' These *must* include a vector argument `balance.covariates` and an integer or character in the `treatment.indicator` argument. The 
#' `balance.covariates` argument should be either an integer vector giving the column numbers of the covariates to balance or a character vector
#' giving the names of the columns to balance. Similarly, the `treatment.indicator` argument should be a integer giving the column number of the 
#' treatment labels or a character giving the column name.
#' 
#' ## Constraints
#' The constraint argument is used by the balancing methods like "SBW". This will specify a tolerance for basis function balance.
#' 
#' If method "COT"/"Wasserstein" is used, will specify the penalty parameter to put on the weights. For "ATT" and "ATC" estimands, must be of the form `list(penalty = ###)`, while for estimand "ATE", must be a list of length 2 specifying penalty first for the controls and then for treated: `list(list(penalty = ###), list(penalty = ###))`.
#' 
#' This argument is not needed if `grid.search` is TRUE.
#' 
#' ## Formula
#' For methods "SBW" or "COT", should be a formula object or character without a response but with the covariate functions desired. e.g., "~." includes all covariates without transformation.
#' 
#' For methods "Logistic" and "Probit", a propensity score model either as a formula object or character: "z ~.".
#' 
#' ## Additional arguments in `...`
#' In addition to the already mentioned arguments, there are several additional
#' optional arguments for the method "COT".
#' * `p`. The power of the Wasserstein distance to use.
#' * `metric`. The metric to use for the ground cost function. See [dist.metrics()] for supported distance metrics.
#' * `penalty`. What type of penalty should be used on the weights? Must be one of "entropy" or "L2".
#' * `add.divergence`. TRUE or FALSE. If TRUE, `penalty` defaults to entropy.
#'   and will calculate the Sinkhorn divergence version of Causal Optimal Transport.
#'   If choosing Sinkhorn divergences, the Python
#'    package `geomloss` must be installed.
#' * `balance.constraints`. The tolerance for the balancing basis function methods.
#' * `cost`. If the cost matrix is already calculated, you can supply this to potentially save time.
#' 
#' Additionally, methods like "SBW" and "COT" need the specification of a solver function if using balancing functions, i.e. if the `formula` argument is specified.
#' * `solver`. Should be one of "mosek" or "osqp".
#'    
#' @seealso [estimate_effect()][causalOT::estimate_effect]
#'  
#' @return An object of class [causalWeights][causalOT::causalWeights-class]
#' @export
#'
#' @examples
#' set.seed(23483)
#' n <- 2^7
#' p <- 6
#' overlap <- "low"
#' design <- "A"
#' estimate <- "ATE"
#' #### get simulation functions ####
#' data <- causalOT::Hainmueller$new(n = n, p = p, 
#'       design = design, overlap = overlap)
#'       data$gen_data()
#' weights <- calc_weight(data = data, 
#'       p = p,
#'       estimand = estimate,
#'       method = "NNM")
#' \dontrun{
#' # Needs Python package GeomLoss
#' COTweights <- calc_weight(data = data, 
#'       p = 2,
#'       constraint = list(list(penalty = 1000),
#'                         list(penalty = 10000)),
#'       estimand = estimate,
#'       method = "COT",
#'       penalty = "entropy",
#'       add.divergence = TRUE,
#'       verbose = TRUE
#'       )
#' # with basis function balancing.
#'  COTweightsBF <- calc_weight(data = data, 
#'       p = 2,
#'       constraint = list(list(penalty = 1000),
#'                         list(penalty = 10000)),
#'       estimand = estimate,
#'       method = "COT",
#'       penalty = "entropy",
#'       add.divergence = TRUE,
#'       formula = "~.",
#'       balance.constraints = 0.2,
#'       solver = "osqp",
#'       verbose = TRUE
#'       )
#' }
calc_weight <- function(data, constraint=NULL,  estimand = c("ATE","ATT", "ATC","cATE", "feasible"), 
                        method = supported.methods(),
                        formula = NULL,
                        transport.matrix = FALSE, grid.search = FALSE,
                        ...) {
  method <- match.arg(method)
  if (method == "COT") method <- "Wasserstein"
  estimand <- match.arg(estimand)
  grid.search <- isTRUE(grid.search)
  args <- list(data = data, constraint = constraint, estimand = estimand, method = method, 
               formula = formula, transport.matrix = transport.matrix, ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  if (isTRUE(args[["penalty"]] == "none") && isTRUE(grid.search == TRUE) && ((
     (args$method %in% c("Constrained Wasserstein","Wasserstein")) &&
      (isFALSE(args$add.margins) && isFALSE(args$joint.mapping))) ||
     args$method == "SCM")  ) {
    grid.search <- FALSE
  }
  
  # if(args$method == "Wasserstein" && isTRUE(args$add.divergence) ) {
  #   grid.search <- FALSE
  # }
  
  output <- if (grid.search & method %in% c("SBW","RKHS.dose","Constrained Wasserstein","Wasserstein","SCM")) {
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
                       "RKHS.dose" = "RKHS_grid_search",
                       "SCM" = "wass_grid_search"
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
  } else if (method == "CBPS") {
    # f.call <- as.call(c(list(as.name("calc_weight_CBPS")), argn))
    # eval(f.call, envir = args)
    calc_weight_CBPS(data = data, estimand = estimand, ...)
  } else if ( method != "Logistic" && method != "Probit") {
    f.call <- as.call(c(list(as.name("calc_weight_bal")), argn))
    eval(f.call, envir = args)
    # do.call("calc_weight_bal", list(data, constraint,  estimand = estimand, 
    #                                 method = method,
    #                                 ...))
  } else {
    if (estimand == "cATE") estimand <- "ATE"
    calc_weight_glm(data = data, constraint = constraint, estimand = estimand, formula = formula, method = method, ...)
  }
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, ...)
  }
  output$estimand <- estimand
  output$method <- method
  class(output) <- "causalWeights"
  return(output)
}

#' Calculate Nearest Neighbor Matching weights
#'
#' @param data 
#' @param estimand 
#' @param transport.matrix Should we calculate the transport matrix (TRUE/FALSE)? Default is FALSE.
#' @param sample_weight The sample weights
#' @param ... 
#'
#' @return a list object returned to main [calc_weight] function
#'
#' @keywords internal
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
  cost <- dots[["cost"]]
  p <- dots[["p"]]
  if (is.null(p)) p <- floor(ncol(x)/2 + 1)
  if(p < ncol(x)/2) warning("power of distance metric is not rate optimal for nearest neighbor matching (p > ncol(x)/2)")
  if (is.null(dots[["metric"]])) dots$metric <- "mahalanobis"
  if (dots[["metric"]] == "RKHS" & is.null(dots$rkhs.args) & is.null(dots[["cost"]])) {
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
                 args = addl.args,
                 estimand = estimand,
                 method = "NNM"
                 )
  
  if (isTRUE(transport.matrix)) {
    output$gamma <- calc_gamma(output, cost = cost, p = p, ...)
  }
  # output <- convert_sol(sol, estimand, method, ns["n0"], ns["n1"])
  
  return(output)
}

#' Calculate IPW using Logistic regression
#'
#' @param data 
#' @param constraint 
#' @param estimand 
#' @param ... 
#'
#' @return a list returned to the main [calc_weight()][causalOT::calc_weight] function.
#'
#' @keywords internal
calc_weight_glm <- function(data, constraint,  estimand = c("ATE","ATT", "ATC"),
                            method = c("Logistic","Probit"),
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
  method <- match.arg(method)
  fam <- switch(method,
                "Logistic" = binomial(link = "logit"),
                "Probit" = binomial(link = "probit"))
  dots$formula <- form_all_squares(dots$formula, colnames(df))
  mod <- glm(dots$formula, data.frame(z = z, df), family = fam)
  pred <- predict(mod, type = "response")
  
  if (isTRUE(constraint > 0) & isTRUE(constraint < 1)) {
    Ks  <- sort(c(constraint, 1 - constraint))
    up  <- Ks[2]
    low <- Ks[1]
    
    pred[pred > up]  <- up
    pred[pred < low] <- low
  }
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL, 
                 args = c(constraint = constraint , dots),
                 estimand = estimand, method = "Logistic"
                 )
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

# calculate balancing weights
calc_weight_bal <- function(data, constraint,  estimand = c("ATE","ATT", "ATC", "cATE", "feasible"), 
                            method = c("SBW",ot.methods()),
                            solver = supported.solvers(),
                            sample_weight = NULL,
                            add.divergence = FALSE,
                            ...) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  sample_weight <- get_sample_weight(sample_weight, get_z(data, ...))
  solver <- match.arg(solver)
  
  solve.method <- solver
  
  if(add.divergence == TRUE && method == "Wasserstein") {
    solve.method <- "div"
  }
  
  if(method != "Wasserstein" && solver == "lbfgs") {
    solver <- "osqp"
  }
  
  if(method == "Wasserstein" && isTRUE(list(...)$penalty == "entropy") && solve.method != "div" && solver != "mosek" && solver != "lbfgs") {
    solve.method <- solver <- "lbfgs"
  }
  
  if (solver == "lbfgs" && method != "Wasserstein") {
    solver <- "osqp"
    solve.method <- "osqp"
  }
  
  solve.method <- switch(solve.method,
                         "lbfgs" = 1L,
                         "div" = 2L,
                         3L)
  
  res <- switch(solve.method,
                calc_weight_lbfgs(data = data, constraint = constraint,  estimand = estimand, 
                                  method = method, sample_weight = sample_weight,
                                  solver = solver,
                                  ...),
                calc_weight_div(data = data, constraint = constraint,  estimand = estimand, 
                                method = method, sample_weight = sample_weight,
                                solver = solver,
                                ...),
                calc_weight_qp(data = data, constraint = constraint,  estimand = estimand, 
                                   method = method, sample_weight = sample_weight,
                               solver = solver,
                                   ...)
                )
  
  ns <- get_n(data, ...)
  
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  output <- convert_sol(res, estimand, method, ns["n0"], ns["n1"], sample_weight)
  output$estimand <- estimand
  output$method <- method
  dots <- list(...)
  if(method %in% ot.methods()) {
    if(is.null(dots[["p"]])) dots[["p"]] <- 2
    if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    if(!is.null(dots$cost)) dots$cost <- NULL
    if(dots$metric == "RKHS") {
      if(is.null(dots$rkhs.args)) dots$rkhs.args <- list()
      if(is.null(dots$rkhs.args$opt.method)) dots$rkhs.args$opt.method <- "stan"
      if(is.null(dots$rkhs.args$kernel)) dots$rkhs.args$kernel <- "RBF"
    }
    dots$metric <- dots$metric
    dots[["power"]] <- dots[["p"]]
    dots[["p"]] <- NULL
    dots[["add.divergence"]] <- add.divergence
  }
  output$args <- c(output$args, list(solver = solver, constraint = constraint),
                        # sol = sol$sol),
                        dots)
  
  return(output)
}

#calculate weights from the quadratic program
calc_weight_qp <- function(data, constraint, estimand,
                           method, sample_weight,
                           solver, ...) {
  solver_fun <- function(qp, solver, ...) {
    tryCatch(QPsolver(qp, solver = solver, ...),
             error = function(e) {
               warning(e$message)
               list(sol = NA_real_, dual = NULL)
             })
  }
  
  
  soc <- switch(solver,
                "mosek" = TRUE,
                FALSE)
  qp <- quadprog(data = data, constraint = constraint,  estimand = estimand, 
                 method = method, sample_weight = sample_weight,
                 soc = soc,
                 ...)
  dots <- list(...)
  
  res <- lapply(qp, solver_fun, solver = solver, ...)
  return(res)
}

# #lbfgs dual solver
# calc_weight_div <- function(data = data, constraint = constraint,  estimand = estimand, 
#                 method = method, sample_weight = sample_weight,
#                 solver = solver, penalty = c("entropy","L2"),
#                 cost = NULL,
#                 p = NULL,
#                 niter = 2000,
#                 ...) {
#   
#   
#   
# }

# calculate penalized divergence based weights
calc_weight_div <- function(data, constraint, estimand,
                            method, sample_weight = NULL,
                            solver, penalty = c("entropy", "L2"), 
                            cost = NULL,
                            add.margins = FALSE,
                            metric = dist.metrics(),
                            p = 2,
                            niter = 2000,
                            tol = 1e-5,
                            search = "LBFGS",
                            stepsize = 1e-2,
                            reach = NULL,
                            diameter = NULL,
                            scaling = 0.5, truncate = 5,
                            kernel = NULL,
                            cluster_scale = NULL, 
                            verbose = FALSE, backend='auto',
                            formula = NULL,
                            balance.constraints = NULL,
                            ...) {
  
  pd <- prep_data(data, ...)
  z <- as.numeric(pd$z)
  x.df <- pd$df[,!(colnames(pd$df) == "y")]
  x <- as.matrix(x.df)
  cn <- colnames(x.df)
  x.df <- as.data.frame(x.df)
  colnames(x) <- colnames(x.df) <- cn
  
  X0 <- x[z==0,,drop = FALSE]
  X1 <- x[z==1,,drop = FALSE]
  
  # penalty <- match.arg(penalty, c("entropy", "L2") )
  penalty <- "entropy" #L2 not work!
  
  optClass <- switch(penalty,
                     L2 = wassDivL2,
                     entropy = wassDivEnt)

  if (estimand == "ATE") {
    
    sw1 <- list(a = sample_weight$b, b = sample_weight$total)
    sw0 <- sw1
    sw0$a <- sample_weight$a
    
    if(is.null(constraint)) {
      stop("Must specify a penalty for optimal transport divergences")
    }
    
    if(is.list(constraint)) {
      if(length(constraint) == 2) {
        constraint0 <- constraint[[1]]
        constraint1 <- constraint[[2]]
      } else {
        constraint0 <- constraint
        constraint1 <- constraint
      }
    } else {
      if(length(constraint) == 1) {
        constraint0 <- constraint1 <- list(penalty = constraint)
      } else {
        constraint0 <- list(penalty = constraint[[1]])
        constraint1 <- list(penalty = constraint[[2]])
      }
      
    }
    
    if (is.list(cost) && length(cost) == 2) {
      cost0 <- cost[[1]]
      cost1 <- cost[[2]]
    } else {
      cost0 <- cost1 <- cost
    }
    
    if( length(balance.constraints) == 1) balance.constraints <- c(balance.constraints, balance.constraints)
    
    # if (is.null(stepsize)) {
    #   stepsize0 <- constraint0/1e5
    #   stepsize1 <- constraint1/1e5
    # }
    
    optimizer0 <- optClass$new(X1 = X0, X2 = x, 
                               cost = cost0,
                               prog_solver = solver,
                               lambda = constraint0$penalty,
                               add.margins = FALSE,
                               metric = metric,
                               power = p,
                               niter = niter,
                               tol = tol,
                               search = search,
                               stepsize = stepsize,
                               sample_weight = sw0,
                               reach = reach,
                               diameter = diameter,
                               scaling = scaling, truncate = truncate,
                               kernel = kernel,
                               cluster_scale = cluster_scale, 
                               debias = TRUE, 
                               verbose = verbose, backend = backend,
                               balance.function.formula = formula,
                               balance.function.delta = balance.constraints[1])
    
    optimizer1 <- optClass$new(X1 = X1, X2 = x, 
                               cost = cost1,
                               prog_solver = solver,
                               lambda = constraint1$penalty,
                               add.margins = FALSE,
                               metric = metric,
                               power = p,
                               niter = niter,
                               tol = tol,
                               search = search,
                               stepsize = stepsize,
                               sample_weight = sw1,
                               reach = reach,
                               diameter = diameter,
                               scaling = scaling, truncate = truncate,
                               kernel = kernel,
                               cluster_scale = cluster_scale, 
                               debias = TRUE, 
                               verbose = verbose, backend = backend,
                               balance.function.formula = formula,
                               balance.function.delta = balance.constraints[2])
    
    
    cg(optimizer0, verbose = verbose)
    cg(optimizer1, verbose = verbose)
    
    return(list(list(sol = optimizer0$return_cw(),
                     dual = optimizer0$get_param()) , 
                list(sol = optimizer1$return_cw(),
                     dual = optimizer1$get_param())))
    
  } else if (estimand == "ATC") {
    sample_weight[c("a","b")] <- sample_weight[c("b","a")]
    if(!is.list(constraint)) constraint <- list(penalty = constraint)
    if(!is.null(cost)) cost <- t(cost)
    optimizer <- optClass$new(X1 = X1, X2 = X0, 
                              cost = cost,
                              prog_solver = solver,
                              lambda = constraint$penalty,
                              add.margins = FALSE,
                              metric = metric,
                              power = p,
                              niter = niter,
                              tol = tol,
                              search = search,
                              stepsize = stepsize,
                              sample_weight = sample_weight,
                              reach = reach,
                              diameter = diameter,
                              scaling = scaling, truncate = truncate,
                              kernel = kernel,
                              cluster_scale = cluster_scale, 
                              debias = TRUE, 
                              verbose = verbose, backend = backend,
                              balance.function.formula = formula,
                              balance.function.delta = balance.constraints)
    
    cg(optimizer, verbose = verbose)
    
    return(list(list(sol = optimizer$return_cw(), 
                     dual = optimizer$get_param())))
           
    
  } else if (estimand == "ATT") {
    if(!is.list(constraint)) constraint <- list(penalty = constraint)
    optimizer <- optClass$new(X1 = X0, X2 = X1, 
                              cost = cost,
                              prog_solver = solver,
                              lambda = constraint$penalty,
                              add.margins = FALSE,
                              metric = metric,
                              power = p,
                              niter = niter,
                              tol = tol,
                              search = search,
                              stepsize = stepsize,
                              sample_weight = sample_weight,
                              reach = reach,
                              diameter = diameter,
                              scaling = scaling, truncate = truncate,
                              kernel = kernel,
                              cluster_scale = cluster_scale, 
                              debias = TRUE, 
                              verbose = verbose, backend = backend,
                              balance.function.formula = formula,
                              balance.function.delta = balance.constraints)
    
    cg(optimizer, verbose = verbose)
    
    return(list(list(sol = optimizer$return_cw(),
                     dual = optimizer$get_param())))
  }
  
  
}

calc_weight_lbfgs <- function(data, constraint,  estimand, 
                              method, sample_weight,
                              solver, penalty = "entropy", metric = "mahalanobis",
                              p = 2, cost = NULL,
                              formula = NULL,
                              balance.constraints = NULL,
                              ...) {
  
  pd <- prep_data(data, ...)
  z <- as.numeric(pd$z)
  x.df <- pd$df[,!(colnames(pd$df) == "y")]
  x <- as.matrix(x.df)
  cn <- colnames(x.df)
  x.df <- as.data.frame(x.df)
  colnames(x) <- colnames(x.df) <- cn
  
  X0 <- x[z==0,,drop = FALSE]
  X1 <- x[z==1,,drop = FALSE]
  
  # penalty <- match.arg(penalty)
  penalty <- match.arg(penalty, c("entropy", "L2") )
  method <- match.arg(method, c("Wasserstein"))
  
  if (estimand == "ATE") {
    
    sw1 <- list(a = sample_weight$b, b = sample_weight$total)
    sw0 <- sw1
    sw0$a <- sample_weight$a
    
    if(is.null(constraint)) {
      stop("Must specify a penalty for optimal transport divergences")
    }
    
    if(is.list(constraint)) {
      if(length(constraint) == 2) {
        constraint0 <- constraint[[1]]$penalty
        constraint1 <- constraint[[2]]$penalty
      } else {
        constraint0 <- constraint$penalty
        constraint1 <- constraint$penalty
      }
    } else {
      constraint0 <- constraint1 <- constraint
    }
    
    if (is.list(cost) && length(cost) == 2) {
      cost0 <- cost[[1]]
      cost1 <- cost[[2]]
    } else {
      cost0 <- cost1 <- cost
    }
    
    sol0 <- dual_opt(x = X0, target = x, 
             init = NULL,
             sample_weights = sw0, 
             method = method,
             penalty = penalty,
             wasserstein = list(metric = metric,
                                power = p,
                                cost = cost0,
                                lambda = constraint0),
             balance = list(balance.functions = NULL,
                            formula = formula,
                            balance.constraints = balance.constraints),
             ...)
    
    
    sol1 <- dual_opt(x = X1, target = x, 
                     init = NULL,
                     sample_weights = sw1, 
                     method = method,
                     penalty = penalty,
                     wasserstein = list(metric = metric,
                                        power = p,
                                        cost = cost1,
                                        lambda = constraint1),
                     balance = list(balance.functions = NULL,
                                    formula = formula,
                                    balance.constraints = balance.constraints),
                     ...)
    
    return(list(list(sol = sol0$weight,
                     dual = sol0$dual,
                     penalty = penalty) , 
                list(sol = sol1$weight,
                     dual = sol1$dual,
                     penalty = penalty))
    )
    
  } else if (estimand == "ATC") {
    sample_weight[c("a","b")] <- sample_weight[c("b","a")]
    if(!is.list(constraint)) constraint <- list(penalty = constraint)
    sol <- dual_opt(x = X1, target = X0, 
                     init = NULL,
                     sample_weights = sample_weight, 
                     method = method,
                     penalty = penalty,
                     wasserstein = list(metric = metric,
                                        power = p,
                                        cost = cost,
                                        lambda = constraint$penalty),
                     balance = list(balance.functions = NULL,
                                    formula = formula,
                                    balance.constraints = balance.constraints),
                     ...)
    
    return(list(list(sol = sol$weight, dual = sol$dual,
                     penalty = penalty)))
    
    
  } else if (estimand == "ATT") {
    if(!is.list(constraint)) constraint <- list(penalty = constraint)
    sol <- dual_opt(x = X0, target = X1, 
                    init = NULL,
                    sample_weights = sample_weight, 
                    method = method,
                    penalty = penalty,
                    wasserstein = list(metric = metric,
                                       power = p,
                                       cost = cost,
                                       lambda = constraint$penalty),
                    balance = list(balance.functions = NULL,
                                   formula = formula,
                                   balance.constraints = balance.constraints),
                    ...)
    
    return(list(list(sol = sol$weight, dual = sol$dual,
                     penalty = penalty)))
    
  }
  
  
  
}

# calculate RKHS weights
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
    output$w1 <- renormalize(sol[[1]]$sol[pd$z == 1])
  } else if (estimand == "ATT") {
    output$w0 <- renormalize(sol[[1]]$sol[pd$z == 0])
    output$w1 <- sample_weight$b
  } else if (estimand == "cATE") {
    output$w0 <- renormalize(sol[[1]]$sol[pd$z == 0])
    output$w1 <- renormalize(sol[[2]]$sol[pd$z == 1])
  } else if (estimand == "ATE") {
    output$w0 <- renormalize(sol[[1]]$sol[pd$z == 0])
    output$w1 <- renormalize(sol[[1]]$sol[pd$z == 1])
  }
  return(output)
}

# calculate CBPS weights
calc_weight_CBPS <- function(data, formula,  estimand = c("ATE","ATT", "ATC"),
                              niter = 1000, sample_weight = NULL, ...) {
  
  ATT.flag <- switch(estimand,
                     "ATT" = 1,
                     "ATC" = 2,
                     "ATE" = 0)
  
  if (missing(formula) || is.null(formula)) {
    formula <- "z ~."
  }

  dat <- prep_data(data, ...)
  
  sw <- get_sample_weight(sample_weight, dat$z)
  
  if ( !is.null(attr(dat$df, "outcome"))) {
    dat$df$y <- NULL
  }
  z <- dat$z
  x <- dat$df
  
  cbp.dat <- cbind(z = z, x)
  
  fit <- CBPS::CBPS(formula = formula, data = cbp.dat, 
             ATT = ATT.flag, iterations = niter,
             standardize = TRUE, method = "over",
             sample.weights = sw$total,
             twostep = TRUE, ...)
  
  pred <- fit$weights
  
  
  output <- list(w0 = pred[z == 0], w1 = pred[z == 1], gamma = NULL, estimand = estimand, method = "CBPS",
                 args = list(...))
  
  return(output)
  
}

# return list if error thrown from solver
calc_weight_error <- function(n0 = NULL, n1 = NULL) {
  
  if (is.null(n0) || is.null(n1)) {
    return(list(w0 = NA_real_,
                w1 = NA_real_))
  }
  return(list(w0 = rep(NA_real_, n0),
              w1 = rep(NA_real_, n1)))
}

# convert solution from solver
convert_sol <- function(res, estimand, method, n0, n1, sample_weight) {
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  stopifnot(is.list(res))
  stopifnot(inherits(sample_weight, "sampleWeights"))
  
  if ( method %in% c("Wasserstein", "Constrained Wasserstein","SCM") ) {
    if (estimand == "ATC") {
      if(inherits(res[[1]]$sol, "causalWeights")) {
        output <- res[[1]]$sol
        output$w0 <- res[[1]]$sol$w1
        output$w1 <- res[[1]]$sol$w0
        output$gamma <- t(res[[1]]$sol$gamma)
        output$estimand <- "ATC"
        
        return(output)
      } else {
        sol <- res[[1]]$sol[1:(n0*n1)]
        output$gamma <- matrix(sol, n0, n1, byrow = TRUE) #matrix(sol[[1]]$result, n0, n1)
        output$w0 <- sample_weight$a
        output$w1 <- colSums(output$gamma)
        dual <- res[[1]]$dual
      }
      
    } else if (estimand == "ATT") {
      if(inherits(res[[1]]$sol, "causalWeights")) return(res[[1]]$sol)
      sol <- res[[1]]$sol[1:(n0*n1)]
      
      output$gamma <- matrix(sol, n0, n1) #matrix(sol[[1]]$result, n0, n1)
      output$w0 <- rowSums(output$gamma)
      output$w1 <- sample_weight$b
      dual <- res[[1]]$dual
    } else if (estimand == "cATE") {
      sol2 <- res[[1]]$sol[1:(n0*n1)]
      sol1 <- res[[2]]$sol[1:(n0*n1)]
      
      output$w0 <- rowSums(matrix(sol1, n0, n1)) #matrix(sol[[2]]$result, n0, n1)
      output$w1 <- colSums(matrix(sol2, n0, n1, byrow = TRUE)) ##matrix(sol[[1]]$result, n0, n1)
      dual <- list(res[[2]]$dual,
                   res[[1]]$dual)
    } else if (estimand == "ATE") {
      if(inherits(res[[1]]$sol, "causalWeights") &&  inherits(res[[2]]$sol, "causalWeights")) {
        output <- res[[1]]$sol
        output$w0 <- res[[1]]$sol$w0
        output$w1 <- res[[2]]$sol$w0
        output$gamma <- list(w0 = res[[1]]$sol$gamma,
                             w1 = res[[2]]$sol$gamma)
        output$args$dual <- list(w0 = res[[1]]$sol$args$dual, 
                                 w1 = res[[2]]$sol$args$dual)
        output$estimand <- "ATE"
      } else {
        N <- n0 + n1
        sol1 <- renormalize(res[[1]]$sol[1:(n0*N)])
        sol2 <- renormalize(res[[2]]$sol[1:(n1*N)])
        output$w0 <-  rowSums(matrix(sol1, n0, N)) #matrix(sol[[1]]$result, n0, n1)
        output$w1 <-  rowSums(matrix(sol2, n1, N)) #matrix(sol[[2]]$result, n0, n1)
        #note both are rowSums here
        output$gamma <- list(w0 =  matrix(sol1, n0, N),
                             w1 =  matrix(sol2, n1, N))
      }
      
      dual <- list(res[[1]]$dual,
                   res[[2]]$dual)
    }
  } else {
    if (estimand == "ATT") {
      output$w0 <- res[[1]]$sol #renormalize(sol[[1]]$result)
      output$w1 <- sample_weight$b
      dual <- res[[1]]$dual
    } else if (estimand == "ATC") {
      output$w0 <- sample_weight$a
      output$w1 <- res[[1]]$sol #renormalize(sol[[1]]$result)
      dual <- res[[1]]$dual
    } else if (estimand == "feasible") {
      output$w0 <- renormalize(res[[1]]$sol[1:n0]) #renormalize(sol[[1]]$result[1:n0])
      output$w1 <- renormalize(res[[1]]$sol[n0 + 1:n1]) #renormalize(sol[[1]][n0 + 1:n1])
      dual <- res[[1]]$dual
    } else if (estimand == "cATE") {
      output$w0 <- renormalize(res[[1]]$sol) #renormalize(sol[[1]]$result)
      output$w1 <- renormalize(res[[2]]$sol) #renormalize(sol[[2]]$result)
      dual <- list(res[[1]]$dual,
                   res[[2]]$dual)
    } else if (estimand == "ATE") {
      output$w0 <- renormalize(res[[1]]$sol) #renormalize(sol[[1]]$result)
      output$w1 <- renormalize(res[[2]]$sol) #renormalize(sol[[2]]$result)
      dual <- list(res[[1]]$dual,
                   res[[2]]$dual)
    }
  }
  penalty <- res[[1]]$penalty
  output$args <- list(dual = dual, penalty = penalty)
  return(output)
}

# combine ATT and ATC estimates to ATE
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
    p <- dots[["p"]]
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
    output$gamma <- calc_gamma(weights = output, cost = cost, p = p, ...)
  }
  
  return(output)
}

# calculate transport matrix
calc_gamma <- function(weights, ...) {
  if (!is.null(weights$gamma)) return(weights$gamma)
  dots <- list(...)
  n1 <- length(weights$w1)
  n0 <- length(weights$w0)
  if(!is.null(dots$cost) & !is.null(dots[["p"]])) {
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
      if( is.null(dots$niter)) {
        niter <- 1e6
      } else {
        niter <- dots$niter
      }
      tplan <- approxOT::transport_plan_given_C(mass_x = a,
                                    mass_y = b,
                                    p = dots[["p"]],
                                    cost = cost,
                                    method = "exact", niter = niter)
      
      temp_gamma[cbind(tplan[[1]], tplan[[2]])] <- tplan[[3]]
      
    }
    gamma[nzero_row, nzero_col] <- temp_gamma
    # transport::transport.default
    
  } else {
    gamma <- NULL
  }
  return(gamma)
}


ate_sample_weight <- function(sw, ...) {
  if (inherits(sw, "sampleWeights")) {
    b <- sw$total
    a0 <- sw$a
    a1 <- sw$b
    
    n0 <- length(a0)
    n1 <- length(a1)
    n <- n0 + n1
    
    output0 <- list(a = a0, b = b,
                    total = renormalize(c(a0 * n0, b * n)))
    
    output1 <- list(a = a1, b = b,
                    total = renormalize(c(a1 * n1, b * n)))
    
    class(output1)  <- class(output0) <- "sampleWeights"
    output <- list(output0, output1)
    return(output)
  } else {
    stop("Must be of class sampleWeights")
  }
}

switch_sample_weight <- function(sw,z) {
  if (inherits(sw, "sampleWeights")) {
    output <- list(a = sw$b, b = sw$a,
                   total = sw$total)
    class(output) <- "sampleWeights"
    return(output)
  } else {
    stop("Must be of class sampleWeights")
  }
}

# Get sample weights 
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

# get treatment indicators for DataSim
get_z.DataSim <- function(data,...) {
  return(data$get_z())
}

# get treatment indicators for data.frame
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

# get treatment indicators for matrix
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

# get number of observations for DataSim
get_n.DataSim <- function(data,...) {
  return(data$get_n())
}

# get number of observations for data.frame
get_n.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z == 0), n1 = sum(df$z == 1))
  return(ns)
}

# get number of observations for matrix
get_n.matrix <- function(data,...) {
  df <- prep_data(data,...)
  ns <- c(n0 = sum(df$z == 0), n1 = sum(df$z == 1))
  return(ns)
}

# get covariate dimension for DataSim
get_p.DataSim <- function(data,...) {
  return(data$get_p())
}

# get covariate dimension for data.frame
get_p.data.frame <- function(data,...) {
  df <- prep_data(data,...)
  if(!is.null(df$df$y)){
    return(ncol(df$df) - 1)
  } else {
    return(ncol(df$df))
  }
}

# get covariate dimension for matrix
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
