# TODO: walk through methods and make everything an R6 class

#' Estimate causal weights
#'
#' @param x A numeric matrix of covariates. You can also pass an object of class [causalOT::dataHolder] or [causalOT::DataSim], which will make argument `z` not necessary,
#' @param z A binary treatment indicator.
#' @param estimand The estimand of interest. One of "ATT","ATC", or "ATE".
#' @param method The method to estimate the causal weights. Must be one of the methods returned by [supported_methods()][supported_methods()].
#' @param options The options for the solver. Specific options depend on the solver you will be using and you can use the solver specific options functions as detailed below..
#' @param weights The sample weights. Should be `NULL` or have a weight for each observation in the data. Normalized to sum to one.
#' @param ... Not used at this time.
#' 
#'
#' @details
#' We detail some of the particulars of the function arguments below.
#' 
#' ### Causal Optimal Transport (COT)
#' This is the.main method of the package. This method relies on various solvers depending on the particular options chosen. Please see [cotOptions()][cotOptions] for more details.
#' 
#' ### Energy Balancing Weights (EnergyBW)
#' This is equivalent to COT with an infinite penalty parameter, `options(lambda = Inf)`. Uses the same solver and options as COT, [cotOptions()][cotOptions].
#' 
#' ### Nearest Neighbor Matching with replacement (NNM)
#' This is equivalent to COT with a penalty parameter = 0, `options(lambda = 0)`. Uses the same solver and options as COT, [cotOptions()][cotOptions].
#' 
#' ### Synthetic Control Method (SCM)
#' The SCM method is equivalent to an OT problem from a different angle. See [scmOptions()][scmOptions()].
#' 
#' ### Entropy Balancing Weights (EntropyBW)
#' This method balances chosen functions of the covariates specified in the data argument, `x`. See [entBWOptions()][entBWOptions()] for more details. Hainmueller (2012).
#' 
#' ### Stable Balancing Weights (SBW)
#' Entropy Balancing Weights with a different penalty parameter, proposed by Zuizarreta (2012). See [sbwOptions()][sbwOptions()] for more details
#' 
#' ### Covariate Balancing Propensity Score (CBPS)
#' The CBPS method of Imai and Ratkovic. Options argument is passed to the function [CBPS()][CBPS::CBPS()].
#' 
#' ### Logistic Regression or Probit Regression
#' The main methods historically for implementing inverse probability weights. Options are passed directly to the `glm` function from `R`.
#'
#' @seealso [estimate_effect()][estimate_effect()]
#'
#' @return An object of class [causalWeights][causalOT::causalWeights-class]
#' @export
#'
#' @examples
#' set.seed(23483)
#' n <- 2^5
#' p <- 6
#' #### get data ####
#' data <- Hainmueller$new(n = n, p = p)
#' data$gen_data()
#' x <- data$get_x()
#' z <- data$get_z()
#' 
#' if (torch::torch_is_installed()) {
#' # estimate weights
#' weights <- calc_weight(x = x,
#'                                  z = z, 
#'                                  estimand = "ATE",
#'                                  method = "COT",
#'                                  options = list(lambda = 0))
#' #we can also use the dataSim object directly
#' weightsDS <- calc_weight(x = data,
#'                                  z = NULL,
#'                                  estimand = "ATE",
#'                                  method = "COT",
#'                                  options = list(lambda = 0))
#' all.equal(weights@w0, weightsDS@w0)
#' all.equal(weights@w1, weightsDS@w1)
#' }
calc_weight <- function(x, z,
                        estimand = c("ATC","ATT","ATE"),
                        method = supported_methods(),
                        options = NULL, weights = NULL, ...) {
  # preserve call
  mc     <- match.call()
  
  # make sure method and estimand are supported
  stopifnot("Estimand must be supplied in `estimand` argument."=!missing(estimand))
  method <- match.arg(method, supported_methods())
  if (method == "Wasserstein") method <- "COT"
  estimand <- match.arg(estimand)

  # Check the labels before dataHolder coerces them to integers.
  z_binary <- !is.matrix(x) ||
    ((is.numeric(z) || is.logical(z)) && !anyNA(z) && all(z %in% c(0, 1)))

  # process the data
  data <- dataHolder(x = x, z = z, y = NA_real_, weights = weights)

  selection <- NULL
  if (method == "COT") {
    if (!inherits(options, "cotOptions")) options <- do.call(cotOptions, as.list(options))
    # "auto" predicts binary treatment from the pooled rows. It needs uniform
    # observation weights and the standard cost; otherwise use the bootstrap.
    n <- length(data@z)
    auto_ok <- z_binary &&
      all(tabulate(data@z + 1L, 2L) >= 2L) &&
      all(is.finite(data@weights)) &&
      max(abs(data@weights - 1 / n)) <= 1e-10 / n &&
      isTRUE(options$p == 2) && is.null(options$cost.function)
    if (!auto_ok) {
      if (identical(options$lambda, "auto")) {
        options["lambda"] <- list(NULL)
        if (is.null(options$grid.length)) options$grid.length <- 7L
      }
      if (identical(options$lambda.selector, "auto")) options$lambda.selector <- "bootstrap"
    }
    if (identical(options$lambda, "auto")) {
      selection <- cot_select_trt_rbf(data, options)
      options$lambda <- selection$lambda
    }
    # The tuning rows are the original pooled observations for every estimand.
    if (identical(options$lambda.selector, "auto")) {
      options$lambda.selection.data <- list(x = data@x, z = data@z)
    }
  }

  # set up the solver
  problem <- cotProblem(data, estimand, method,
                             options)

  # estimate the weights
  cw <- cot_solve(problem)
  
  # save the call and data
  # these involve copying the object which isn't ideal. 
  # may change in the future
  if (is.null(selection) && method == "COT") {
    metrics <- if (estimand == "ATE") cw@info$w0$metric else cw@info$metric
    if (is.list(metrics)) selection <- metrics$lambda.details
  }
  if (!is.null(selection)) cw@info$lambda_selection <- selection
  cw@call <- mc
  cw@data <- data
  
  # return solution
  return(cw)
}


cotProblem <- function(data, estimand, method,
                       options) {
  if (!is.list(options)) options <- list(options)
  
  if(method %in% grid_search_methods()) {
    prob <- gridSearch(data, estimand, method, options) #defined in gridSearch.R
  } else if (method %in% likelihood_methods()) {
    prob <- likelihoodMethods(data, estimand, method, options) #defined in likelihoodClass.R
  } else {
    stop("method not found!")
  }
  return(prob)
}

# S3 method for R6 objects
cot_solve <- function(object) { UseMethod("cot_solve") }

# # S4 method to solve the respective problems
setGeneric("cot_solve", function(object) standardGeneric("cot_solve"))



