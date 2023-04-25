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
  method <- match.arg(method, supported_methods())
  if (method == "Wasserstein") method <- "COT"
  estimand <- match.arg(estimand)

  # process the data
  data <- dataHolder(x = x, z = z, y = NA_real_, weights = weights)

  # set up the solver
  problem <- cotProblem(data, estimand, method,
                             options)

  # estimate the weights
  cw <- cot_solve(problem)
  
  # save the call and data
  # these involve copying the object which isn't ideal. 
  # may change in the future
  cw@call <- mc
  cw@data <- data
  
  # return solution
  return(cw)
}


cotProblem <- function(data, estimand, method,
                       options) {
  if (!is.list(options)) options <- list(options)
  
  if(method %in% grid_search_methods()) {
    prob <- gridSearch(data, estimand, method, options)
  } else if (method %in% likelihood_methods()) {
    prob <- likelihoodMethods(data, estimand, method, options)
  } else {
    stop("method not found!")
  }
  return(prob)
}

# S4 method to solve the respective problems
setGeneric("cot_solve", function(object) standardGeneric("cot_solve"))

setOldClass("likelihoodMethods")
setMethod("cot_solve", signature(object = "likelihoodMethods"),
function(object) {
  
  fit <- object@solver(object)
  
  cw <- methods::new("causalWeights", 
            w0 = fit$w0,
            w1 = fit$w1,
            estimand = fit$estimand,
            method = fit$method,
            penalty = list(NULL),
            info = list(NULL),
            data = object@data,
            call = call("calc_weight")
            )
  
  if(!inherits(cw, "causalWeights")) stop("Solver didn't return object of class causalWeights!")
  return(cw)
}
)
#TODO: self object solver object@solver$solve(), like COT
setOldClass("gridSearch")
setMethod("cot_solve", signature(object = "gridSearch"),
          function(object) {
            # browser()
            if (inherits(object@solver, "COT") ) {
              object@solver$solve()
              res <- object@solver$grid_search()
              cw <- causalWeights(object, res$weight, res)
              
              return(cw)
            }
            
            # set up terms for the loop
            n_penalty <- length(object@penalty_list)
            w <- vector("list", n_penalty)
            penalty <- NULL
            
            # run solver on each penalty set
            for(k in 1:n_penalty) {
              penalty <- object@penalty_list[k]
              # print(penalty)
              w[[k]] <- object@solver$solve(penalty, w[k-1])
            }
            
            # boot strap to find optimal penalty parameters
            grid_info <- grid_select(object, w)
            grid_info$penalty.grid <- object@penalty_list
              
            # store final weight and selected penalty parameters
            w_final <- grid_info$weight
            pen_final <- object@penalty_list[grid_info$idx]
            
            # create causalWeights object
            cw <- causalWeights(object, w_final, grid_info)
            
            return(cw)
          }
)

setOldClass("ateClass")
# grid search for ATE class of gridSearch objects
setMethod("cot_solve", signature(object = "ateClass"),
          function(object) {
            
            # control weights targeting full sample
            cw_w0 <- cot_solve(object@control)
            
            # treated weights targeting full sample
            cw_w1 <- cot_solve(object@treated)
            
            # combine objects
            cw <- causalWeights(cw_w0, cw_w1)
            
            return(cw)
          }
)
