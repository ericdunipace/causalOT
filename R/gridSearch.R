# TODO: make gridSearch class R6 method so plays more nicely with other R6
# grid search methods

#' @include calc_weight.R
#' @include weightsClass.R
#' @include balanceFunctions.R
#' @include scmClass.R
#' @include cotClass.R

setOldClass(c("COT","R6"))
setOldClass(c("SCM","R6"))
setOldClass(c("balanceFunction", "R6"))
setOldClass(c("EntropyBW", "balanceFunction", "R6"))
setOldClass(c("SBW", "balanceFunction", "R6"))


#' gridSearch S4 class
#'
#' @slot penalty_list numeric. 
#' @slot nboot integer. 
#' @slot solver R6. 
#' @slot method character. 
#' @slot estimand character. 
#'
#' @keywords internal
setClass("gridSearch",
         slots = c(
           penalty_list = "numeric",
           nboot = "integer",
           solver = "R6",
           method = "character",
           estimand = "character"
           # ns = "integer",
           # nt = "integer"
         )
         # contains = "cotProblem"
)


setClass("ateClass",
         slots = c(
           control = "gridSearch",
           treated = "gridSearch"
         ))

gridSearch <- function(data, estimand, method, options = NULL) {
  # maybe arrange data into source and targets
  stopifnot(method %in% grid_search_methods())
  stopifnot(estimand %in% c("ATT","ATC", "ATE", "ATE.C", "ATE.T"))
  
  if (estimand == "ATE") {
   return(
     methods::new("ateClass",
        control = gridSearch(data, estimand = "ATE.C",
                             method, options),
        treated = gridSearch(data, estimand = "ATE.T",
                             method, options))
   )
  } else if (estimand %in% c("ATT","ATC", "ATE.C","ATE.T") ) {
    source_target <- data_separate(data, estimand)
    if (method %in% balancing_distributions_methods()) {
      # if(!is.null(options$balance.functions)) {
      #   options$balance.functions$data <- data_separate(options$balance.functions$data,
      #                                              estimand)
      # }
      
      problem_list <- balanceDistributions(source = source_target$source,
                                   target = source_target$target,
                                   a = source_target$a,
                                   b = source_target$b,
                                   method = method,
                                   options = options)
      prob <- problem_list$problem
      options <- problem_list$options
      
      grid <- list(lambda = options$lambda, 
                   delta = options$delta)
      grid.length <- options$grid.length
      nboot <- options$nboot
      
    } else if (method %in% balancing_function_methods()) {
      
      if(!inherits(options, "gridSearchOptions")) {
        grid_specific_opts <- do.call(gridSearchOptions, options)
      } else {
        grid_specific_opts <- options
      }
      # browser()
      prob <- bfMethods(source = source_target$source,
                target = source_target$target,
                a = source_target$a,
                b = source_target$b,
                method = method, 
                options = options)
      
      grid <- grid_specific_opts$delta
      grid.length <- length(grid)
      if(grid.length == 0) grid.length <- grid_specific_opts$grid.length
      nboot <- grid_specific_opts$nboot
      
      if(!is.null(grid_specific_opts$delta) && length(grid) == 1) nboot <- 0L
    } else {
      stop("method not found!")
    }
  } else {
    stop("estimand not found!")
  }

  return(methods::new("gridSearch",
      penalty_list = prob$gridInit(grid, grid.length),
      nboot = nboot,
      solver = prob,
      method = method,
      estimand = estimand
      ))

}

gridSearchOptions <- function(nboot = 1000L, grid.length = 20L, ...) {
  
  output <- list()
  if(missing(nboot) || is.na(nboot) || is.null(nboot)) {
    output$nboot <- 1000L
  } else {
    output$nboot <- as.integer(nboot)
  }
  if(output$nboot < 0) stop("nboot must be >= 0")
  
  if(missing(grid.length) || is.na(grid.length) || is.null(grid.length)) {
    output$grid.length <- 20L
  } else {
    output$grid.length <- as.integer(grid.length)
  }
  if(output$grid.length < 0) stop("grid.length must be >= 0")
  
  output <- c(output, list(...))
  class(output) <- c("gridSearchOptions", class(output))
  
  return(output)
}

# setGeneric("data_separate", function(data, estimand) standardGeneric("data_separate"))

data_separate <- function(data, estimand) UseMethod("data_separate")
#' Title
#'
#' @param data dataHolder. 
#' @param estimand character. 
#'
#' @keywords internal
#' @include dataHolder.R
data_separate.dataHolder <- function(data, estimand) {
  if(estimand == "ATC") {
    source = get_x1(data)
    target = get_x0(data)
    a = get_w1(data)
    b = get_w0(data)
  } else if (estimand == "ATT") {
    source = get_x0(data)
    target = get_x1(data)
    a = get_w0(data)
    b = get_w1(data)
  } else if (estimand == "ATE.C") {
    source = get_x0(data)
    target = get_x(data)
    a = get_w0(data)
    b = get_w(data)
  } else if (estimand == "ATE.T") {
    source = get_x1(data)
    target = get_x(data)
    a = get_w1(data)
    b = get_w(data)
  } else {
    stop("estimand not found!")
  }
  return(
    list(source = source, target = target,
         a = a, b = b)
  )
}
# setMethod("data_separate", signature(data = "dataHolder", estimand = "character"),
# function(data, estimand) {
#   if(estimand == "ATC") {
#     source = get_x1(data)
#     target = get_x0(data)
#     a = get_w1(data)
#     b = get_w0(data)
#   } else if (estimand == "ATT") {
#     source = get_x0(data)
#     target = get_x1(data)
#     a = get_w0(data)
#     b = get_w1(data)
#   } else if (estimand == "ATE.C") {
#     source = get_x0(data)
#     target = get_x(data)
#     a = get_w0(data)
#     b = get_w(data)
#   } else if (estimand == "ATE.T") {
#     source = get_x1(data)
#     target = get_x(data)
#     a = get_w1(data)
#     b = get_w(data)
#   } else {
#     stop("estimand not found!")
#   }
#   return(
#     list(source = source, target = target,
#          a = a, b = b)
#   )
# }
# )

# create a generic right now but don't need necessarily...
setGeneric("grid_select", function(object, w) standardGeneric("grid_select"))
setMethod("grid_select", signature(object = "ANY", w = "list"),
function(object,w) {

  if (length(w) == 1) {
    return(list(weight = w[[1]],
                idx = 1))
  }
  
  # how many bootstrap iterations
  nboot <- object@nboot
  
  # run special COT selections
  if (inherits(object@solver, "COT") ) {
    res <- object@solver$grid_search(nboot)
    
    return(res)
  }
  
  # remove weights that are infeasible
  null_tf <- sapply(w, is.null)
  if (any(null_tf)) w <- w[!null_tf]
  
  na_tf <- sapply(w, function(ww) all(is.na(ww)))
  if (any(na_tf)) w <- w[!na_tf]

  
  # C++ for loop
  metric <- bootStrap_(w,
             as.integer(nboot),
             object@solver)

  # get index of minimizing penalty parameter
  min.idx <- which.min(metric)
  
  # return values
  return(list(weight = w[[min.idx]] #selected weight
              , idx = min.idx # selected index,
              , metric = metric
              , penalty.grid = NULL
              ))
}
)


# causalWeights def -------------------------------------------------------

setMethod("causalWeights", signature(object1 = "gridSearch", object2 = "numeric"), 
          function(object1, object2, ...) {
            penalty <- list(...)[[1]]
            estimand <- object1@estimand
            
            if (estimand == "ATE.T" || estimand == "ATC") {
              w0 <- object1@solver$b
              w1 <- as.numeric(object2)
            } else if (estimand == "ATE.C" || estimand == "ATT") {
              w0 <- as.numeric(object2)
              w1 <- object1@solver$b
            } else {
              stop("estimand not found!")
            }
            if(!is.list(penalty$penalty)) penalty$penalty <- as.list(penalty$penalty)
            
            methods::new("causalWeights",
                         w0 = w0,
                         w1 = w1,
                         estimand = object1@estimand,
                         method  = object1@method,
                         penalty = penalty$penalty,
                         info    = list(metric = penalty$metric,
                                        penalty.grid = penalty$penalty.grid,
                                        gradients = NULL,
                                        hessian = NULL
                         ),
                         data = dataHolder(x = matrix(0,0,0), z = numeric(0)),
                         call = call("calc_weight"))
          }
)



# cot_solve function ------------------------------------------------------

#TODO: self object solver object@solver$solve(), like COT

#' cot_solve for gridSearch
#'
#' @param object gridSearch. 
#'
#' @return returns object of class [causalWeights][causalOT::causalWeights-class]
#' @keywords internal
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

#' cot_solve method for ateClass objects
#'
#' @param object ateClass. 
#'
#' @return object of class [causalWeights][causalOT::causalWeights-class]
#' @keywords internal
setMethod("cot_solve", signature(object = "ateClass"),
          function(object) {
            
            # "auto" shares one score across arms when both fit the same grid
            share_rbf <- inherits(object@control@solver, "COT") &&
              identical(object@control@solver$lambda_selector, "auto")
            if (share_rbf) {
              g0 <- object@control@solver$rbf_grid()
              g1 <- object@treated@solver$rbf_grid()
              share_rbf <- !is.null(g0) && !is.null(g1) &&
                identical(sort(g0), sort(g1))
            }
            if (share_rbf) {
              control <- object@control; treated <- object@treated
              control@solver$solve(); treated@solver$solve()
              p0 <- control@solver$rbf_path_eligibility()
              p1 <- treated@solver$rbf_path_eligibility()
              arm_order <- match(p0$original_index, p1$original_index)
              if (anyNA(arm_order) || !identical(p0$lambda, p1$lambda[arm_order]))
                stop("ATE arm fits have different lambda grids.")
              common <- p0$eligible & p1$eligible[arm_order] & is.finite(p0$lambda) & p0$lambda > 0
              if (!any(common)) stop("No common usable lambda candidate across ATE arms.")
              common_treated <- common[match(p1$original_index, p0$original_index)]
              shared <- control@solver$rbf_shared_scores()
              r0 <- control@solver$grid_search(.rbf.eligible = common, .rbf.shared = shared)
              r1 <- treated@solver$grid_search(.rbf.eligible = common_treated, .rbf.shared = shared)
              return(causalWeights(causalWeights(control, r0$weight, r0),
                                   causalWeights(treated, r1$weight, r1)))
            }
            # control weights targeting full sample
            cw_w0 <- cot_solve(object@control)
            
            # treated weights targeting full sample
            cw_w1 <- cot_solve(object@treated)
            
            # combine objects
            cw <- causalWeights(cw_w0, cw_w1)
            
            return(cw)
          }
)