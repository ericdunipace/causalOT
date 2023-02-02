
# grid search methods
setOldClass("cotProblem")
setOldClass(c("SBW","R6"))
setOldClass(c("COT","R6"))
setOldClass(c("EntropyBW","R6"))
setOldClass(c("SCM","R6"))
setClass("gridSearch",
         slots = c(
           penalty_list = "numeric",
           nboot = "integer",
           solver = "R6",
           method = "character",
           estimand = "character"
           # ns = "integer",
           # nt = "integer"
         ),
         contains = "cotProblem"
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
     new("ateClass",
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

  return(new("gridSearch",
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

setGeneric("data_separate", function(data, estimand) standardGeneric("data_separate"))
setMethod("data_separate", signature(data = "dataHolder", estimand = "character"),
function(data, estimand) {
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
)

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
  return(list(weight = w[[min.idx]], #selected weight
              idx = min.idx # selected index
              ))
}
)
