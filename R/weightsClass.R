
#' causalWeights class
#'
#' @slot w0 A slot with the weights for the control group with \eqn{n_0} entries. Weights sum to 1.
#' @slot w1 The weights for the treated group with \eqn{n_1} entries. Weights sum to 1.
#' @slot estimand A character denoting the estimand targeted by the weights. One of "ATT","ATC", or "ATE".
#' @slot method A character denoting the method used to estimate the weights.
#' @slot penalty A list or the selected penalty parameters, if relevant.
#' @slot data The dataHolder object containing the original data.
#' @slot call The call used to construct the weights.
#'
#' @docType class
#' 
#' @details This object is returned by the `calc_weight` function in this package. The slots can be accessed as any S4 object. There is no publicly accessible constructor function. 
#' 
#'
#' @export
setClass("causalWeights", slots = c(w0 = "numeric", w1 = "numeric",
                                    estimand = "character",
                                    method = "character", 
                                    penalty = "list",
                                    data = "dataHolder",
                                    call = "call"))


setGeneric("causalWeights", function(object1, object2, ...) standardGeneric("causalWeights"))

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
  
  methods::new("causalWeights",
      w0 = w0,
      w1 = w1,
      estimand = object1@estimand,
      method = object1@method,
      penalty = list(penalty),
      data = dataHolder(x = matrix(0,0,0), z = numeric(0)),
      call = call("calc_weight"))
}
)

setMethod("causalWeights", signature(object1 = "causalWeights", object2 = "causalWeights"), 
function(object1, object2, ...) {
  stopifnot(object1@estimand == "ATE.C")
  stopifnot(object2@estimand == "ATE.T")
  
  methods::new("causalWeights",
    w0 = object1@w0,
    w1 = object2@w1,
    estimand = "ATE",
    method = object1@method,
    penalty = list(w0 = object1@penalty[[1]],
                   w1 = object2@penalty[[1]]),
    data = object1@data,
    call = object1@call
  )
  
}
)
