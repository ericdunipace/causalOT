# setOldClass("causalWeights")

#' Effective Sample Size
#'
#' @param x Either a vector of weights summing to 1 or an object of class
#' [causalWeights][causalOT::causalWeights-class]
#' 
#' @details Calculates the effective sample size as described by Kish (1965). 
#' However, this calculation has some problems and the [PSIS()]
#' function should be used instead.
#'
#' @return Either a number denoting the effective sample size or if `x` is of class
#' [causalWeights][causalOT::causalWeights-class], then returns a list of both values in the treatment
#' and control groups.
#' 
#' @seealso [PSIS()][causalOT::PSIS]
#' 
#' @export
#' 
#' @docType methods
#' @rdname ESS
#'
#' @examples
#' x <- rep(1/100,100)
#' ESS(x)
setGeneric("ESS", def = function(x) standardGeneric("ESS"))


#' @describeIn ESS default ESS method for numeric vectors
#' @method ESS default
setMethod("ESS", signature(x = "numeric"), 
    function(x) {
      return(1/sum(x^2))
    }          
)

#' @describeIn ESS ESS method for objects of class [causalWeights][causalOT::causalWeights-class]
#' @method ESS causalWeights
#' @include weightsClass.R
setMethod("ESS", signature = signature(x = "causalWeights"),
function(x) {
  return(c("Control" = ESS(x@w0), "Treated" = ESS(x@w1)))
}
)
