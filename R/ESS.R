ESS.default <- function(x) {
  return(1/sum(x^2))
}

ESS.causalWeights <- function(x) {
  return(c("Control" = ESS(x$w0), "Treated" = ESS(x$w1)))
}

setOldClass("causalWeights")

#' Effective Sample Size
#'
#' @param x Either a vector of weights summing to 1 or an object of class
#' [causalWeights][causalWeights]
#' 
#' @details Calculates the effective sample size as described by Kish (1965). 
#' However, this calculation has some problems and the [PSIS][PSIS] 
#' function should be used instead.
#'
#' @return Either a number denoting the effective sample size or if `x` is of class
#' [causalWeights][causalWeights], then returns a list of both values in the treatment
#' and control groups.
#' 
#' @export
#'
#' @examples
#' x <- rep(1/100,100)
#' ESS(x)
setGeneric("ESS", function(x) UseMethod("ESS"))
setMethod("ESS", signature(x = "numeric"), ESS.default)
setMethod("ESS", signature(x = "causalWeights"), ESS.causalWeights)