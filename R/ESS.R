ESS.default <- function(x) {
  return(1/sum(x^2))
}

ESS.causalWeights <- function(x) {
  return(c("Control" = ESS(x$w0), "Treated" = ESS(x$w1)))
}

setOldClass("causalWeights")
setGeneric("ESS", function(x) UseMethod("ESS"))
setMethod("ESS", signature(x = "numeric"), ESS.default)
setMethod("ESS", signature(x = "causalWeights"), ESS.causalWeights)