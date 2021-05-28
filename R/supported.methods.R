#' Title
#'
#' @return
#' @export
#'
#' @examples
supported.methods <- function() {
  return(c("Logistic","SBW", "SCM", 
           # "EBW",
           "CBPS",
           "RKHS", "RKHS.dose", "NNM", 
           "Wasserstein",
           "Constrained Wasserstein",
           "None"))
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
ot.methods <- function() {
  return(c("NNM"
           ,"Wasserstein"
           ,"Constrained Wasserstein"
           ,"SCM"
           # ,"EBW"
           ))
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
dist.metrics <- function() {
  c("sdLp", "mahalanobis",  "Lp", "RKHS" )
}