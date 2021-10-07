#' Supported weighting methods
#'
#' @return A character vector with values "Logistic","SBW", "SCM", 
#' "CBPS", "RKHS", "RKHS.dose", "NNM", "Wasserstein",
#' "Constrained Wasserstein", and "None".
#' @export
#'
#' @examples
#' supported.methods()
supported.methods <- function() {
  return(c("Wasserstein",
           "Logistic", "Probit",
           "SBW", "SCM", 
           # "EBW",
           "CBPS",
           "RKHS", "RKHS.dose", "NNM", 
           
           # "Constrained Wasserstein",
           "None"))
}

#' Supported optimal transport methods
#'
#' @return A character vector with values "NNM","Wasserstein", and "SCM".
#' @export
#'
#' @examples
#' ot.methods()
ot.methods <- function() {
  return(c("NNM"
           ,"Wasserstein"
           ,"Constrained Wasserstein"
           ,"SCM"
           # ,"EBW"
           ))
}

#' Supported distance metrics
#'
#' @return a character vector with values "sdLp", "mahalanobis",  "Lp", and "RKHS"
#' 
#' @details The "sdLp" method uses a metric with distances normalized by the standard
#' deviation of each of the covariates `((x[,j]-y[,j])/sd(c(x[,j],y[,j])))^p`, where `x` and `y` are the data matrices in each group and `j` is a column in each matrix.
#' 
#' The "mahalanobis" metric is related except it normalizes by the full 
#' variance-covariance matrix. Be  warned that neither "sdLp" or "mahalanobis"
#' may make sense for binary covariates and care should be taken. 
#' 
#' The "Lp" method uses the simple \math{L_p} norm, while "RKHS" 
#' calculates the kernel for a reproducible kernel Hilbert space (RKHS).
#' 
#' @export
#'
#' @examples
#' dist.metrics()
dist.metrics <- function() {
  c("sdLp", "mahalanobis",  "Lp", "RKHS" )
}