#' Supported weighting methods
#'
#' @return A character vector with values "Logistic", "Probit", "SBW", "SCM", 
#' "CBPS",  "NNM", 
#' "Wasserstein" or equivalently "COT",
#'  and "None".
#' @export
#'
#' @examples
#' supported.methods()
supported.methods <- function() {
  return(c("Wasserstein", "COT",
           "Logistic", "Probit",
           "SBW", "SCM", 
           # "EBW",
           "CBPS",
           # "RKHS", "RKHS.dose", 
           "NNM", 
           "None"))
}

#' Supported optimal transport methods
#' 
#' Lists the supported OT methods. Note "COT" and "Wasserstein" are equivalent.
#'
#' @return A character vector with values "NNM","Wasserstein", "COT", and "SCM".
#' @export
#'
#' @examples
#' ot.methods()
ot.methods <- function() {
  return(c("NNM"
           ,"Wasserstein",
           "COT"
           # ,"Constrained Wasserstein"
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
#' The "Lp" method uses the simple \eqn{L_p} norm, while "RKHS" 
#' calculates the kernel for a reproducible kernel Hilbert space (RKHS).
#' 
#' @export
#'
#' @examples
#' dist.metrics()
dist.metrics <- function() {
  c("sdLp", "mahalanobis",  "Lp", "RKHS" )
}


#' Supported solvers
#'
#' @return a character vector with values "lbfgs","mosek","gurobi", and "osqp"
#' 
#' @details The solvers "mosek" and "gurobi" are commercial solvers that require software licenses. "quadprog" uses the `osqp` R package and "lbfgs" will either use `pytorch` in Python or the `lbfgs3c` package in R, which are both free.
#' 
#' @export
#'
#' @examples
#' supported.solvers()
supported.solvers <- function() {
  c("lbfgs","mosek",
    # "gurobi",
    "osqp")
}