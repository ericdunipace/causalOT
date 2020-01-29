cost_calc_lp <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise")) {
  
  dir <- match.arg(direction)
  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if(dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("Dim of X and Y should be equal to have same dimension. They can different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  } else {
    if (nrow(X) != nrow(Y)) {
      stop("Rows of X and Y should be equal to have same dimension. Observations should be unique by column")
    }
  }
  stopifnot(ground_p > 0)
  
 return(causalOT::cost_calculation_(X,Y,as.double(ground_p))) 
}