cost_calc_lp <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise")) {
  
  dir <- match.arg(direction)
  
  if(dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("Dim of X and Y should be equal to have same dimension. They can different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  }
 return(limbs::cost_calc(X,Y,ground_p)) 
}