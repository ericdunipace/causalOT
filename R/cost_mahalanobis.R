cost_mahalanobis <- function(X, Y, ground_p = 2, direction = c("rowwise", "colwise")) {
  
  dir <- match.arg(direction)
  
  if(dir == "rowwise") {
    if (ncol(X) != ncol(Y)) {
      stop("Dim of X and Y should be equal to have same dimension. They can different numbers of observations")
    }
    X <- t(X)
    Y <- t(Y)
  }
  
  if (!is.double(ground_p) ) ground_p <- as.double(ground_p)
  if (nrow(X) != nrow(Y)) {
    stop("Rows of X and Y should be equal to have same dimension. Observations should be unique by column")
  }
  
  return( cost_mahal_(X, Y, ground_p) )
}