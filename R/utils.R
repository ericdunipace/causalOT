sqrt_mat <- function(X) {
  p <- ncol(X)
  decomp <- eigen(X)
  return(tcrossprod(decomp$vectors %*% diag(sqrt(decomp$values), p, p), decomp$vectors))
}

mahal_transform <- function(X, Y) {
  p <- ncol(X)
  decomp <- eigen(0.5 * cov(X) + 0.5 * cov(Y))
  L_inv <- tcrossprod(decomp$vectors %*% diag(1/sqrt(decomp$values), p, p), decomp$vectors)
  
  return(list(X %*% L_inv, Y %*% L_inv))
}