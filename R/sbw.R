sbw_dual <- function(x, target, constraint) {
  
  n  <- nrow(x)
  x_m <- colMeans(x)
  x_v <- matrixStats::colVars(x)
  
  if (!is.null(nrow(target)) ) {
    if (nrow(target) > 1 )  {
      S <- sqrt(0.5 * x_v + 0.5 * matrixStats::colVars(target))
      target <- colMeans(target)
    }
  } else {
    S <- sqrt(x_v)
  }
  
  if ( constraint == 0) {
    # equivalent to (x'x)^(-1) (x' 1_n - target)
    QR <- qr(x)
    R <- qr.R(QR)
    
    beta <- 2 * (qr.coef(QR, rep(1,ncol(x))) - backsolve(R, forwardsolve(l = R, target,
                                                         transpose = TRUE)))
    
  } else {
    # use lasso in oem that can take xtx and xty
    
    XtY <- (x_m - target) #* 1/S)
    # XtX <- crossprod(scale(x, center = FALSE, scale = S))
    XtX <- crossprod(x)/4
    
    fit <- oem::oem.xtx(xtx = XtX, xty = XtY, family = "gaussian",
                        penalty = "lasso", lambda = constraint
                        , penalty.factor = S
    )
    beta <- fit$beta$lasso
  }
  
 
  
  unconst_wt  <- dual_to_wt(x, lambda = beta, n = n)
  
  return(list(weight = simplex_proj(unconst_wt), 
              unconstrained_weight = unconst_wt,
              lambda = beta))
}

dual_to_wt <- function(x, lambda, n) {
  eta <- (x %*% lambda )
  return(-(eta)/2.0 + 1.0 / n)
}
