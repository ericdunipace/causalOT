update_wp_tol <- function(x, new_val = NULL, p=2) {
  stopifnot(inherits(x, "OP"))
  stopifnot(is.numeric(new_val))
  LC <- ROI::constraints(x)
  LC$rhs[1] <- new_val^p
  ROI::constraints(x) <- LC
  return(x)
}