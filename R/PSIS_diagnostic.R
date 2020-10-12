PSIS <- function(x, r_eff = NULL,...) {UseMethod("PSIS")}

PSIS.default <- function(x, r_eff = NULL, ...) {
  
  pos <- x>0
  lx.p <- log(x[pos])
  
  if(is.null(r_eff)) r_eff <- 1
  
  res <- loo::psis(lx.p, r_eff = r_eff, ...)
  res$weights <- rep(0, length(x))
  res$weights[pos] <- exp(res$log_weights)
  return(res)
}

PSIS.causalWeights <- function(x, r_eff = NULL, ...) {

  if(!is.null(r_eff)) {
    if(length(r_eff) == 1) r_eff <- rep(r_eff, 2)
    stopifnot(length(r_eff) == 2)
  } else {
    r_eff <- c(1,1)
  }
  res <- list(w0 = NULL,
              w1 = NULL)
  
  if(is.null(x$estimand)) {
    res$w0 = PSIS.default(x$w0, r_eff = r_eff[1],...)
    res$w1 = PSIS.default(x$w1, r_eff = r_eff[2],...)
  } else if(x$estimand == "ATT") {
    res$w0 <- PSIS.default(x$w0, r_eff = r_eff[1],...)
    res$w1 <- list(diagnostics = list(pareto_k = NA, n_eff = length(x$w1)))
  } else if (x$estimand == "ATC") {
    res$w0 <- list(diagnostics = list(pareto_k = NA, n_eff = length(x$w0)))
    res$w1 = PSIS.default(x$w1, r_eff = r_eff[2],...)
  } else if (x$estimand == "ATE" | x$estimand == "feasible" | x$estimand == "cATE") {
    res$w0 = PSIS.default(x$w0, r_eff = r_eff[1],...)
    res$w1 = PSIS.default(x$w1, r_eff = r_eff[2],...)
  }
  
  class(res) <- "causalPSIS"
                   
  return(res)
}

PSIS.list <- function(x, r_eff = NULL, ...) {
  
  if(!is.null(r_eff)) {
    if(length(r_eff) == 1) r_eff <- rep(r_eff, length(x))
    stopifnot(length(r_eff) == length(x))
  } else {
    r_eff <- lapply(1:length(x), function(r) 1)
  }
  check.class <- !all(sapply(x, function(x) inherits(x, "causalWeights") |
                          inherits(x, "numeric")))
  if(check.class) stop("Members of list must be of class causalWeights or
                       a numeric vector")
  
  res <- mapply(function(w, r) {PSIS(w, r, ...)}, w = x, r = r_eff,
                SIMPLIFY = FALSE)
  
  return(res)
  
  
}


# setGeneric("PSIS", function(x, r_eff = NULL, ...) UseMethod("PSIS"))
# setMethod("PSIS", signature(x = "numeric"), PSIS.default)
# setMethod("PSIS", signature(x = "causalWeights"), PSIS.causalWeights)
# setMethod("PSIS", signature(x = "list"), PSIS.list)

PSIS_diag <- function(x, ...) UseMethod("PSIS_diag")

PSIS_diag.numeric <- function(x) {

  y   <- PSIS(x)
  res <- y$diagnostics
  return(res)
}

PSIS_diag.causalWeights <- function(x, r_eff = NULL, ...) {
  
  y   <- PSIS.causalWeights(x, r_eff = r_eff, ...)
  res <- list(w0 = y$w0$diagnostics,
              w1 = y$w1$diagnostics)
  
  return(res)
}

PSIS_diag.causalPSIS <- function(x) {


  res <- list(w0 = x$w0$diagnostics,
              w1 = x$w1$diagnostics)
  
  return(res)
}

PSIS_diag.list <- function(x, r_eff = NULL) {
  
  if(!is.null(r_eff)) {
    if(length(r_eff) == 1) r_eff <- rep(r_eff, length(x))
    stopifnot(length(r_eff) == length(x))
  } else {
    r_eff <- lapply(1:length(x), function(r) 1)
  }
  
  res <- mapply(function(w, r) {PSIS_diag(x=w, r_eff = r, ...)}, w = x, r = r_eff,
                SIMPLIFY = FALSE)

  # res <- lapply(x, function(w) PSIS_diag(w))

  return(res)

}




