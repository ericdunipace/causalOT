extract_x.DataSim <- function(data, ...) {
  x1 <- data$get_x1()
  x0 <- data$get_x0()
  return(list(x0 = x0, x1 = x1))
}

extract_x.data.frame <- function(data, ...) {
  #create data.frame
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  cov_bal <- dots$balance.covariates
  
  if(is.null(cov_bal)) {
    stop("must specify covariates for balance 'balance.covariates' either by name or column number")
  }
  
  if(is.null(tx_ind)) {
    stop("must specify treatment indicator 'treatment.indicator' either by name or column number")
  }
  
  tx.var <- if(is.character(tx_ind)) {
    match(tx_ind, colnames(data))
  } else {
    tx_ind
  }
  x.vars <- if(is.character(cov_bal)) {
    match(cov_bal, colnames(data))
  } else {
    cov_bal
  }
  z <- as.integer(data[ , tx.var])
  
  x1 <- data[z == 1, x.vars]
  x0 <- data[z == 0, x.vars]
  return(list(x0 = x0, x1 = x1))
}

setGeneric("extract_x", function(data, ...){UseMethod("extract_x")})
setMethod("extract_x", "DataSim", extract_x.DataSim)
setMethod("extract_x", "data.frame", extract_x.data.frame)
setMethod("extract_x", "matrix", extract_x.data.frame)