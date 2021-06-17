prep_data.data.frame <- function(data,...) {
  
  #create data.frame
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  cov_bal <- dots$balance.covariates
  outcome <- dots$outcome
  
  if (tx_ind %in% cov_bal) stop("treatment.indicator in balance.covariates!")
  if (!is.null(outcome) ) {
    if (tx_ind %in% outcome) stop("treatment.indicator is outcome!")
    if (outcome %in% cov_bal) stop("outcome in balance.covariates!")
  }

  if (is.null(cov_bal)) {
    stop("must specify covariates for balance 'balance.covariates' either by name or column number")
  }
  
  if (is.null(tx_ind)) {
    stop("must specify treatment indicator 'treatment.indicator' either by name or column number")
  }
  
  # if(is.null(outcome)) {
  #   stop("must specify outcome 'outcome' either by name or column number")
  # }
  
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
  if(!is.null(outcome)) {
    y.var <- if(is.character(outcome)) {
      match(outcome, colnames(data))
    } else {
      outcome
    }
    df <- data.frame(y = data[[y.var]], data[x.vars])
    attr(df, "outcome") <- "y"
    attr(df, "balance.covariates") <- colnames(df)[2:ncol(df)]
  } else {
    df <- data.frame(data[x.vars])
    attr(df, "balance.covariates") <- colnames(df)
  }
  
  
  
  # df <- data[c(y.var, x.vars)]
  z <- as.integer(data[[tx.var]])
  attr(z, "treatment.indicator") <- "z"
  return(list(df = df, z = z))
}

prep_data.matrix <- function(data,...) {
  
  #create data.frame
  dots <- list(...)
  tx_ind <- dots$treatment.indicator
  cov_bal <- dots$balance.covariates
  outcome <- dots$outcome
  
  if (tx_ind %in% cov_bal) stop("treatment.indicator in balance.covariates!")
  if (!is.null(outcome) ) {
    if (tx_ind %in% outcome) stop("treatment.indicator is outcome!")
    if (outcome %in% cov_bal) stop("outcome in balance.covariates!")
  }
  if (tx_ind %in% cov_bal) stop("treatment.indicator in balance.covariates!")
  
  
  if(is.null(cov_bal)) {
    stop("must specify covariates for balance 'balance.covariates' either by name or column number")
  }
  
  if(is.null(tx_ind)) {
    stop("must specify treatment indicator 'treatment.indicator' either by name or column number")
  }
  
  # if(is.null(outcome)) {
  #   stop("must specify outcome 'outcome' either by name or column number")
  # }
  
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
  if (!is.null(outcome)) {
    y.var <- if(is.character(outcome)) {
      match(outcome, colnames(data))
    } else {
      outcome
    }
    df <- data.frame(y = data[,y.var], data[,x.vars])
    attr(df, "outcome") <- "y"
    attr(df, "balance.covariates") <- colnames(df)[2:ncol(df)]
  } else {
    df <- data.frame(data[,x.vars])
    attr(df, "balance.covariates") <- colnames(df)
  }
  
  # df <- data[c(y.var, x.vars)]
  z <- as.integer(data[,tx.var])
  attr(z, "treatment.indicator") <- "z"
  return(list(df = df, z = z))
}

prep_data.list <- function(data, ...) {
  #create data.frame
  df <- data.frame(y = data$y, data$x)
  attr(df, "outcome") <- "y"
  attr(df, "balance.covariates") <- colnames(df)[2:ncol(df)]
  z <- as.integer(data$z)
  attr(z, "treatment.indicator") <- "z"
  return(list(df = df, z = z))
}

prep_data.DataSim <- function(data, ...) {
  #create data.frame
  df <- data.frame(y = data$get_y(), data$get_x())
  z <- data$get_z()
  attr(df, "outcome") <- "y"
  attr(df, "balance.covariates") <- colnames(df)[2:ncol(df)]
  attr(z, "treatment.indicator") <- "z"
  return(list(df = df, z = z))
}

#' prep_data function to prepare data for weighting and outcome estimation
#'
#' @param data Either a matrix, a data.frame, or a DataSim class
#' @param ... Additional arguments are necessary if a matrix or data.frame is given. These are:
#' \itemize{
#' \item{balance.covariates: }{ The names of the covariates to include for balancing}
#' \item{treatment.indicator: }{ The name of the covariate denoting the treatment indicator}
#' \item{outcome: }{ (Optional) The name of the covariate that denotes the outcome.}
#' }
#'
#' @return A list with slots "df" and "z". "df" is a data.frame with the covariates and outcome and "z" is the treatment indicator as an integer vector.
#'
#' @keywords internal
setGeneric("prep_data", function(data, ...) UseMethod("prep_data"))
setMethod("prep_data", signature(data = "Hainmueller"), prep_data.DataSim)
setMethod("prep_data", signature(data = "Sonabend2020"), prep_data.DataSim)
setMethod("prep_data", signature(data = "DataSim"), prep_data.DataSim)
setMethod("prep_data", signature(data = "data.frame"), prep_data.data.frame)
setMethod("prep_data", signature(data = "matrix"), prep_data.matrix)
setMethod("prep_data", signature(data = "list"), prep_data.list)