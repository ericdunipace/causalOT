#' @include DataSimClass.R
setOldClass(c("DataSim","R6"))
setOldClass(c("Hainmueller", "DataSim","R6"))

check_dataHolder <- function(object) {
  errors <- character()
  if(!is.matrix(object@x)){
    msg <- paste("x is not a matrix! ")
    errors <- c(errors, msg)
  } else {
    n_x = nrow(object@x)
  }
  if(!is.integer(object@z)) {
    msg <- paste("z is not an integer! ")
    errors <- c(errors, msg)
  } else {
    n_z = length(object@z)
  }
  if(!all(object@z %in% c(0,1))) {
    msg <- paste("z must be in {0,1}! ")
    errors <- c(errors, msg)
  }
  if(!is.numeric(object@y)) {
    msg <- paste("y is not a numeric vector! ")
    errors <- c(errors, msg)
  } else {
    n_y = length(object@y)
  }
  if(!is.numeric(object@weights)) {
    msg <- paste("sample weights are not a numeric vector! ")
    errors <- c(errors, msg)
  } else {
    n_w = length(object@weights)
  }
  
  n <- (object@n1 + object@n0)
  if(n_z != n) {
    msg <- paste("Length of z not equal total sample size! ")
    errors <- c(errors, msg)
  }
  
  if(n_z != n_x) {
    msg <- paste("Length of z not equal number of rows of x! ")
    errors <- c(errors, msg)
  }
  
  if(n_x !=  n){
    msg <- paste("Number of rows of x not equal total sample size! ")
    errors <- c(errors, msg)
  }
  
  if(!all(is.na(object@y)) && n_y !=  n){
    msg <- paste("Length of y not equal total sample size! ")
    errors <- c(errors, msg)
  }
  
  if(n_w !=  n){
    msg <- paste("Length of sample weights not equal total sample size! ")
    errors <- c(errors, msg)
  }
  
  if(sum(object@z) !=  object@n1){
    msg <- paste("Total number of treated not equal n1! ")
    errors <- c(errors, msg)
  }
  
  if(sum(1- object@z) !=  object@n0){
    msg <- paste("Total number of control not equal n0! ")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}

#' dataHolder-class
#' 
#' @slot x matrix. A matrix of confounders.
#' @slot z integer. The treatment indicator, \eqn{z_i \in \{0,1\}}.
#' @slot y numeric. The outcome data.
#' @slot n0 integer. The number of observations where z==0
#' @slot n1 integer. The number of observations where z==1
#' @slot weights numeric. The empirical distribution of the full sample.
#' @keywords internal
#' @rdname dataHolder-class
setClass("dataHolder",
         slots = c(x = "matrix",
                   z = "integer",
                   y = "numeric",
                   n0 = "integer",
                   n1 = "integer",
                   weights = "numeric"),
         validity = check_dataHolder
         )

# functions to pull out slots and subset x
setGeneric("get_z", function(object) standardGeneric("get_z"))
setMethod("get_z", signature(object = "dataHolder"), function(object) {
  object@z
})
setGeneric("get_x", function(object) standardGeneric("get_x"))
setMethod("get_x", signature(object = "dataHolder"), function(object) {
  object@x
})
setGeneric("get_x0", function(object) standardGeneric("get_x0"))
setMethod("get_x0", signature(object = "dataHolder"), function(object) {
  object@x[object@z == 0,,drop = FALSE]
})
setGeneric("get_x1", function(object) standardGeneric("get_x1"))
setMethod("get_x1", signature(object = "dataHolder"), function(object) {
  object@x[object@z == 1,,drop = FALSE]
})
setGeneric("get_y", function(object) standardGeneric("get_y"))
setMethod("get_y", signature(object = "dataHolder"), function(object) {
  object@y
})
setGeneric("update_x", function(object, x) standardGeneric("update_x"))
setMethod("update_x", signature(object = "dataHolder", x = "ANY"), function(object, x) {
  x <- as.matrix(x)
  n_new <- nrow(x)
  if(object@n0+object@n1 != n_new) stop("new x doesn't have the same number of rows as the data")
  object@x <- x
  return(object)
})
setGeneric("update_y", function(object, y) standardGeneric("update_y"))
setMethod("update_y", signature(object = "dataHolder", y = "ANY"), function(object, y) {
  y <- as.numeric(y)
  n_new <- length(y)
  if(object@n0+object@n1 != n_new) stop("new y isn't the same length as data")
  object@y <- y
  return(object)
})
setGeneric("get_y1", function(object) standardGeneric("get_y1"))
setMethod("get_y1", signature(object = "dataHolder"), function(object) {
  object@y[object@z == 1]
})
setGeneric("get_y0", function(object) standardGeneric("get_y0"))
setMethod("get_y0", signature(object = "dataHolder"), function(object) {
  object@y[object@z == 0]
})
setGeneric("get_n", function(object) standardGeneric("get_n"))
setMethod("get_n", signature(object = "dataHolder"), function(object) {
  object@n1 + object@n0
})
setGeneric("get_n0", function(object) standardGeneric("get_n0"))
setMethod("get_n0", signature(object = "dataHolder"), function(object) {
  object@n0
})
setGeneric("get_n1", function(object) standardGeneric("get_n1"))
setMethod("get_n1", signature(object = "dataHolder"), function(object) {
  object@n1
})
setGeneric("get_w", function(object) standardGeneric("get_w"))
setMethod("get_w", signature(object = "dataHolder"), function(object) {
  object@weights
})
setGeneric("get_w0", function(object) standardGeneric("get_w0"))
setMethod("get_w0", signature(object = "dataHolder"), 
function(object) {
  return(renormalize(object@weights[object@z==0]))
})
setGeneric("get_w1", function(object) standardGeneric("get_w1"))
setMethod("get_w1", signature(object = "dataHolder"), 
          function(object) {
            return(renormalize(object@weights[object@z==1]))
          })

# constructor function

dataHolder <- function(x,z,y = NA_real_,weights=NA_real_) UseMethod("dataHolder")

#' dataHolder
#' 
#' @details Creates an object used internally by the `causalOT` package for data management. 
#'
#' @param x the covariate data. Can be a matrix, an object of class `dataHolder` or a [causalOT::DataSim] object. The latter two object types won't need arguments `z` or `y`.
#' @param z the treatment indicator
#' @param y the outcome data
#' @param weights the empirical distribution of the sample
#'
#' @return Returns an object of class dataHolder with slots
#' \itemize{
#' \item `x` matrix. A matrix of confounders.
#' \item `z` integer. The treatment indicator, \eqn{z_i \in \{0,1\}}.
#' \item `y` numeric. The outcome data.
#' \item `n0` integer. The number of observations where z==0
#' \item `n1` integer. The number of observations where z==1
#' \item `weights` numeric. The empirical distribution of the full sample.
#' }
#' @export
#'
#' @examples
#' x <- matrix(0, 100, 10)
#' z <- stats::rbinom(100, 1, 0.5)
#' 
#' # don't need to provide outcome
#' # function will assume each observation gets equal mass
#' dataHolder(x = x, z = z) 
setGeneric("dataHolder", function(x,z,y = NA_real_,weights=NA_real_) standardGeneric("dataHolder"))

#' @rdname dataHolder-methods
#' @keywords internal
setMethod("dataHolder", signature(x = "dataHolder", z = "ANY", y = "ANY",weights = "ANY" ), function(x, z = NA_integer_, y = NA_real_) {
  return(x)
})


# x is matrix, z is numeric, no balance subset
#' dataHolder-methods
#' @rdname dataHolder-methods
#' @keywords internal
setMethod("dataHolder", signature(x = "matrix", z = "ANY", y = "ANY", weights = "ANY"), function(x, z, y = NA_real_, weights = NA_real_) {
  z <- as.integer(z)
  if(length(z) > 0  && max(z) > 1) {
    z <- as.integer(as.numeric(factor(z)) - 1L)
  }
  n1 <- as.integer(sum(z))
  n0 <- as.integer(sum(1-z))
  y <- as.numeric(y)
  
  if(length(y) == 0) y <- NA_real_
  
  
  if(any(is.na(weights)) || is.null(weights)) {
    n <- n1 + n0
    weights <- rep(1.0/n,n)
  }
  weights <- renormalize(as.numeric(weights))
  
  dH <- methods::new("dataHolder", x = x, z = z, y = y, n0 = n0, n1 = n1, weights = weights )
  return(dH)
})

# x is DataSim, z is numeric, no balance subset
# setOldClass(c("DataSim","R6"))
# setOldClass(c("Hainmueller", "DataSim","R6"))

#' dataHolder-methods
#' @rdname dataHolder-methods
#' @keywords internal
dataHolder.DataSim <- function(x, z, y = NA_real_, weights = NA_real_) {
  
  ds <- x
  x_mat <- as.matrix(ds$get_x())
  z <- as.integer(ds$get_z())
  if(length(z) > 0 && max(z) > 2) {
    z <- as.integer(as.numeric(factor(z)) - 1L)
  }
  ns <- ds$get_n()
  n1 <- as.integer(ns[2])
  n0 <- as.integer(ns[1])
  y <- as.numeric(ds$get_y())
  
  if(length(y) == 0) y <- NA_real_
  
  if(any(is.na(weights)) || is.null(weights)) {
    n <- n1 + n0
    weights <- rep(1.0/n,n)
  }
  weights <- renormalize(as.numeric(weights))
  
  dH <- methods::new("dataHolder", x = x_mat, z = z, y = y, n0 = n0, n1 = n1, weights = weights )
  return(dH)
}

# setMethod("dataHolder", signature(x = "DataSim", z = "ANY", y = "ANY", weights = "ANY"), dataHolder.DataSim)

#' dataHolder-methods
#' @rdname dataHolder-methods
#' @keywords internal
setMethod("dataHolder", signature(x = "ANY", z = "ANY", y = "ANY", weights = "ANY" ), function(x, z = NA_integer_, y = NA_real_, weights = NA_real_) {
  if ( inherits(x, "DataSim") ) return(dataHolder.DataSim(x, z, y))
  
  stop("must supply one of the following combinations: (x =matrix, z = integer), (x = result of df2dataHolder, z = ANY), (x = dataHolder, z = ANY), or (x = DataSim, z = ANY)")
})


# method to go from data.frame and a formula to a dataHolder object
#' df2dataHolder
#' 
#' @description Function to turn a data.frame into a dataHolder object.
#'
#' @param treatment.formula a formula specifying the treatment indicator and covariates. Required.
#' @param outcome.formula an optional formula specifying the outcome function.
#' @param data a data.frame with the data
#' @param weights optional vector of sampling weights for the data
#'
#' @return Returns an object of class [causalOT::dataHolder()]
#' 
#' @details This will take the formulas specified and transform that data.frame into a dataHolder object that is used internally by the causalOT package. Take care if you do not specify an outcome formula that you do not include the outcome in the data.frame. If you are not careful, the function may include the outcome as a covariate, which is not kosher in causal inference during the design phase. 
#' 
#' If both outcome.formula and treatment.formula are specified, it will assume you are in the design phase, and create a combined covariate matrix to balance on the assumed treatment and outcome models. 
#' 
#' If you are in the outcome phase of estimation, you can just provide a dummy formula for the treatment.formula like "z ~ 0" just so the function can identify the treatment indicator appropriately in the data creation phase.
#' 
#' @export
#'
#' @examples
#' 
#' set.seed(20348)
#' n <- 15
#' d <- 3
#' x <- matrix(stats::rnorm(n*d), n, d)
#' z <- rbinom(n, 1, prob = 0.5)
#' y <- rnorm(n)
#' weights <- rep(1/n,n)
#' df <- data.frame(x, z, y)
#' dh <- df2dataHolder(
#'   treatment.formula = "z ~ .",
#'   outcome.formula = "y ~ ." ,
#'   data = df,
#'   weights = weights)
setGeneric("df2dataHolder", function(treatment.formula, outcome.formula = NA_character_, data, weights = NA_real_) standardGeneric("df2dataHolder"))

#' df2dataHolder-methods
#' @rdname df2dataHolder-methods
#' @keywords internal
setMethod("df2dataHolder", signature(treatment.formula = "ANY", outcome.formula = "ANY", data = "data.frame",weights = "ANY"), function(treatment.formula = NA_character_, outcome.formula = NA_character_, data, weights = NA_real_) {
  # browser()
  tf <- !missing(treatment.formula) && (inherits(treatment.formula, "formula") || !is.na(treatment.formula))
  of <- !missing(outcome.formula) && (inherits(outcome.formula, "formula") || !is.na(outcome.formula))
  
  if (of) {
    outcome.formula <- stats::formula(outcome.formula)
    mf.y <- stats::model.frame(outcome.formula, data)
    y <- stats::model.response(mf.y, "numeric")
    mt.y <- attr(mf.y, "terms")
    y.name <- attr(mf.y,"names")[1]
    attr(mt.y,"intercept") <- 0
    # x <- stats::model.matrix(mt.y, mf.y)
  } else {
    y <- NA_real_
  }
  if (tf) {
    treatment.formula <- stats::formula(treatment.formula)
    mf.z <- stats::model.frame(formula = treatment.formula, data = data)
    z <- stats::model.response(mf.z, "any")
    if (is.factor(z)) {
      if(length(levels(z)) > 2) stop("Treatment indicator must be binary.")
      vals.z <- levels(z)
      z <- as.integer(z != vals.z[1L])
    }
    if (!all(z %in% c(0L, 1L))) {
      vals.z <- sort(unique(z))
      z <- as.integer( z != vals.z[1L] )
    } else {
      vals.z <- c(0L, 1L)
    }
    z  <- as.integer(z)
    mt.z <- attr(mf.z, "terms")
    z.name <- attr(mf.z,"names")[1]
    
    if(of) {
      tf2 <- stats::update(mt.z, paste0(". ~ . -", y.name))
      of2 <- stats::update(mt.y, paste0("z ~ . - ", z.name))
      tf_new <- stats::update(tf2, paste0(". ~. + ", as.character(of2[3])))
      mf.z <- stats::model.frame(tf_new, data)
      
      mt.z <- attr(mf.z, "terms")
    }
    attr(mt.z,"intercept") <- 0
    x <- stats::model.matrix(mt.z, mf.z)
    
  } else {
    stop("Must specify treatment formula to identify treatment indicator 'z' in data, even if it's just `z ~ 0'")
  }
  
  
  n1 <- as.integer(sum(z))
  n0 <- as.integer(sum(1-z))
  
  if(any(is.na(weights)) ) {
    n <- n1 + n0
    weights <- rep(1.0/n, n)
  }
  dH <- methods::new("dataHolder", x = x, z = z, y = y, n0 = n0, n1 = n1, weights = weights)
  
  # setup terms objects in the attributes
  attr(mt.z, "tx_vals") <- vals.z
  if(of) {
    termsDH <- list(treatment = mt.z,
                    outcome = mt.y)
    
  } else {
    termsDH <- list(treatment = mt.z)
  }
  
  
  attributes(dH) <- c(attributes(dH), list(terms = termsDH))
  return(dH)
})

#' dataHolder-methods
#'
#' @param x dataHolder object constructed from a formula
#' @param ...  Not used at this time
#'
#' @return a list with the formula terms for treatment and, if present, outcome formulae.
#' @method terms dataHolder
#' @keywords internal
#' @rdname dataHolder-methods
terms.dataHolder <- function(x, ...) {
  attr(x, "terms")
}


# setMethod("terms", signature(x = "dataHolder"),
#     function(x, ...) {
#       attr(x, "terms")
#     }
# )



#' print.dataHolder
#'
#' @param x dataHolder object
#' @param ... Not used
#'
#' @return NULL
#' @method print dataHolder
#' @export 
print.dataHolder <- function(x, ...) {
  utils::str(x)
}


#' @keywords internal
setMethod("show", signature(object = "dataHolder"),
      function(object) {
        attributes(object)$terms <- NULL
        utils::str(object)
      })


