setClass("likelihoodMethods",
         slots = c(
           data = "dataHolder",
           estimand = "character",
           method = "character",
           solver = "function",
           solver.options = "list"
         )
)

likelihoodMethods <- function(data, estimand, method, options = NULL) {
  
  stopifnot(inherits(data, "dataHolder"))
  
  glm_methods <- c("Logit", "Logistic", "Probit")
  allowed_methods <- c(glm_methods, "CBPS")
  
  estimand <- match.arg(estimand, c("ATE", "ATT","ATC"))
  method <- match.arg(method, allowed_methods)
  
  stopifnot(estimand %in% c("ATE", "ATT","ATC"))
  stopifnot(method %in% allowed_methods)
  if(method == "Logit") method <- "Logistic"
  
  optimizer <- if(method %in% glm_methods) {
    glm_optimizer
  } else if (method == "CBPS") {
    cbps_optimizer
  } else {
    stop("Method not found")
  }
  
  if(is.null(options)) options <- list(NULL)
  if(!is.list(options)) stop("options must be a list")
  
  prob <- methods::new("likelihoodMethods",
      data = data,
      estimand = estimand,
      method = method,
      solver = optimizer,
      solver.options = options
  )
  return(prob)
}

glm_optimizer <- function(object) {
  stopifnot(inherits(object, "likelihoodMethods"))
  dH <- object@data
  method <- object@method
  estimand <- object@estimand
  solver.options <- object@solver.options
  
  stopifnot(method %in% c("Logistic", "Probit"))
  
  x <- get_x(dH)
  z <- get_z(dH)

  n1 <- get_n1(dH)
  n0 <- get_n0(dH)
  n  <- get_n(dH)
  
  weights <- get_w(dH) * n
  
  fam <- switch(method,
                "Logistic" = binomial(link = "logit"),
                "Probit" = binomial(link = "probit"))
  formula <- "z ~ ."
  
  mod <- do.call("glm", c(list(formula = formula, family = fam, 
                          data = data.frame(z = z, x),
                          weights = weights),
                          solver.options))
  predicted_prob <- stats::predict(mod, type = "response")

  output <- list(w0 = NULL, w1 = NULL,
                 estimand = estimand, method = method
                 )
  
  if (estimand == "ATT") {
    output$w1 <- rep(1/n1,n1)
    output$w0 <- predicted_prob[z == 0]/(1 - predicted_prob[z == 0]) * 1/n1
  } else if (estimand == "ATC") {
    output$w1 <- (1 - predicted_prob[z == 1])/predicted_prob[z == 1] * 1/n0
    output$w0 <- rep(1/n0,n0)
  } else if (estimand == "ATE") {
    output$w1 <- 1/predicted_prob[z == 1] * 1/n
    output$w0 <- 1/(1 - predicted_prob[z == 0]) * 1/n
  }

  return(output)
}

cbps_optimizer <- function(object) {
  stopifnot(inherits(object, "likelihoodMethods"))
  dH <- object@data
  method <- object@method
  estimand <- object@estimand
  solver.options <- object@solver.options
  
  stopifnot(method == "CBPS")
  
  x <- get_x(dH)
  z <- get_z(dH)
  
  n1 <- get_n1(dH)
  n0 <- get_n0(dH)
  n  <- get_n(dH)
  
  weights <- get_w(dH)
  
  ATT.flag <- switch(estimand,
                      "ATT" = 1,
                      "ATC" = 2,
                      "ATE" = 0)
  
  formula <- "z ~."
  
  cbp.dat <- data.frame(z = z, x)
  
  fit <- do.call(CBPS::CBPS, c(
                  list(formula = formula, data = cbp.dat,
              ATT = ATT.flag, sample.weights = weights),
              solver.options))
  
  predicted_prob <- fit$weights
  
  
  output <- list(w0 = predicted_prob[z == 0], 
                  w1 = predicted_prob[z == 1], 
                  estimand = estimand, method = "CBPS")
  
  return(output)

 }

setMethod("print", signature(x = "likelihoodMethods"), 
function(x,...) {
  utils::str(x)
})


# cot_solve method --------------------------------------------------------

#' cot_solve method for likelihoodMethods
#'
#' @param object likelihoodMethods. 
#' 
#' @include calc_weight.R
#'
#' @return object of class [causalWeights][causalOT::causalWeights-class]
#' @keywords internal
setMethod("cot_solve", signature(object = "likelihoodMethods"),
          function(object) {
            
            fit <- object@solver(object)
            
            cw <- methods::new("causalWeights", 
                               w0 = fit$w0,
                               w1 = fit$w1,
                               estimand = fit$estimand,
                               method = fit$method,
                               penalty = list(NULL),
                               info = list(NULL),
                               data = object@data,
                               call = call("calc_weight")
            )
            
            if(!inherits(cw, "causalWeights")) stop("Solver didn't return object of class causalWeights!")
            return(cw)
          }
)
