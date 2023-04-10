# TODO add grid_solve method, update for lambda and warm start
# TODO remove dependency on lbfgsb3c in favor of torch??? could get GPU support...


# balance function super class
balanceFunction <- R6::R6Class("balanceFunction",
  public = list(
    source = "matrix",
    target = "matrix",
    target_vector = "vector",
    A = "matrix",
    a = "vector",
    b = "vector",
    n = "integer",
    m = "integer",
    k = "integer",
    delta = "numeric",
    delta_idx = "numeric",
    grid_length = "integer",
    initialize = function(source, target, a = NULL, b = NULL, delta) {
      # browser()
      source <- as.matrix(source)
      target <- as.matrix(target)
      
      a <- check_weights(a, source)
      b <- check_weights(b, target)
      
      self$a <- a
      self$b <- b
      self$n <- length(a)
      self$m <- length(b)
      
      # drop 0 var covariates
      sd_source  <- matrixStats::colWeightedSds(source, w = a*self$n)
      non_zero <- which(sd_source > 0)
      
      if(length(non_zero) == 0) stop("All source covariates have 0 variance.")
      
      # other SD
      
      sd_target  <- matrixStats::colWeightedSds(target, w = b*self$m)
      sd_total   <- matrixStats::colWeightedMeans(
        rbind(sd_source, sd_target),
        w = c(self$n, self$m)
      )
      
      self$source <- scale(source[,non_zero, drop = FALSE], center = FALSE, scale = sd_total[non_zero])
      self$target <- scale(target[,non_zero, drop = FALSE], center = FALSE, scale = sd_total[non_zero])
      
      self$target_vector <- c(matrixStats::colWeightedMeans(self$target, w = b))
      # setup scaled and center matrix A
      self$A <- scale(self$source, 
                             center = self$target_vector, 
                             scale = FALSE)
      self$k <- ncol(self$A)
      
      if ( missing(delta) || is.null(delta) || all( is.na(delta) ) ) {
        self$delta <- max(abs(matrixStats:::colWeightedMeans(self$A, self$a)))
      } else {
        self$delta <- as.numeric(delta)
      }
      return(invisible(self))
    }
  ),
)

balanceFunction$set("public", "evalBoot",
 function(a, b) {
  
  # get weighted means for each group    
  mean_s <- c(matrixStats::colWeightedMeans(self$source, a))
  mean_t <- c(matrixStats::colWeightedMeans(self$target, b))
  
  # average absolute differences
  return(mean(abs(mean_s - mean_t)))
 }
)

balanceFunction$set("public", "gridInit", 
        function(grid, length) {
          
          if (missing(grid) || is.null(grid) || all(is.na(grid)) ) {
            if(missing(length) || is.null(length) || all(is.na(length)) ) {
              stop("One of grid.length or delta values must be provided in options")
            }
            max_delta <- max(abs(matrixStats::colWeightedMeans(self$A,self$a)))
            delta <- seq(max_delta, 1e-4, length.out = length)
          } else if (is.numeric(grid)) {
            delta <- sort(grid, decreasing = TRUE)
            if(!all(delta >= 0)) stop("Balancing function constraints for COT must be numbers >= 0.")
          } else {
            stop("Grid values weren't provided but I can't initialize any. Report this bug!")
          }
          
          return(delta)
        }
)

bfMethods <- function(source, target, a, b, method, options) {
  
  method <- match.arg(method, balancing_function_methods())
  # set up problems
  prob <- if (method == "SBW") {
    SBW$new(source, target, a, b, options)
  } else if (method == "EntropyBW") {
    EntropyBW$new(source, target, a, b,options)
  } else {
    stop("method not found!")
  }
 
  return(prob)
}

SBW <- R6::R6Class("SBW",
                   inherit = balanceFunction,
  public = list(
    solve = "function",
    initialize = function(source, target, a = NULL, b = NULL, options = list()) {
      
      # browser()
      
      if(!inherits(options, "sbwOptions")) {
        if(!is.list(options)) {
          options <- as.list(options)
        }
        options <- do.call("sbwOptions", options)
      }
      delta <- options$delta
      super$initialize(source,target, a, b, delta)
      
      # get lengths  
      n <- self$n
      k <- self$k
      self$grid_length <- options$grid.length
      
      # set linear constraints
      l_bounds <- rep(0,n)
      u_bounds <- rep(Inf,n)
      A_bounds <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1)
      if (options$sumto1) { # right now always sum to 1
        l_sum_const <- u_sum_const <- 1
        A_sum_const <- Matrix::sparseMatrix(i = rep(1,n),
                                            j = 1:n,
                                            x = 1)
      } else {
        l_sum_const <- u_sum_const <- A_sum_const <-NULL
      }
      A_delta <- Matrix::Matrix(data = t(self$A), 
                                         sparse = TRUE)
      l_delta <- rep(-self$delta[1], nrow(A_delta))
      u_delta <- rep(self$delta[1], nrow(A_delta))
      
      # set final params
      P <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1)
      q <- rep(0, n)
      l <- c(l_bounds, l_sum_const, l_delta)
      u <- c(u_bounds, u_sum_const, u_delta)
      
      self$delta_idx <- (length(c(l_bounds, l_sum_const))+1):length(l)
      A <- rbind(A_bounds, A_sum_const, A_delta)
      
      private$osqp_args <- list(P = P, q = q, A = A, l = l, u = u,
                                pars = options$solver.options)
      # browser()
      private$solver <- do.call(osqp::osqp, private$osqp_args)
      self$delta <- delta[1]
      self$solve <- function(penalty, w = NULL) {
        if (penalty < 0) stop("Penalty must be greater than or equal to 0")
        
        self$delta <- as.numeric(penalty)
        # if (is.null(w) || (is.list(w) && length(w) == 0)) {
        #   w <- self$a
        # } else {
        #   w <- if(is.list(w)) {
        #     w[[1]]
        #   } else {
        #     as.numeric(w)
        #   }
        # }
        private$osqp_args$u[self$delta_idx] <- self$delta
        private$osqp_args$l[self$delta_idx] <- -self$delta
        private$solver <- do.call(osqp::osqp, private$osqp_args)
        
        tryCatch(osqp_R6_solve(private$solver, NULL, self$delta_idx, self$a),
                      error = function(e) {NULL})
      }
    }
  ),
  private = list(
    osqp_args = "list",
    solver = "function"
  )
)

#' Options for the SBW method
#'
#' @param delta A number or vector of tolerances for the balancing functions. Default is NULL which will use a grid search
#' @param grid.length The number of values to try in the grid search
#' @param nboot The number of bootstrap samples to run during the grid search.
#' @param ... Arguments passed on to [osqpSettings()][osqp::osqpSettings()]
#'
#' @return A list of class `sbwOptions` with slots
#' \itemize{
#'  \item `delta` Delta values to try
#'  \item `grid.length` The number of parameters to try
#'  \item `sumto1` Forced to be TRUE. Weights will always sum to 1.
#'  \item `nboot` Number of bootstrap samples
#'  \item `solver.options`A list with arguments passed to [osqpSettings()][osqp::osqpSettings()]
#' }
#' @export
#'
#' @details
#' # Function balancing
#' This method will balance  functions of the covariates within some tolerance, \eqn{\delta}. For these functions \eqn{B}, we will desire
#' \deqn{\frac{\sum_{i: Z_i = 0} w_i B(x_i) - \sum_{j: Z_j = 1} B(x_j)/n_1}{\sigma} \leq \delta}, where in this case we are targeting balance with the treatment group for the ATT. \eqn{\sigma} is the pooled standard deviation prior to balancing.
#' 
#' @examples 
#' opts <- sbwOptions(delta = 0.1)
sbwOptions <- function(
                   delta = NULL,
                   grid.length = 20L,
                   nboot = 1000L,
                   # verbose = FALSE,
                   ...) {
  output <- list()
  # browser()
  if(arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  if ( arg_not_used(grid.length) ) {
    output$grid.length <- 20L
  } else {
    output$grid.length <- as.integer(grid.length)
    if(grid.length <= 0) stop("grid.length must be greater than 0")
  }
  output$sumto1 <- TRUE
  
  
  if ( arg_not_used(nboot) ) {
    output$nboot <- 1000L
  } else {
    output$nboot <- as.integer(nboot)
  }
  
  output$solver.options <- list(...)[...names() %in% methods::formalArgs(osqp::osqpSettings)]
  class(output) <- "sbwOptions"
  return(output)
  
}

EntropyBW <- R6::R6Class("EntropyBW",
 inherit = balanceFunction,
 public = list(
   solve = function(penalty, w = NULL) {
     # init <- private$weight_to_coefficient(w)
     k <- self$k
     beta_full <- lbfgs3c_R6_solve(init = rep(0, 2 * k), options = self$solver.options, 
                           bounds = private$bounds,
                           objective = private$objective,
                           gradient = private$gradient, 
                           A = self$A, 
                           delta = penalty)
    
     beta <- beta_full[1:k] - beta_full[-c(1:k)]
     lp <- c(self$A %*% beta)
     w <- exp(lp - logSumExp(c(lp)))
     return(w)
   },
   solver.options = "list",
   initialize = function(source, target, a = NULL, b = NULL, options = list()) {
     # browser()
     
     if(!inherits(options, "entBWOptions")) {
       if(!is.list(options)) {
         options <- as.list(options)
       }
       options <- do.call("entBWOptions", options)
     }
     delta <- options$delta
     self$solver.options <- options$solver.options
     
     super$initialize(source,target, a, b, delta)
     
     self$grid_length <- options$grid.length
     k <- self$k
     self$delta <- delta
     private$objective <- entBW_obj_
     private$gradient  <- entBW_grad_
     private$bounds <- cbind(rep(0,   2 * k),
                             rep(Inf, 2 * k))
     
   }
   
 ),
 private = list(
   bounds = "matrix",
   objective = "function",
   gradient = "function",
   weight_to_coefficient = function(w = NULL) {
     if(is.null(w)) w <- self$a
     y <- log(w)
     beta <- coef(lm(y ~ self$A + 0))
     beta_pos <- beta_neg <- rep(0, length(beta))
     beta_pos[beta>0] <- beta[beta>0]
     beta_neg[beta<0] <- abs(beta[beta<0])
     return(c(beta_pos, beta_neg))
   }
 )
)

#' Options for the Entropy Balancing Weights
#'
#' @param delta A number or vector of tolerances for the balancing functions. Default is NULL which will use a grid search
#' @param grid.length The number of values to try in the grid search
#' @param nboot The number of bootstrap samples to run during the grid search.
#' @param ... Arguments passed on to [lbfgsb3c()][lbfgsb3c::lbfgsb3c()]
#'
#' @return A list of class `entBWOptions` with slots
#' \itemize{
#'  \item `delta` Delta values to try
#'  \item `grid.length` The number of parameters to try
#'  \item `nboot` Number of bootstrap samples
#'  \item `solver.options` A list of options passed to `[lbfgsb3c()][lbfgsb3c::lbfgsb3c()]
#' }
#' @export
#'
#' @details
#' # Function balancing
#' This method will balance  functions of the covariates within some tolerance, \eqn{\delta}. For these functions \eqn{B}, we will desire
#' \deqn{\frac{\sum_{i: Z_i = 0} w_i B(x_i) - \sum_{j: Z_j = 1} B(x_j)/n_1}{\sigma} \leq \delta}, where in this case we are targeting balance with the treatment group for the ATT. \eqn{\sigma} is the pooled standard deviation prior to balancing.
#' 
#' @examples 
#' opts <- entBWOptions(delta = 0.1)
entBWOptions <- function(
    delta = NULL,
    grid.length = 20L,
    nboot = 1000L,
    ...) {
  output <- list()
  
  if(arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  if ( arg_not_used(grid.length) ) {
    output$grid.length <- 20L
  } else {
    output$grid.length <- as.integer(grid.length)
    if(grid.length <= 0) stop("grid.length must be greater than 0")
  }
  
  if ( arg_not_used(nboot) ) {
    output$nboot <- 1000L
  } else {
    output$nboot <- as.integer(nboot)
  }
  
  output$solver.options <- lbfgs3c_control(...)
  
  class(output) <- "entBWOptions"
  return(output)
}


# for use inside of OTProblem classes internally
SBW_4_oop <- R6::R6Class("SBW",
 public = list(
   source = "torch_tensor",
   target = "torch_tensor",
   target_vector = "vector",
   target_scale = "torch_tensor",
   source_scale = "torch_tensor",
   n = "integer",
   d = "integer",
   delta_idx = "numeric",
   solve = function(penalty) {
     if (penalty < 0) stop("Penalty must be greater than or equal to 0")
     
     model <- private$solver
     
     if (!missing(penalty) && !is.null(penalty)) {
       l <- model$GetData(element = "l")
       u <- model$GetData(element = "u")
       u[self$delta_idx] <- penalty + self$target_vector
       l[self$delta_idx] <- -penalty + self$target_vector
       model$Update(l = l, u = u)
     }
     
     res <- model$Solve()
     
     if (res$info$status_val == -3 || res$info$status_val == -4) {
       return(NULL) 
     } else {
       return(res$x)
     }
     
   },
   initialize = function(source, target, prob.measure = TRUE,
                         osqp_opts) {
     
     if (!inherits(source, "torch_tensor")) {
       source <- torch::torch_tensor(as.matrix(source),
                                     dtype = torch::torch_double())
     }
     
     if (!inherits(target, "torch_tensor")) {
       target <- torch::torch_tensor(target,
                                     dtype = torch::torch_double())
     }
     
     stopifnot("source and target must have same number of columns" = ncol(source) == length(target))
     
     self$source <- source$detach()
     self$target <- target$detach()
     
     self$source_scale  <- self$source/self$source$std(1)
     self$target_scale  <- self$target/self$source$std(1)
     self$target_vector <- as.numeric(self$target_scale$to(device = "cpu"))
     
     self$d <- ncol(self$source)
     self$n <- nrow(self$source)
     
     # quadratic
     P <- Matrix::Diagonal(self$n, 1)
     
     # set linear constraints
     l_bounds <- rep(0,   self$n)
     u_bounds <- rep(Inf, self$n)
     A_bounds <- Matrix::Diagonal(self$n, 1)
     
     if (prob.measure) { # right now always sum to 1
       l_sum_const <- u_sum_const <- 1
       A_sum_const <- Matrix::sparseMatrix(i = rep(1,self$n),
                                           j = 1:self$n,
                                           x = 1)
     } else {
       l_sum_const <- u_sum_const <- A_sum_const <-NULL
     }
     A_delta <- Matrix::Matrix(data = t(as.matrix(self$source_scale$to(device = "cpu"))), 
                               sparse = TRUE)
     
     l_delta <- self$target_vector
     u_delta <- self$target_vector
     
     # set final params
     q <- rep(0, self$n)
     l <- c(l_bounds, l_sum_const, l_delta)
     u <- c(u_bounds, u_sum_const, u_delta)
     
     self$delta_idx <- (length(c(l_bounds, l_sum_const))+1):length(l)
     A <- rbind(A_bounds, A_sum_const, A_delta)
     private$solver <- osqp::osqp(P = P, q = q, A = A, l = l, u = u,
                                  pars = osqp_opts)
     
   }
 ),
 private = list(
   solver = "function"
 )
)

SBW_4_oop$set("public", "evalBoot",
                    function(a) {
                      
                      # get weighted means for each group    
                      # mean_s <- c(matrixStats::colWeightedMeans(self$source_scale, a))
                      mean_s <- self$source_scale$transpose(1,2)$matmul(a)
                      mean_t <- self$target_scale
                      
                      # average absolute differences
                      return(as.numeric(mean(abs(mean_s - mean_t))))
                    }
)