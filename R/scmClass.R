#### SCM ####
SCM <- R6::R6Class("SCM",
                   public = list(
                     source = "matrix",
                     target = "matrix",
                     a = "vector",
                     b = "vector",
                     A = "matrix",
                     bf = "R6",
                     gridInit = function(...) {
                       return(c(NA_real_))
                     },
                     solve = function(penalty = NULL, w = NULL) {
                       
                       n_t <- nrow(self$target)
                       n_s <- nrow(self$source)
                       d   <- ncol(self$source)
                       
                       w_list <- vector("list", n_t)
                       res <- vector("list", 1)
                       addl_0 <- rep(0,d)
                       w_init <- c(rep(1/n_s,n_s), addl_0)
                       
                       # update BF penalty function
                       if(inherits(self$bf, "balanceFunction") && !is.null(penalty)) {
                         if(is.list(penalty)) penalty <- penalty[[1]]
                         if(is.numeric(penalty) && penalty >= 0) {
                           self$bf$delta <- as.numeric(penalty)
                         } else {
                           warning("provided penalty for balanceFunction method is not a number > 0. Not used!")
                         }
                       }
                       
                       #update target
                       t_idx <- private$target_idx
                       l <- private$solver$GetData(element = "l")
                       u <- private$solver$GetData(element = "u")
                       
                       for (j in 1:n_t) {
                         l[t_idx] <- u[t_idx] <- c(self$target[j,])
                         private$solver$Update(l = l, u = u)
                         res[[1]] <- osqp_R6_solve(private$solver, self$bf$delta, 
                                                   self$bf$delta_idx, w_init, normalize = FALSE)
                         res[[1]][res[[1]] < 0] <- 0
                         w_list[[j]] <- renormalize(res[[1]][1:n_s]) * self$b[j]
                       }
                       return(Reduce("+", w_list))
                     },
                     initialize = function(source, target,
                                           a = NULL, b = NULL,
                                           options = list()) {
                       # browser()
                       if(!inherits(options, "scmOptions")) {
                         if(!is.list(options)) stop("options must be a list or output of scmOptions function")
                         options <- do.call(scmOptions, options)
                       }
                       
                       self$source <- as.matrix(source)
                       self$target <- as.matrix(target)
                       self$a <- check_weights(a, self$source)
                       self$b <- check_weights(b, self$target)
                       
                       n <- nrow(self$source)
                       d <- ncol(self$source)
                       if (!is.null(options$balance.formula)) {
                         if (is.null(colnames(source))) colnames(source) <- paste0("X", 1:ncol(source))
                         if (is.null(colnames(target))) colnames(target) <- colnames(source)
                         tbf <- terms(formula(options$balance.formula))
                         attr(tbf, "intercept") <- 0
                         source.bf <- model.matrix(tbf, data.frame(source))
                         target.bf <- model.matrix(tbf, data.frame(target))
                         self$bf <- balanceFunction$new(source = source.bf, 
                                                        target = target.bf,
                                                        a = as.numeric(self$ot$a),
                                                        b = as.numeric(self$ot$b),
                                                        delta = options$delta)
                         # self$runbf <- TRUE
                         k <- ncol(target.bf)
                       } else {
                         # self$runbf <- FALSE
                         self$bf <- list(NULL)
                         k <- 0
                       }
                       
                       nvars <- n + d
                       
                       # q <- c(self$source %*% c(self$target[1,]), rep(0, d))
                       q <- NULL
                       P <- Matrix::sparseMatrix(i = n + 1:d, j = n + 1:d, x = 1)
                       
                       A_quad <- cbind(Matrix::Matrix(t(self$source)),
                                       Matrix::Diagonal(d, x = 1)
                       )
                       u_quad <- l_quad <- c(self$target[1,])
                       
                       A_bounds <- rbind(cbind(Matrix::Diagonal(n, x = 1),
                                               Matrix::Matrix(matrix(0, n, d))
                       ),
                       Matrix::sparseMatrix(i = rep(1,n), j = 1:n, x = 1, dims = c(1, nvars)))
                       l_bounds <- c(rep(0,   n), 1)
                       u_bounds <- c(rep(Inf, n), 1)
                       
                       if(k > 0) {
                         A_bf <- cbind(
                           Matrix::Matrix(
                             t(self$bf$A)
                           ),
                           Matrix::Matrix(matrix(0, k, d)))
                         l_bf <- rep(-self$bf$delta, k)
                         u_bf <- rep(self$bf$delta, k)
                         self$bf$delta_idx <- 1:length(l_bf)
                       } else {
                         A_bf <- l_bf <- u_bf <- NULL
                       }
                       
                       l <- c(l_bf, l_bounds, l_quad)
                       u <- c(u_bf, u_bounds, u_quad)
                       A <- rbind(A_bf, A_bounds, A_quad)
                       
                       private$target_idx <- (k + n + 2):length(l)
                       
                       private$solver <- osqp::osqp(P = P, q = q,
                                                    A = A, l = l, u = u,
                                                    pars = options$solver.options)
                       
                     }
                   ),
                   private = list(
                     target_idx = "numeric",
                     solver = "R6"
                   )
)


# @param lambda Penalty parameter on the weights. Not currently used but here because a plan is to add it.
# @param delta The constraint parameter for the balancing functions. Not currently used.
# @param grid.length The number of penalty parameters to try. Not currently used.
# @param nboot The number of bootstrap samples. Not currently used.
# @param balance.formula The formula that denotes the covariate functions to balance. Not currently used.

#' Options for the SCM Method
#'
#' @param ... Arguments passed to the [osqpSettings()][osqp::osqpSettings()] function which solves the problem.
#'
#' @return A list with arguments to pass to [osqpSettings()][osqp::osqpSettings()]
#' @export
#' 
#' @details Options for the solver used in the optimization of the Synthetic Control Method of Abadie and Gardeazabal (2003).
#'
#' @examples
#' opts <- scmOptions()
scmOptions <- function(
    # lambda = NULL,
  #                    delta = NULL,
  #                    grid.length = 7L,
  #                    nboot = 1000L,
  # balance.formula = NULL,
  ...) { # dots are the osqp args
  mc <- match.call()
  used.args <- as.list(mc)[-1]
  # browser()
  
  grid.length <- 7L
  nboot <- 1000L
  delta <- NULL
  balance.formula <- NULL
  gsOpts <- gridSearchOptions(nboot = nboot, grid.length = grid.length)
  
  nboot <- gsOpts$nboot
  grid.length <- gsOpts$grid.length
  
  output <- list()
  # lambda currently not used but may consider in future with following code
  # if(arg_not_used(lambda)) {
  #   output["lambda"] <- list(NULL)
  # } else {
  #   if(any(lambda < 0)) stop("lambda must be >= 0")
  #   output$lambda <- sort(lambda, decreasing = TRUE)
  # }
  output$lambda <- NULL
  
  if(arg_not_used(delta)) {
    output["delta"]<- list(NULL)
  } else {
    if(any(delta < 0)) stop("delta must be >= 0")
    output$delta <- sort(delta, decreasing = TRUE)
  }
  
  # also only L2 penalty at this time but may consider in future
  # if ( arg_not_used(penalty) ) {
  #   output$penalty <- "entropy"
  # } else {
  #   output$penalty <- match.arg( penalty, c("entropy", "L2") )
  # }
  
  if( arg_not_used(balance.formula) ) {
    output["balance.formula"] <- list(NULL)
  } else {
    balance.formula <- as.character(balance.formula)
    bf_split <- strsplit(balance.formula, "~")
    output$balance.formula <- paste0("~ 0 +", bf_split[[1]][2])
  }
  
  # only for delta at this time
  if ( arg_not_used(grid.length) ) {
    output$grid.length <- 7L
  } else {
    output$grid.length <- as.integer(grid.length)
    if(grid.length <= 0) stop("grid.length must be greater than 0")
  }
  
  if (!is.null(nboot)) {
    output$nboot <- as.integer(nboot)
  } else {
    output$nboot <- 1000L
  }
  
  # if (!is.null(output$lambda) && !is.null(output$delta)) {
  #   output["grid.length"] <- list(NULL)
  # }
  if ( !is.null(output$delta)) {
    output["grid.length"] <- list(NULL)
  }
  
  # if (output$penalty == "L2") {
  output$solver.options <- list(...)[...names() %in% methods::formalArgs(osqp::osqpSettings)]
  # }
  
  class(output) <- "scmOptions"
  return(output)
}