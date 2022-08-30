#### parent class ####

#' R6 Data Generating Parent Class
#' 
#' @rdname DataSim
#' 
#' @return An [R6][R6::R6Class] object
#' 
#' @details Can be used to make your own data simulation class. 
#' Should have the same slots listed in this class at a minimum, 
#' but you can add your own, of course. An easy way to do this is
#' to make your class inherit from this one. See the example.
#' 
#' @export
#' 
#' @examples 
#' MyClass <- R6::R6Class("MyClass", 
#' inherit = DataSim,
#' public = list(),
#' private = list())
DataSim <- R6::R6Class("DataSim",
        public = list(
                      #' @description
                      #' Gets the covariate data
                      get_x = function() { return( private$x)},
                      #' @description
                      #' Gets the outcome vector
                      get_y = function() { return( private$y)},
                      #' @description
                      #' Gets the treatment indicator
                      get_z = function() { return( private$z)},
                      #' @description
                      #' Gets the number of observations
                      get_n = function() {return(c("n0" = private$n0, "n1" = private$n1))},
                      #' @description
                      #' Gets the covariate data for the treated individuals
                      get_x1 = function() {
                        if(!is.character(private$x)) {
                          if(is.character(private$x1)) {
                            private$x1 <- private$x[private$z == 1,,drop=FALSE]
                          }
                          return(private$x1)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      #' @description
                      #' Gets the covaraiate data for the control individuals
                      get_x0 = function(){
                        if(!is.character(private$x)) {
                          if(is.character(private$x0)) private$x0 <- private$x[private$z == 0,,drop=FALSE]
                          return(private$x0)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      #' @description
                      #' Gets the dimensionality covariate data
                      get_p = function() {
                        return(private$p)
                      },
                      #' @description
                      #' Gets the individual treatment effects
                      get_tau = function() {
                        if(is.character(private$mu1) | is.character(private$mu0)) {
                          stop("Need to generate outcome data first")
                        }
                        return(private$mu1 - private$mu0)
                      },
                      #' @description
                      #' Generates the  data. Default is an empty function
                      gen_data = function(){NULL},
                      #' @description
                      #' Gets the optimal weights to get the correct expectation
                      #' @param estimand One of "ATT","ATC","ATE"
                      #' @param augment Should we use an augmented estimator? TRUE or FALSE.
                      #' @param solver One of "mosek" or "gurobi"
                      opt_weight = function(estimand = "ATE", augment = FALSE, solver = "mosek") {
                          if(estimand == "cATE") estimand <- "ATE"
                          estimand <- match.arg(estimand, choices = c("ATT","ATC","ATE"))
                          aug <- isTRUE(augment)
                          solver <- match.arg(solver, choices = c("mosek", "gurobi", "cplex"))
                          design0 <- cbind(1,private$x[private$z==0,,drop = FALSE])
                          design1 <- cbind(1,private$x[private$z==1,,drop = FALSE])
                          design <- rbind(design0,design1)
                          # the fitted model
                          m0 <- .lm.fit(design0, private$y[private$z==0,drop=FALSE])
                          m1 <- .lm.fit(design1, private$y[private$z==1,drop=FALSE])
                          
                          # the values of y or the residuals
                          f = switch(as.character(aug),
                                     "FALSE" =switch(estimand,
                                                     "ATT" = private$y[private$z==0,drop=FALSE],
                                                     "ATC" = private$y[private$z==1,drop=FALSE],
                                                     "ATE" = list(private$y[private$z==0,drop=FALSE],
                                                                  private$y[private$z==1,drop=FALSE])),
                                     "TRUE" = switch(estimand,
                                                     "ATT" = m0$residuals,
                                                     "ATC" = m1$residuals,
                                                     "ATE" = list(m0$residuals,
                                                                  m1$residuals)))
                          const <- switch(as.character(aug),
                                          "FALSE" = switch(estimand,
                                                           "ATT" = 0,
                                                           "ATC" = 0,
                                                           "ATE" = c(0,
                                                                     0)),
                                          "TRUE" = switch(estimand,
                                                          "ATT" = mean(design %*% m0$coefficients),
                                                          "ATC" = mean(design %*% m1$coefficients),
                                                          "ATE" = c(mean(design %*% m0$coefficients),
                                                                    mean(design %*% m1$coefficients))))
                          
                          
                          mu <- switch(estimand,
                                           "ATE" = c(mean(private$mu0), mean(private$mu1)),
                                           "ATT" = mean(private$mu0[private$z==1]),
                                           "ATC" = mean(private$mu1[private$z==0]))
                          
                          # f <- switch(estimand,
                          #             "ATE" = f,
                          #             "ATT" = f[private$z == 0],
                          #             "ATC" = f[private$z == 1])
                          if(estimand != "ATE" ) {
                            Q <- Matrix::Matrix(tcrossprod(f), sparse = TRUE)
                            if(solver == "mosek") Q <- as(as(Q, "generalMatrix"), "CsparseMatrix")
                            L <- -2 * (mu - const) * t(f)
                            A <- Matrix::sparseMatrix(i = rep(1, length(f)),
                                                      j = 1:length(f),
                                                      x = 1) # sum to 1 constraint
                            problem <- list(obj = list(Q = Q, L = L),
                                            LC = list(A = A,
                                                      lc = rep(1, nrow(A)),
                                                      uc = rep(1, nrow(A))),
                                            nvar = length(L),
                                            bounds = list(lb = rep(0, length(L)),
                                                          ub = rep(Inf, length(L))))
                            
                            weights <- switch(solver,
                                              "mosek" =  mosek_solver(problem),
                                              "gurobi" = gurobi_solver(problem),
                                              "cplex"  = cplex_solver(problem),
                                              "quadprog" = quadprog_solver(problem)
                            )
                            weights <- renormalize(weights$sol)
                          } else {
                            f1 <- f[[2]]
                            f0 <- f[[1]]
                            
                            Q0 <- Matrix::Matrix(tcrossprod(f0), sparse = TRUE)
                            if(solver == "mosek") Q0 <- as(as(Q0, "generalMatrix"), "CsparseMatrix")
                            
                            L0 <- -2 * (mu[1] - const[1]) * t(f0)
                            A0 <- Matrix::sparseMatrix(i = rep(1, length(f0)),
                                                       j = 1:length(f0),
                                                       x = 1)
                            problem0 <- list(obj = list(Q = Q0, L = L0),
                                             LC = list(A = A0,
                                                       lc = rep(1, nrow(A0)),
                                                       uc = rep(1, nrow(A0))),
                                             nvar = length(L0),
                                             bounds = list(lb = rep(0,length(L0)),
                                                           ub = rep(Inf, length(L0))))
                            
                            Q1 <- Matrix::Matrix(tcrossprod(f1), sparse = TRUE)
                            if(solver == "mosek") Q1 <- as(as(Q1,"generalMatrix"), "CsparseMatrix")
                            
                            L1 <- -2 * (mu[2] - const[2]) * t(f1)
                            A1 <- Matrix::sparseMatrix(i = rep(1, length(f1)),
                                                       j = 1:length(f1),
                                                       x = 1)
                            problem1 <- list(obj = list(Q = Q1, L = L1),
                                             LC = list(A = A1,
                                                       lc = rep(1, nrow(A1)),
                                                       uc = rep(1, nrow(A1))),
                                             nvar = length(L1),
                                             bounds = list(lb = rep(0,length(L1)),
                                                           ub = rep(Inf, length(L1))))
                            
                            weights <- list(w0 = switch(solver,
                                                        "mosek" =  mosek_solver(problem0),
                                                        "gurobi" = gurobi_solver(problem0),
                                                        "cplex"  = cplex_solver(problem0),
                                                        "quadprog" = quadprog_solver(problem0)
                            ),
                            w1 = switch(solver,
                                        "mosek" =  mosek_solver(problem1),
                                        "gurobi" = gurobi_solver(problem1),
                                        "cplex"  = cplex_solver(problem1),
                                        "quadprog" = quadprog_solver(problem1)
                            ) )
                            weights <- lapply(weights, function(w) renormalize(w$sol))
                          }
                          return(weights)
                          
                        },
                      #' @description
                      #' Gets the distance of the weights from the optimal weights
                      #' @param weight The estimated weights 
                      #' @param estimand One of "ATT","ATC","ATE"
                      #' @param augment Should we use an augmented estimator? TRUE or FALSE.
                      #' @param solver One of "mosek", "gurobi", or "quadprog"
                      opt_weight_dist = function(weight, estimand = "ATE", augment = FALSE, solver = "mosek") {
                        if(estimand == "cATE" | estimand == "feasible") estimand <- "ATE"
                        estimand <- match.arg(estimand, choices = c("ATT","ATC","ATE"))
                        
                        w_star <- self$opt_weight(estimand = estimand, augment = augment, solver = solver)
                        # norm_weight <- list(w0 = renormalize(weight$w0),
                        #                     w1 = renormalize(weight$w1))
                        if(estimand == "ATE") {
                          dist <- mean(abs(unlist(w_star) - c(weight$w0, weight$w1)))
                        } else if (estimand == "ATT") {
                          dist <- mean(abs((w_star - c(weight$w0))))
                        } else if (estimand == "ATC") {
                          dist <- mean(abs((w_star - c(weight$w1))))
                        }
                        return(dist)
                      }
                        
                      ),
         private = list(n = "numeric",
                       p = "numeric",
                       x = "matrix",
                       y = "vector",
                       z = "vector",
                       # index = "vector",
                       # tx_index = "vector",
                       # ctrl_index = "vector",
                       param = "list",
                       n1 = "numeric",
                       n0 = "numeric",
                       x0 = "vector",
                       x1 = "vector",
                       mu1 = "vector",
                       mu0 = "vector",
                       check_data = function() {
                         complete <- all(is.matrix(private$x) & is.vector(private$z) )
                         if (complete) {
                           private$n1 <- sum(private$z == 1)
                           private$n0 <- sum(private$z == 0)
                           private$x1 <- private$x[private$z == 1,,drop = FALSE]
                           private$x0 <- private$x[private$z == 0,,drop = FALSE]
                         }
                       }#,
                       )
)


#### Hainmueller: 0 tx effect and mix of covariates (normal, binary, etc.) ####

#' Hainmueller data example
#' 
#' @details Generates the data as described in Hainmueller (2012).
#' 
#' @return An [R6][R6::R6Class] object of class [DataSim][DataSim]
#' @export
#' @rdname DataSim-Hainmueller
Hainmueller <- R6::R6Class("Hainmueller", 
                           inherit = DataSim,
                           public = list(
                             #' @description
                             #' Generates the data
                             gen_data = function() {
                               self$gen_x()
                               self$gen_z()
                               self$gen_y()
                               # private$\check_data()
                               invisible(self)
                             },
                             #' @description
                             #' Generates the covaraiate data
                             gen_x = function() {
                               stopifnot(length(private$n) > 0 )
                               x13 <- matrix(private$param$param_x$x_13$mean, nrow = private$n,
                                             ncol = 3, byrow = TRUE) + 
                                 matrix(rnorm(private$n * 3), 
                                        nrow = private$n, 
                                        ncol = 3) %*% chol(private$param$param_x$x_13$covar)
                               x4 <- runif(private$n, private$param$param_x$x4$lower, private$param$param_x$x4$upper)
                               x5 <- rchisq(private$n, df = private$param$param_x$x5$df)
                               x6 <- rbinom(private$n, size = 1, prob = private$param$param_x$x6$p)
                               
                               private$x <- cbind(x13, x4, x5, x6)
                               colnames(private$x) <- paste0("X",1:6)
                               private$check_data()
                               invisible(self)
                             },
                             #' @description
                             #' Generates the outcome data
                             gen_y = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               mean_y <- private$mu0 <- private$mu1 <- 
                                 if(private$design =="A" || private$design == 1) {
                                   private$x %*% private$param$beta_y
                                 } else if (private$design =="B" || private$design == 2) {
                                   (private$x[,c(1,2,5)] %*% private$param$beta_y[1:3])^2
                                 } else if (private$design == "C" || private$design == 3) {
                                   cos(private$x[,1] * private$x[,3])^2 * private$param$beta_y[1] +
                                     private$x[,1]^2 * private$x[,5] * private$param$beta_y[2] +
                                     private$x[,3]^2 * private$param$beta_y[3]
                                 } else {
                                   stop("design must be one of 'A' or 'B' or 'C'")
                                 }
                               private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = private$param$sigma_y))
                               # private$check_data()
                               invisible(self)
                             },
                             #' @description
                             #' Generates the treatment indicator
                             gen_z = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               mean_z <- private$x %*% private$param$beta_z
                               latent_z <- if(private$overlap != "medium") {
                                 mean_z + rnorm(private$n, mean=0, sd = private$param$sigma_z)
                               } else {
                                 mean_z + (rchisq(private$n, df = 5) - 5) * private$param$sigma_z + 0.5
                               }
                               
                               private$pscore <- if(private$overlap != "medium") {
                                 pnorm(0, mean = mean_z, sd = private$param$sigma_z, lower.tail = FALSE)
                               } else {
                                 pchisq(q = ((-mean_z - 0.5) / private$param$sigma_z + 5), 
                                        df = 5, lower.tail = FALSE)
                               }
                               private$z <- c(ifelse(latent_z < 0, 0, 1))
                               private$check_data()
                               invisible(self)
                             },
                              #' @description
                              #' Generates the the Hainmueller `R6` class
                              #' @param n The number of observations
                              #' @param p The dimensions of the covariates. 
                              #' Fixed to 6.
                              #' @param param The data generating parameters
                              #' fed as a list.
                              #' @param design One of "A" or "B". See details.
                              #' @param overlap One of "high", "low", or "medium".
                              #' See details.
                              #' @param ... Extra arguments. Currently unused.
                              #' 
                              #' @details
                              #' ## Design
                              #' Design "A"
                              #' is the setting where the outcome is generated 
                              #' from a linear model, 
                              #' \eqn{Y(0) = Y(1) = X_1 + X_2 + X_3 - X_4 + X_5 + X_6 + \eta} 
                              #' and design "B" is where the outcome is 
                              #' generated from the non-linear model 
                              #' \eqn{Y(0) = Y(1) = (X_1 + X_2 +X_5 )^2 + \eta}.
                              #' 
                              #' ## Overlap
                              #' The treatment indicator is generated from
                              #' \eqn{Z = 1(X_1 + 2 X_2 - 2 X_3 - X_4 - 0.5 X_5 + X_6 + \nu > 0)}, where \eqn{\nu} 
                              #' depends on the overlap selected. If overlap is "high",
                              #' then \eqn{\nu \sim N(0, 100).} If overlap is 
                              #' "low", then \eqn{\nu \sim N(0, 30).} Finally,
                              #' if overlap is "medium", then \eqn{\nu} is drawn
                              #' from a \eqn{\chi^2} with 5 degrees of freedom
                              #' that is scaled and centered to have mean 0.5 and 
                              #' variance 67.6.
                              #'
                              #' @return An object of class [DataSim][DataSim].
                              #' @export
                              #'
                              #' @examples
                              #' data <- Hainmueller$new(n = 100, p = 6, design = "A", overlap = "low")
                              #' data$gen_data()
                              #' print(data$get_x()[1:2,])
                             initialize = function(n = 100, p = 6, param = list(), design = "A", overlap = "low", ...) {
                               
                               if(p != 6) warning("'p' set to 6 automatically")
                               private$p <- 6 # p is always 6 for this guy
                               
                               if(missing(n) | is.null(n)) {
                                 private$n <- 100
                               } else {
                                 private$n <- n
                               }
                               if(missing(design ) | is.null(design) ) {
                                 private$design <- "A"
                               } else {
                                 private$design <- match.arg(design, c("A","B","C"))
                               }
                               if( missing(overlap) | is.null(overlap) ) {
                                 private$overlap <- "low"
                               } else {
                                 private$overlap <- match.arg(overlap, c("low","medium","high"))
                               }
                               private$
                                 set_param(beta_z = param$beta_z, beta_y = param$beta_y,
                                           sigma_z = param$sigma_z, sigma_y = param$sigma_y,
                                           param_x = param$param_x)
                               
                             },
                             #' @description
                             #' Returns the chosen design parameters
                             get_design = function() {
                               return(c(design = private$design, overlap = private$overlap))
                             },
                             #' @description
                             #' Returns the true propensity score
                             get_pscore = function() {
                               return(private$pscore)
                             }
                           ),
                           private = list(design = "character",
                                          overlap = "character",
                                          pscore = "numeric",
                                          set_param = function(beta_z, beta_y, sigma_z, sigma_y, param_x) {
                                            miss.null <- function(xx) {
                                              return(missing(xx) | is.null(xx))
                                            }
                                            if(is.null(private$design) & (miss.null(beta_y) ) ) {
                                              private$design <- "A"
                                            }
                                            if(is.null(private$overlap) & (miss.null(sigma_z) )) {
                                              private$overlap <- "low"
                                            }
                                            default_param <- list(
                                              beta_z = c(1,2,-2,-1,-0.5,1),
                                              beta_y = list(A = c(1,1,1,-1,1,1),
                                                            B = c(1,1,1),
                                                            C = c(1, -1, 1)),
                                              sigma_z= list(low = sqrt(30),
                                                            medium = sqrt(67.6/10),
                                                            high = sqrt(100)),
                                              sigma_y = 1,
                                              param_x = list(x_13 = list(mean = rep(0, 3),
                                                                         covar = matrix(c(2,1,-1,1,1,-0.5, -1, -0.5, 1), nrow = 3,ncol = 3)),
                                                             x4 = list(lower = -3, upper = 3),
                                                             x5 = list(df = 1),
                                                             x6 = list(p = 0.5))
                                            )
                                            temp_param <- list()
                                            if(miss.null(beta_z)) {
                                              temp_param$beta_z <- default_param$beta_z
                                            } else {
                                              stopifnot(is.vector(param$beta_z))
                                              temp_param$beta_z <- param$beta_z
                                            }
                                            if(miss.null(sigma_z)) {
                                              temp_param$sigma_z <- default_param$sigma_z[[private$overlap]]
                                            } else {
                                              stopifnot(is.numeric(sigma_z))
                                              temp_param$sigma_z <- param$sigma_z
                                            }
                                            if(miss.null(beta_y)) {
                                              temp_param$beta_y <- default_param$beta_y[[private$design]]
                                            } else {
                                              stopifnot(is.vector(beta_y))
                                              temp_param$beta_y <- param$beta_y
                                            }
                                            if(miss.null(sigma_y)) {
                                              temp_param$sigma_y <- default_param$sigma_y
                                            } else {
                                              temp_param$sigma_y <-param$sigma_y
                                            }
                                            if(miss.null(param_x)) {
                                              temp_param$param_x <- default_param$param_x
                                            } else {
                                              names_param_x <- c('x_13','x4','x5','x6')
                                              if(miss.null(param_x$x_13)) {
                                                temp_param$param_x$x_13 <- default_param$param_x$x_13
                                              } else {
                                                stopifnot(length(param_x$x_13$mean)==3)
                                                stopifnot(all(dim(param_x$x_13$covar) %in% 3 ))
                                                temp_param$param_x$x_13 <-param$param_x$x_13
                                              }
                                              if(miss.null(param_x$x_4)) {
                                                temp_param$param_x$x_4 <- default_param$param_x$x_4
                                              } else {
                                                stopifnot(is.numeric(param_x$x_4$lower))
                                                stopifnot(is.numeric(param_x$x_4$upper))
                                                temp_param$param_x$x_4 <-param$param_x$x_4
                                              }
                                              if(miss.null(param_x$x_5)) {
                                                temp_param$param_x$x_5 <- default_param$param_x$x_5
                                              } else {
                                                stopifnot(is.numeric(param_x$x_5$df))
                                                stopifnot(is.numeric(param_x$x_5$df>0))
                                                temp_param$param_x$x_5 <-param$param_x$x_5
                                              }
                                              if(miss.null(param_x$x_6)) {
                                                temp_param$param_x$x_6 <- default_param$param_x$x_6
                                              } else {
                                                stopifnot(is.numeric(param_x$x_6$p))
                                                stopifnot(is.numeric(param_x$x_6$p>0))
                                                stopifnot(is.numeric(param_x$x_6$p<1))
                                                temp_param$param_x$x_6 <-param$param_x$x_6
                                              }
                                            }
                                            private$param <- temp_param
                                          }
                           )
)

#### Kallus2019: continuous tx ####

Kallus2019 <- R6::R6Class("Kallus2019", 
                           inherit = DataSim,
                           public = list(
                             gen_data = function() {
                               self$gen_x()
                               self$gen_z()
                               self$gen_y()
                               invisible(self)
                             },
                             gen_x = function() {
                               stopifnot(length(private$n) >0 )
                               private$x <- matrix(rnorm(private$n * private$p, 
                                                         mean = private$param$param_x$mean, 
                                                         sd =  private$param$param_x$sd), 
                                                   nrow = private$n, ncol = private$p)
                               colnames(private$x) <- paste0("X",1:5)
                               private$check_data()
                               invisible(self)
                             },
                             gen_y = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               if(all(dim(private$z) == 0)) gen_z()
                               mean_y <- 0.75 * private$z + 
                                 0.05 * private$z^2 + 0.001 * private$z^3 +
                                 1.5 * rowSums(private$x) + 
                                 1.125 * rowSums(private$x)
                               
                               private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = private$param$sigma_y))
                               private$check_data()
                               invisible(self)
                             },
                             gen_z = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               design_z <- cbind(1, rowSums(private$x)^private$d)
                               mean_z <- design_z %*% private$param$beta_z
                               private$z <- c(mean_z + rnorm(private$n, 0, private$param$sigma_z))
                               private$check_data()
                               invisible(self)
                             },
                             initialize = function(n = 1024, p = 5, param = list(), design = "A", ...) {
                               
                               if(p != 5) warning("'p' set to 5 automatically")
                               private$p <- 5 # p is always 6 for this guy
                               
                               if(missing(n) | is.null(n)) {
                                 private$n <- 1024
                               } else {
                                 private$n <- n
                               }
                               if(missing(design ) | is.null(design) ) {
                                 private$design <- "A"
                               } else {
                                 private$design <- match.arg(design, c("A","B","C"))
                               }
                               private$
                                 set_param(
                                           beta_z = param$beta_z, #beta_y = param$beta_y,
                                           sigma_z = param$sigma_z, 
                                           sigma_y = param$sigma_y,
                                           param_x = param$param_x)
                               
                             },
                             get_design = function() {
                               return(switch(private$design,
                                      A = "A: linear",
                                      B = "B: quadratic",
                                      C = "C: cubic"))
                             }
                           ),
                           private = list(d = "numeric",
                                          design = "character",
                                          # overlap = "character",
                                          set_param = function(beta_z, sigma_z, sigma_y, param_x) {
                                            miss.null <- function(xx) {
                                              return(missing(xx) | is.null(xx))
                                            }
                                            if (is.null(private$design) ) {
                                              private$design <- "A"
                                            }
                                            private$d <- switch(private$design,
                                                                A = 1, 
                                                                B = 2,
                                                                C = 3)
                                            default_param <- list(
                                              # beta_z = c(1,2,-2,-1,-0.5,1),
                                              beta_z = switch(private$design,
                                                              A = c(0,1),
                                                              B = c(-3, 0.25),
                                                              C = c(-2.5, 0.05)),
                                              sigma_z= sqrt(5),
                                              sigma_y = 0,
                                              param_x = list(mean = 0,
                                                              sd = sqrt(5))
                                            )
                                            temp_param <- list()
                                            if(miss.null(beta_z)) {
                                              temp_param$beta_z <- default_param$beta_z
                                            } else {
                                              stopifnot(is.vector(param$beta_z))
                                              temp_param$beta_z <- param$beta_z
                                            }
                                            if(miss.null(sigma_z)) {
                                              temp_param$sigma_z <- default_param$sigma_z[[private$overlap]]
                                            } else {
                                              stopifnot(is.numeric(sigma_z))
                                              temp_param$sigma_z <- param$sigma_z
                                            }
                                            # if(miss.null(beta_y)) {
                                            #   temp_param$beta_y <- default_param$beta_y[[private$design]]
                                            # } else {
                                            #   stopifnot(is.vector(beta_y))
                                            #   temp_param$beta_y <- param$beta_y
                                            # }
                                            if(miss.null(sigma_y)) {
                                              temp_param$sigma_y <- default_param$sigma_y
                                            } else {
                                              temp_param$sigma_y <-param$sigma_y
                                            }
                                            if(miss.null(param_x)) {
                                              temp_param$param_x <- default_param$param_x
                                            } else {
                                              if(is.null(param$param_x$mean) | is.null(param$param_x$sd)) stop("Must specify parameters of x as list(mean = , sd = )")
                                              temp_param$param_x$mean <- param$param_x$mean
                                              temp_param$param_x$sd <- param$param_x$sd
                                            }
                                            private$param <- temp_param
                                          }
                           )
)


#### Kallus2018: binary tx ####
  Kallus2018 <- R6::R6Class("Kallus2018", 
                            inherit = DataSim,
                            public = list(
                              gen_data = function() {
                                self$gen_x()
                                self$gen_z()
                                self$gen_y()
                                invisible(self)
                              },
                              gen_x = function() {
                                stopifnot(length(private$n) >0 )
                                private$x <- matrix(rnorm(private$n * private$p, 
                                                          mean = private$param$param_x$mean, 
                                                          sd =  private$param$param_x$sd), 
                                                    nrow = private$n, ncol = private$p)
                                colnames(private$x) <- paste0("X",1:4)
                                private$check_data()
                                invisible(self)
                              },
                              gen_y = function() {
                                if(all(dim(private$x) == 0)) gen_x()
                                if(is.character(private$z)) gen_z()
                                
                                private$mu0 <- if(private$design == "A"){ 
                                  -mean(private$pi) + 0 + rowSums(private$x[,1:2])
                                } else if(private$design == "B"){
                                  -mean(private$pi) + 0 + rowSums(private$x[,1:2]) + 
                                  rowSums(private$x[,1:2]^2) + apply(private$x[,1:2],1,prod)
                                }
                                
                                private$mu1 <- if(private$design == "A"){ 
                                  -mean(private$pi) + 1 + rowSums(private$x[,1:2])
                                } else if(private$design == "B"){
                                  -mean(private$pi) + 1 + rowSums(private$x[,1:2]) + 
                                    rowSums(private$x[,1:2]^2) + apply(private$x[,1:2],1,prod)
                                }
                                mean_y <- private$mu1 * private$z + private$mu0 *(1-private$z)
                                private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = private$param$sigma_y))
                                private$check_data()
                                invisible(self)
                              },
                              gen_z = function() {
                                if(all(dim(private$x) == 0)) gen_x()
                                if(private$design == "A") {
                                  design.z <- rowSums(private$x[,1:2])
                                } else if (private$design == "B") {
                                  design.z <- rowSums(private$x[,1:2]) + 
                                    rowSums(private$x[,1:2]^2) + apply(private$x[,1:2],1,prod)
                                }
                                if(private$overlap == "low") {
                                  private$pi <- plogis(3 * design.z)
                                } else if (private$overlap == "high") {
                                  private$pi <- plogis(0.1 * design.z)
                                }
                                private$z <- rbinom(private$n, size = 1, prob = private$pi)
                                invisible(self)
                              },
                              initialize = function(n = 1024, param = list(), design = "A", overlap = "high", ...) {
                                
                                # if(p != 4) warning("'p' set to 4 automatically")
                                private$p <- 4 # p is always 4 for this guy
                                
                                if(missing(n) | is.null(n)) {
                                  private$n <- 1024
                                } else {
                                  private$n <- n
                                }
                                if(missing(design ) | is.null(design) ) {
                                  private$design <- "A"
                                } else {
                                  private$design <- match.arg(design, c("A","B"))
                                }
                                if(missing(overlap ) | is.null(overlap) ) {
                                  private$overlap <- "high"
                                } else {
                                  private$overlap <- match.arg(overlap, c("high","low"))
                                }
                                
                              },
                              get_design = function() {
                                return(paste0(switch(private$design,
                                              A = "A: linear",
                                              B = "B: quadratic"),
                                              ", overlap: ",
                                              switch(private$overlap,
                                                     high = "high",
                                                     low = "low"))
                                )
                              }
                            ),
                            private = list(d = "numeric",
                                           design = "character",
                                           overlap = "character",
                                           pi = "numeric"
                            )
  )

# Sonabed2020: indicator function mean model. added
# skew normal covariates to reduce overlap without specifying
# propensity score function
{
  Sonabend2020 <- R6::R6Class("Sonabend2020", 
                           inherit = DataSim,
                           public = list(
                             gen_data = function() {
                               self$gen_x()
                               self$gen_z()
                               self$gen_y()
                               invisible(self)
                             },
                             gen_x = function() {
                               stopifnot(length(private$n) >0 )
                               rskew <- function(n, alpha) {
                                 x1 <- rnorm(n)
                                 x2 <- rnorm(n)
                                 a <- (alpha*abs(x1) + x2)/sqrt(1 + alpha^2)
                                 b <- (1+alpha)/sqrt(2*(1+alpha^2)) * pmax(x1,x2) +
                                   (1-alpha)/sqrt(2*(1+alpha^2)) * pmin(x1,x2)
                                 return(cbind(a,b))
                               }
                               
                               if(private$overlap == "high") {
                                 alpha <- 0.1
                               } else if (private$overlap == "low") {
                                 alpha <- 1
                               }
                               n0 <- floor(private$n/2)
                               n1 <- ceiling(private$n/2)
                               pp <- private$p - 2
                               pp_skew <- floor(pp / 2)
                               alpha_extra <- runif(pp_skew, -0.5,.5)
                               pp_nonskew <- pp - pp_skew * 2
                               
                               v <- rskew(n0, alpha)
                               w <- rskew(n1, -alpha)
                               extra_skew0 <- matrix(replicate(pp_skew, rskew(n0, alpha_extra)), 
                                                     nrow = n0,
                                                     ncol = pp_skew * 2)
                               extra_skew1 <- matrix(replicate(pp_skew, rskew(n1, alpha_extra)),
                                                        nrow = n0,
                                                        ncol = pp_skew * 2)
                               
                               x0 <- cbind( 
                                           # rnorm(n0, mean = v[,1]^2),
                                           v,
                                           extra_skew0,
                                           matrix(rnorm(n0 * pp_nonskew), n0, pp_nonskew))
                               x1 <- cbind(
                                           # rnorm(n1, mean = (w[,1]^2)),
                                           w, 
                                           extra_skew1,
                                           matrix(rnorm(n0 * pp_nonskew), n1, pp_nonskew))
                               
                               private$x <- rbind(x0,x1)
                               colnames(private$x) <- paste0("X",1:private$p)
                               invisible(self)
                             },
                             gen_y = function() {
                               if(all(dim(private$x) == 0)) self$gen_x()
                               if(is.character(private$z)) self$gen_z()
                               private$mu1 <- if(private$design == "A"){
                                 1 + (private$x[,2] > 0)
                               } else if (private$design == "B") {
                                 ytemp <- (round(cos(private$x[,2])))
                                 # ytemp <- 0
                                 10*sign(1/(.01*private$x[,1]+.05*private$x[,2])+ ytemp*100) + private$x[,3]^2/10 
                                 # 5*(sign(1/(.01*private$x[,1]+.05*private$x[,2])) + private$x[,1]^2/10 + ytemp)
                               } else if (private$design == "C") {
                                 private$x[,1:3] %*% c(-1,1,-2) + 0.5 *private$x[,3]^2
                               } else {
                                 stop("Design must be one of A, B, or C")
                               }
                               private$mu0 <- if(private$design == "A"){
                                 -1 - (1.5*abs(private$x[,1]) > 1)
                               }else if (private$design == "B") {
                                 ytemp <- -(1 + (1.5*abs(private$x[,2]) > 1))
                                 # ytemp <- 0
                                 # 5*(sign(-1/(.01*private$x[,1]+.05*private$x[,2])) - private$x[,1]^2/10 + ytemp)
                                 
                                 10*sign(-.05*private$x[,1]+.01*private$x[,2]) + ytemp*5 - private$x[,3]^3/10
                               } else if (private$design == "C") {
                                 private$x[,1:3] %*% c(-1,-1,-2) - 2
                               }
                               mean_y <- private$mu1 * private$z + 
                                 (1 - private$z) * private$mu0
                               private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = 1))
                               private$check_data()
                               invisible(self)
                             },
                             gen_z = function() {
                               if(all(dim(private$x) == 0)) self$gen_x()
                               private$z <- c(rep(0, floor(private$n/2)),
                                                rep(1, ceiling(private$n/2)))
                               private$check_data()
                               invisible(self)
                             },
                             initialize = function(n = 1024, p = 4, design = "A", overlap = "high", ...) {
                               
                               if(p < 3) {
                                 warning("'p' must be at least 4")
                                 private$p <- 4
                               } else {
                                 private$p <- p
                               }
                               
                               if(missing(n) | is.null(n)) {
                                 private$n <- 1024
                               } else {
                                 private$n <- n
                               }
                               if(missing(design ) | is.null(design) ) {
                                 private$design <- "A"
                               } else {
                                 private$design <- match.arg(design, c("A","B"))
                               }
                               if(missing(overlap ) | is.null(overlap) ) {
                                 private$overlap <- "high"
                               } else {
                                 private$overlap <- match.arg(overlap, c("high","low"))
                               }
                               
                             },
                             get_design = function() {
                               return(paste0(switch(private$design,
                                                    A = "A: indicator",
                                                    B = "B: sign function",
                                                    C = "C: linear"),
                                             ", overlap: ",
                                             switch(private$overlap,
                                                    high = "high",
                                                    low = "low"))
                               )
                             }
                           ),
                           private = list(design = "character",
                                          overlap = "character"
                           )
)
}


#### Kang and Schafer data ####

  KangSchafer <- R6::R6Class("KangSchafer", 
                             inherit = DataSim,
                             public = list(
                               gen_data = function() {
                                 self$gen_x()
                                 self$gen_z()
                                 self$gen_y()
                                 invisible(self)
                               },
                               gen_x = function() {
                                 stopifnot(length(private$n) > 0 )
                                 
                                 private$u <- matrix(stats::rnorm(private$n * private$p),
                                                     nrow = private$n, ncol = private$p)
                                 if (private$design %in% c("A", 1) ) {
                                   private$x <- private$u
                                 } else {
                                   private$x <- cbind(
                                     exp(private$u[,1]/2),
                                     private$u[,2]/(1 + exp(private$u[,1])) + 10,
                                     (private$u[,1] * private$u[,3] / 25 + 0.6)^3,
                                     (private$u[,2] + private$u[,4] + 20)^2
                                   )
                                   if (private$p > 4)
                                     private$x <- cbind(private$x, private$u[,5:private$p,drop = FALSE])
                                 }
                                 colnames(private$u) <- paste0("U",1:private$p)
                                 colnames(private$x) <- paste0("X",1:private$p)
                                 private$check_data()
                                 invisible(self)
                               },
                               gen_y = function() {
                                 if (all(dim(private$u) == 0)) gen_x()
                                 if (length(private$x) == 0) gen_z()
                                 
                                 
                                 private$mu0 <- 
                                   210 + 27.4*private$u[,1] + 13.7*private$u[,2] + 
                                   13.7*private$u[,3] + 13.7*private$u[,4]
                                 
                                 private$mu1 <- private$mu0 + private$param$tau 
                                 mean_y <- private$mu0 * (private$z == 0) + 
                                   private$mu1 * (private$z == 1)
                                 private$y <- c(mean_y + rnorm(private$n, 
                                                               mean = 0, 
                                                               sd = private$param$sigma_y))
                                 # private$check_data()
                                 invisible(self)
                               },
                               gen_z = function() {
                                 if (all(dim(private$x) == 0)) gen_x()
                                 
                                 #set coefficient
                                 beta_z <- if ( private$overlap == "high") {
                                   c(-1, 0.5, -0.25, -0.1)
                                 } else if ( private$overlap == "low" ) {
                                   c(-1, -0.5, -0.25, -0.1)
                                 }
                                 
                                 #probability
                                 prob_z <- plogis(private$u[,1:4] %*% beta_z)
                                 
                                 #assign z
                                 private$z <- rbinom(private$n, 
                                                     1, 
                                                     prob = prob_z)
                                 private$check_data()
                                 invisible(self)
                               },
                               get_u = function() {
                                 return(private$u)
                               },
                               initialize = function(n = 100, p = 4, design = "A", overlap = "low", 
                                                     param = list(), ...) {
                                 
                                 if (isTRUE(is.null(p)) | isTRUE(p < 4)) {
                                   warning("'p' must be >= 4. Set to 4.")
                                  private$p <- 4
                                 } else {
                                   private$p <- p
                                 }
                                 
                                 if (missing(n) | is.null(n)) {
                                   private$n <- 100
                                 } else {
                                   private$n <- n
                                 }
                                 if (missing(design ) | is.null(design) ) {
                                   private$design <- "A"
                                 } else {
                                   private$design <- match.arg(design, c("A","B"))
                                 } 
                                 if ( missing(overlap) | is.null(overlap) ) {
                                   private$overlap <- "low"
                                 } else {
                                   private$overlap <- match.arg(overlap, c("low","high"))
                                 }
                                 
                                 if ( missing(param) | is.null(param)) {
                                   private$param <- list(tau = 0,
                                                         sigma_y = 1)
                                 } else {
                                   private$param <- param[c("tau", "sigma_y")]
                                   if (is.null(param$sigma_y) ) param$sigma_y <- 1
                                   if (is.null(param$tau) ) param$tau <- 0.0
                                   
                                   if (!is.numeric(private$param$tau)) stop("treatment effect `tau` must be a real number")
                                   if (!is.numeric(private$param$sigma_y) | private$param$sigma_y <= 0) stop("variance of outcome `sigma_y` must be a positive real number")
                                 }
                               },
                               get_design = function() {
                                 return(c(design = private$design, overlap = private$overlap))
                               }
                             ),
                             private = list(design = "character",
                                            overlap = "character",
                                            u = "matrix"
                                            
                             )
  )


#### LaLonde data ####
  #' LaLonde data example
  #' 
  #' @details Returns the LaLonde data as used by Dehjia and Wahba. Note the data
  #' is fixed and `gen_data()` will just initialize the fixed data.
  #' 
  #' @return An [R6][R6::R6Class] object of class [DataSim][DataSim]
  #' @export
  #' @rdname DataSim-LaLonde
  LaLonde <- R6::R6Class("LaLonde", 
                             inherit = DataSim,
                             public = list(
                               #' @description
                               #' Sets up the data
                               gen_data = function() {
                                 self$gen_x()
                                 self$gen_z()
                                 self$gen_y()
                                 invisible(self)
                               },
                               #' @description
                               #' Returns the experimental treatment effect, $1794
                               get_tau = function() {
                                 return(1794)
                               },
                               #' @description
                               #' Sets up the covariate data
                               gen_x = function() {
                                 form <- as.formula("~. + I(as.numeric(re74 == 0)) + I(as.numeric(re75 == 0)) + + 0")
                                 if (private$design == "Full") {
                                   cn <- colnames(lalonde_full)
                                   
                                   private$x <- model.matrix(object = form, 
                                                             data = lalonde_full[,-match(c("data_id","treat","re78"),cn)])
                                 } else if (private$design == "NSW") {
                                   cn <- colnames(lalonde_nsw)
                                   private$x <- model.matrix(object = form, 
                                                             data = lalonde_nsw[,-match(c("data_id","treat","re78"),cn)])
                                 } else {
                                   stop("Design not found")
                                 }
                                 colnames(private$x) <- c(colnames(private$x)[1:8], "u74", "u75")
                                 private$p <- ncol(private$x)
                                 private$n <- nrow(private$x)
                                 private$check_data()
                                 invisible(self)
                               },
                               #' @description
                               #' Sets up the outcome data
                               gen_y = function() {
                                 if (private$design == "Full") {
                                   cn <- colnames(lalonde_full)
                                   
                                   private$y <- c( unlist(lalonde_full[,match("re78",cn)]) )
                                 } else if (private$design == "NSW") {
                                   cn <- colnames(lalonde_nsw)
                                   private$y <- c( unlist(lalonde_nsw[,match("re78",cn)]) )
                                 } else {
                                   stop("Design not found")
                                 }
                                 private$mu1 <- private$mu0  <- rep(NA, private$n)
                                 private$mu1[private$z == 1] <- private$y[private$z == 1]
                                 private$mu0[private$z == 1] <- private$mu1[private$z == 1] - 1794
                                 
                                 private$mu0[private$z == 0] <- private$y[private$z == 0]
                                 private$mu1[private$z == 0] <- private$mu0[private$z == 0] + 1794
                                 invisible(self)
                               },
                               #' @description
                               #' Sets up the treatment indicator
                               gen_z = function() {
                                 if (private$design == "Full") {
                                   cn <- colnames(lalonde_full)
                                   
                                   private$z <- as.integer(unlist(lalonde_full[,match("treat",cn)]))
                                 } else if (private$design == "NSW") {
                                   cn <- colnames(lalonde_nsw)
                                   private$z <- as.integer(unlist(lalonde_nsw[,match("treat",cn)]))
                                 } else {
                                   stop("Design not found")
                                 }
                                 private$check_data()
                                 invisible(self)
                               },
                               #' @description
                               #' Initializes the LaLonde object.
                               #' 
                               #' @param n Not used. 
                               #' Maintained for symmetry with other 
                               #' DataSim objects.
                               #' @param p Not used.
                               #' Maintained for symmetry with other 
                               #' DataSim objects.
                               #' @param param Not used.
                               #' Maintained for symmetry with other 
                               #' DataSim objects.
                               #' @param design One of "NSW" or "Full". 
                               #' "NSW" uses the original experimental data
                               #' from the job training program while option "Full"
                               #' uses the treated individuals from
                               #' LaLonde's study and compares them to 
                               #' individuals from the 
                               #' Current Population Survey as controls.
                               #' @param ... Not used.
                               #' 
                               #' @examples 
                               #' nsw <- LaLonde$new(design = "NSW")
                               #' nsw$gen_data()
                               #' nsw$get_n()
                               #' 
                               #' obs.study <-  LaLonde$new(design = "Full")
                               #' obs.study$gen_data()
                               #' obs.study$get_n()
                               initialize = function(n = NULL, p = NULL, param = list(), design = "NSW", ...) {
                                 
                                 if (missing(design ) | is.null(design) ) {
                                   private$design <- "NSW"
                                 } else {
                                   private$design <- match.arg(design, c("NSW","Full"))
                                 }
                                 
                               },
                               #' @description
                               #' Returns the chosen design parameters
                               get_design = function() {
                                 return(c(design = private$design))
                               }
                             ),
                             private = list(design = "character"
                             )
  )


#### Fan et al ####

  FanEtAl <- R6::R6Class("FanEtAl", 
                             inherit = DataSim,
                             public = list(
                               gen_data = function() {
                                 self$gen_x()
                                 self$gen_z()
                                 self$gen_y()
                                 # private$\check_data()
                                 invisible(self)
                               },
                               gen_x = function() {
                                 stopifnot(length(private$n) > 0 )
                                 
                                 private$x <- matrix(
                                   stats::rnorm(private$n * private$p),
                                   nrow = private$n,
                                   ncol = private$p
                                 ) %*% sqrt_mat(private$param$sigma_x)
                                 colnames(private$x) <- paste0("X",1:private$p)
                                 private$check_data()
                                 invisible(self)
                               },
                               gen_y = function() {
                                 if (all(dim(private$x) == 0)) gen_x()
                                 if (length(private$z) == 0) gen_z()
                                 
                                 private$mu0 <- rep(0.0, private$n)
                                 private$mu1 <- c(10 + private$x %*% private$param$beta_y)
                                 
                                 private$y <- private$mu1 * private$z + rnorm(private$n, sd = private$param$sigma_y)
                                 # private$check_data()
                                 invisible(self)
                               },
                               gen_z = function() {
                                 if (all(dim(private$x) == 0)) gen_x()
                                 latent_z <- private$x %*% private$param$beta_z
                                 private$U <- runif(1)
                                 
                                 private$z <- c(ifelse(plogis(latent_z) < private$U, 0, 1))
                                 private$check_data()
                                 invisible(self)
                               },
                               initialize = function(n = 100, p = 100, numActive = 4, param = list(), 
                                                     design = "B", ...) {
                                 
                                 if (isTRUE(missing(numActive)) | isTRUE(is.null(numActive)) | isTRUE(numActive < 1) ) {
                                   warning("'numActive' must be greater than 1. Set to 4")
                                   private$p1 <- 4 # p is always 6 for this guy
                                 } else {
                                   private$p1 <- numActive
                                 }
                                 
                                 if (isTRUE(missing(p)) | isTRUE(is.null(p)) | isTRUE(p < private$p1) ) {
                                   warning("'p' must be greater than numActive. Set to 100.")
                                  private$p <- if (100 > private$p1) {
                                    100L
                                  } else {
                                    as.integer(private$p1) + 100L
                                  }
                                 } else {
                                   private$p <- as.integer(p)
                                 }
                                 
                                 if (isTRUE(missing(n)) | isTRUE(is.null(n))) {
                                   private$n <- 500
                                 } else {
                                   private$n <- n
                                 }
                                 if (missing(design ) | is.null(design) ) {
                                   private$design <- "A"
                                 } else {
                                   private$design <- match.arg(design, c("A","B"))
                                 }
                                 private$set_param(beta_z = param$beta_z, beta_y = param$beta_y,
                                             sigma_y = param$sigma_y,
                                             sigma_x = param$sigma_x)
                                 
                               },
                               get_design = function() {
                                 return(c(Design = private$design, "Active Coefficient Number" = private$param$p1))
                               }
                             ),
                             private = list( design = "character",
                                             p1 = "integer",
                                             U = "numeric",
                               set_param = function(beta_z, beta_y, sigma_y, sigma_x) {
                                              miss.null <- function(xx) {
                                                return(missing(xx) | is.null(xx))
                                              }
                                              pm1 <- private$p - 1
                                              Sigma_x <- switch(private$design,
                                                                "A" = diag(1, private$p, private$p),
                                                                "B" = rbind(
                                                                  cbind(1,t(rep(0, pm1))),
                                                                  cbind(rep(0, pm1), 
                                                                        matrix(0.5^(abs(matrix(1:pm1, pm1, pm1, byrow = TRUE) -
                                                                                          matrix(1:pm1,pm1,pm1))), 
                                                                               pm1, pm1)
                                                                  )
                                                                )
                                              )
                                              default_param <- list(
                                                beta_y = switch(private$design,
                                                           "A" = c(rep(1, private$p1), 
                                                                    rep(0, private$p - private$p1)),
                                                           "B" = c(sqrt(0.5^2/c((1 - 0.5^2) * 
                                                                                          (t(1/(1:private$p)^2) %*%
                                                                                             Sigma_x %*% 
                                                                                             (1/(1:private$p)^2) )))) *
                                                             1/(1:private$p)^2
                                                ),
                                                sigma_y = 1,
                                                beta_z = switch(private$design,
                                                                "A" = c(rep(0.5, private$p1), 
                                                                        rep(0, private$p - private$p1)),
                                                                "B" = c(sqrt(pi^2/3 * 0.5^2/c((1 - 0.5^2) * 
                                                                                               (t(1/(1:private$p)^2) %*%
                                                                                                  Sigma_x %*% 
                                                                                                  (1/(1:private$p)^2) )))) *
                                                                  1/(1:private$p)^2
                                                ),
                                                sigma_x = Sigma_x
                                              )
                                              temp_param <- list()
                                              if (miss.null(beta_z)) {
                                                temp_param$beta_z  <- default_param$beta_z
                                              } else {
                                                stopifnot(is.vector(param$beta_z))
                                                temp_param$beta_z  <- param$beta_z
                                              }
                                              if (miss.null(beta_y)) {
                                                temp_param$beta_y  <- default_param$beta_y
                                              } else {
                                                stopifnot(is.vector(beta_y))
                                                temp_param$beta_y  <- param$beta_y
                                              }
                                              if (miss.null(sigma_y)) {
                                                temp_param$sigma_y <- default_param$sigma_y
                                              } else {
                                                temp_param$sigma_y <- param$sigma_y
                                              }
                                              if (miss.null(sigma_x)) {
                                                temp_param$sigma_x <- default_param$sigma_x
                                              } else {
                                                temp_param$sigma_x <- param$sigma_x
                                              }
                                              private$param <- temp_param
                                            }
                             )
  )


#### Li and Li ####

  LiLi <- R6::R6Class("LiLi", 
         inherit = DataSim,
         public = list(
           gen_data = function() {
             self$gen_x()
             self$gen_z()
             self$gen_y()
             # private$\check_data()
             invisible(self)
           },
           gen_x = function() {
             stopifnot(length(private$n) > 0 )
             
             x4 <- rbinom(private$n, size = 1, 
                          prob = 0.5)
             x3 <- rbinom(private$n, size = 1, 
                          prob = 0.6 * x4 + 0.4 (1 - x4))
             
             m1_2 <- cbind(rep(-x3 + x4 + 0.5 * x3 * x4,
                               private$n),
                           rep(x3 - x4 + x3 * x4,
                               private$n))
             half_vc1_2 <- sqrt_mat(x3 * matrix(c(1,0.5,0.5,1), 2, 2) + 
               (1 - x3) * matrix(c(2,0.25,0.25,2), 2, 2))
             x1_2 <- matrix(rnorm(private$n * 2), private$n, 2) %*% half_vc1_2 + m1_2
                            
             private$x <- cbind(
               x1_2, x3, x4
             )
             
             K <- private$p - 4
             if (K > 0) {
               x5_k2 <- matrix(rnorm(private$n * (4 + K/2), sd = sqrt(2)),
                               nrow = private$n,
                               ncol = 4 + K/2)
               x5pk2_4pK <- matrix(rbinomial(private$n * ((4 + K) - (4 + K/2)),
                                             size = 1, prob = 0.5),
                                   nrow = private$n,
                                   ncol = ((4 + K) - (4 + K/2)))
               
               private$x <- cbind(private$x,
                                  x5_k2,
                                  x5pk2_4pk)
             }
             
             colnames(private$x) <- paste0("X",1:private$p)
             private$check_data()
             invisible(self)
           },
           gen_y = function() {
             if (all(dim(private$x) == 0)) self$gen_x()
             if (length(private$z) == 0) self$gen_z()
             
             delta_0 <- private$param$tx_effect(rep(0, private$n))
             delta_1 <- private$param$tx_effect(rep(1, private$n))
             
             if (private$p == 4) {
               design <- switch(private$pscore_design,
                                "A" = private$x,
                                "B" = cbind(private$x, private$x[,1]^2, private$x[,1] * private$x[,2], private$x[,2]^2))
             } else {
               design <- private$x
               K <- private$p - 4
               if (K > 0) {
                 x5_k2 <- rowSums(private$x[, 5:(4 + K/2), drop = FALSE])
                 
                 x5pk2_4pK <- rowSums(private$x[, (5 + K/2):(4 + K), drop = FALSE])
                 
                 design <- cbind(design,
                                 x5_k2,
                                 x5pk2_4pk)
               }
               
             }
             
             tx_indep_mean <- cbind(1, design) %*% private$param$beta_y
             
             private$mu0 <- delta_0 + tx_indep_mean
             private$mu1 <- delta_1 + tx_indep_mean
             
             private$y <- private$mu1 * private$z + private$mu0 * ( 1 - private$z ) + rnorm(private$n, 
                                                                                            sd = private$param$sigma_y)
             # private$check_data()
             invisible(self)
           },
           gen_z = function() {
             if (all(dim(private$x) == 0)) self$gen_x()
             
             if (private$p == 4) {
               design <- switch(private$pscore_design,
                      "A" = private$x,
                      "B" = cbind(private$x, private$x[,1]^2, private$x[,1] * private$x[,2], private$x[,2]^2))
             } else {
               design <- private$x
               K <- private$p - 4
               if (K > 0) {
                 x5_k2 <- rowSums(private$x[, 5:(4 + K/2), drop = FALSE])
                                 
                 x5pk2_4pK <- rowSums(private$x[, (5 + K/2):(4 + K), drop = FALSE])
                 
                 design <- cbind(design,
                                    x5_k2,
                                    x5pk2_4pk)
               }
               
             }
             
             eta <- cbind(1, design) %*% private$param$beta_z
             private$pscore <- plogis(eta)
             private$z <- rbinom(private$n, 1, prob = private$pscore)
             private$check_data()
             invisible(self)
           },
           initialize = function(n = 100, p = 4,
                                 param = list(), 
                                 pscore_design = "A",
                                 outcome_design = "A",
                                 overlap = "high", 
                                 treatment_effect = "A", ...) {
             
             
             if (isTRUE(missing(p)) | isTRUE(is.null(p))  | p < 4) {
               warning("'p' must be greater than of equal to 4. Set to 4.")
               private$p <- 4
             } else {
               private$p <- as.integer(p)
             }
             
             if (isTRUE(missing(n)) || isTRUE(is.null(n))) {
               private$n <- 512
             } else {
               private$n <- n
             }
             if (missing(outcome_design ) || is.null(outcome_design) ) {
               private$outcome_design <- "A"
             } else {
               private$outcome_design <- match.arg(design, c("A","B"))
             }
             if (missing(pscore_design ) || is.null(pscore_design) ) {
               private$pscore_design <- "A"
             } else {
               private$pscore_design <- match.arg(design, c("A","B"))
             }
             if (missing(treatment_effect ) || is.null(treatment_effect) ) {
               private$treatment_effect <- "A"
             } else {
               private$treatment_effect <- match.arg(treatment_effect, c("A","B"))
             }
             private$set_param(beta_z = param$beta_z, 
                               beta_y = param$beta_y,
                               sigma_y = param$sigma_y)
             
           },
           get_design = function() {
             return(c("Outcome Design" = private$outcome_design, "Overlap" = private$overlap))
           }
         ),
         private = list( outcome_design = "character",
                         pscore = "numeric",
                         pscore_design = "character",
                         set_param = function(beta_z, beta_y, sigma_y, sigma_x) {
                           miss.null <- function(xx) {
                             return(missing(xx) | is.null(xx))
                           }
                           extra.terms <- switch(private$overlap,
                                                 "high" = switch(as.character(private$p),
                                                                 "4" = NULL,
                                                                 "20" = c(-.3, .5),
                                                                 c(-.1, .15)),
                                                 "low" = switch(as.character(private$p),
                                                                "4" = NULL,
                                                                "20" = c(-1,1),
                                                                c(-.3, .5)) )
                           default_param <- list(
                             tx_effect = switch(private$treatment_effect,
                                                "A" = function(x){x},
                                                "B" = function(x){private$pscore^2 + 2 * private$pscore + 1},
                                                "C" = function(x){10 * (privatea$pscore - 0.5)^2 + 0.5}),
                             beta_y = switch(as.character(private$p),
                                             "4" = switch(private$outcome_design,
                                                   "A" = c(0.5, 1, .6, 2.2, -1.2),
                                                   "B" = c(0.5, 1, .6, 2.2, -1.2, 0.3)),
                                             "12" = switch(private$outcome_design,
                                                           "A" = c(0.5, 1, .6, 2.2, -1.2, 1, -1, 1, -1),
                                                           "B" = c(0.5, 1, .6, 2.2, -1.2, 1, -1, 1, -1,
                                                                   .5, .6, .5)),
                                             switch(private$outcome_design,
                                                    c(0.5, 1, .6, 2.2, -1.2, 1, -1, 1, -1))
                             ),
                             sigma_y = 1,
                             beta_z = switch(private$p,
                                             "4" = switch(private$pscore_design,
                                                          "A" = switch(private$overlap,
                                                                       "high" = c(-1.5, 0.5, -0.75, 2, -0.5),
                                                                       "low" = c(-1, 0.4, -1.5, 2, -1.5)),
                                                          "B" = switch(private$overlap,
                                                                       "high" = c(-1.5, 0.5, -0.75, 2, -0.5, .43, 1.4, .25, .3),
                                                                       "low" = c(-1, 0.4, -1.5, 2, -1.5, .65, 2, .4, .45))),
                                             # "12" = switch(),
                                             c(-1.5, 0.5, -0.75, 2, -0.5, rep(extra.terms,2))
                             )
                           )
                           temp_param <- list()
                           if (miss.null(beta_z)) {
                             temp_param$beta_z  <- default_param$beta_z
                           } else {
                             stopifnot(is.vector(param$beta_z))
                             temp_param$beta_z  <- param$beta_z
                           }
                           if (miss.null(beta_y)) {
                             temp_param$beta_y  <- default_param$beta_y
                           } else {
                             stopifnot(is.vector(beta_y))
                             temp_param$beta_y  <- param$beta_y
                           }
                           if (miss.null(sigma_y)) {
                             temp_param$sigma_y <- default_param$sigma_y
                           } else {
                             temp_param$sigma_y <- param$sigma_y
                           }
                           if (miss.null(sigma_x)) {
                             temp_param$sigma_x <- default_param$sigma_x
                           } else {
                             temp_param$sigma_x <- param$sigma_x
                           }
                           private$param <- temp_param
                         },
                         treatment_effect = "character"
         )
  )



  HulingMak2020_univariate <- R6::R6Class("HulingMak2020_univariate", 
                             inherit = DataSim,
                             public = list(
                               gen_data = function() {
                                 self$gen_x()
                                 self$gen_z()
                                 self$gen_y()
                                 # private$\check_data()
                                 invisible(self)
                               },
                               gen_x = function() {
                                 stopifnot(length(private$n) > 0 )
                                 private$x <- as.matrix(sort(rnorm(private$n)))
                                 colnames(private$x) <- "X"
                                 private$xdensity <- pnorm(private$x)
                                 private$check_data()
                                 invisible(self)
                               },
                               gen_y = function() {
                                 if (all(is.character(private$x))) self$gen_x()
                                 if (all(is.character(private$z))) self$gen_z()
                                 
                                 mean_y <- private$x + private$x^3 - 1/(0.1 + 0.1* private$x^2)
                                 private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = sqrt(2)))
                                 # private$check_data()
                                 invisible(self)
                               },
                               gen_z = function() {
                                 if (all(is.character(private$x))) self$gen_x()
                                 
                                 latent_z <- switch(private$design,
                                                    "A" = private$x - 1,
                                                    "B" = 2 /3  * private$x^2 + private$x - 1,
                                                    "C" = -1 / 3 * private$x^3 + 2 /3  * private$x^2 + private$x - 1)
                                 private$pscore <- plogis(latent_z)
                                 private$z <- rbinom(private$n, 
                                                     size = 1, 
                                                     prob = private$pscore)
                                 private$check_data()
                                 private$tx_density <- private$pscore * private$xdensity
                                 private$cn_density <- (1 - private$pscore) * private$xdensity
                                 
                                 invisible(self)
                               },
                               get_density = function(which = c("full","pscore", "treated", "control")) {
                                 which.dens <- match.arg(which)
                                 return( switch(which.dens,
                                                "full" = private$xdensity,
                                                "pscore" = private$pscore,
                                                "treated" = private$tx_density,
                                                "control" = private$cn_density) )
                               },
                               initialize = function(n = 128, design = "A", ...) {
                                 
                                 
                                 if (missing(n) | is.null(n)) {
                                   private$n <- 128
                                 } else {
                                   private$n <- n
                                 }
                                 if (missing(design ) | is.null(design) ) {
                                   private$design <- "A"
                                 } else {
                                   private$design <- match.arg(design, c("A","B","C"))
                                 }
                                 private$p <- 1
                               },
                               get_design = function() {
                                 return(c(design = private$design))
                               }
                             ),
                             private = list(design = "character",
                                            overlap = "character",
                                            pscore = "numeric",
                                            xdensity = "numeric",
                                            tx_density = "numeric",
                                            cn_density = "numeric")
  )

