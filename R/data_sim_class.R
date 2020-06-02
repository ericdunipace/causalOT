DataSim <- R6::R6Class("DataSim",
        public = list(get_x = function() { return( private$x)},
                      get_y = function() { return( private$y)},
                      get_z = function() { return( private$z)},
                      get_n = function() {return(c("n0" = private$n0, "n1" = private$n1))},
                      get_x1 = function() {
                        if(!is.character(private$x)) {
                          if(is.character(private$x1)) private$x1 <- private$x[private$z == 1,,drop=FALSE]
                          return(private$x1)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      get_x0 = function(){
                        if(!is.character(private$x)) {
                          if(is.character(private$x0)) private$x0 <- private$x[private$z == 0,,drop=FALSE]
                          return(private$x0)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      get_p = function() {
                        return(private$p)
                      },
                      gen_data = function(){NULL}),
         private = list(n = "numeric",
                       p = "numeric",
                       x = "matrix",
                       y = "vector",
                       z = "vector",
                       param = "list",
                       n1 = "numeric",
                       n0 = "numeric",
                       x0 = "vector",
                       x1 = "vector",
                       check_data = function() {
                         complete <- all(is.matrix(private$x) & is.vector(private$z) )
                         if(complete) {
                           private$n1 <- sum(private$z == 1)
                           private$n0 <- sum(private$z == 0)
                           private$x1 <- private$x[private$z == 1,,drop=FALSE]
                           private$x0 <- private$x[private$z == 0,,drop=FALSE]
                         }
                       }#,
                       )
)

Hainmueller <- R6::R6Class("Hainmueller", 
                           inherit = DataSim,
                           public = list(
                             gen_data = function() {
                               self$gen_x()
                               self$gen_z()
                               self$gen_y()
                               # private$check_data()
                               invisible(self)
                             },
                             gen_x = function() {
                               stopifnot(length(private$n) >0 )
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
                             gen_y = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               mean_y <- if(private$design =="A" | private$design == 1) {
                                 private$x %*% private$param$beta_y
                               } else if (private$design =="B" | private$design == 2) {
                                 (private$x[,c(1,2,5)] %*% private$param$beta_y[1:3])^2
                               } else {
                                 stop("design must be one of 'A' or 'B'")
                               }
                               private$y <- c(mean_y + rnorm(private$n, mean = 0, sd = private$param$sigma_y))
                               # private$check_data()
                               invisible(self)
                             },
                             gen_z = function() {
                               if(all(dim(private$x) == 0)) gen_x()
                               mean_z <- private$x %*% private$param$beta_z
                               latent_z <- mean_z + rnorm(private$n, mean=0, sd = private$param$sigma_z)
                               private$z <- c(ifelse(latent_z < 0, 0, 1))
                               private$check_data()
                               invisible(self)
                             },
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
                                 private$design <- match.arg(design, c("A","B"))
                               }
                               if( missing(overlap) | is.null(overlap) ) {
                                 private$overlap <- "low"
                               } else {
                                 private$overlap <- match.arg(overlap, c("low","high"))
                               }
                               private$
                                 set_param(beta_z = param$beta_z, beta_y = param$beta_y,
                                           sigma_z = param$sigma_z, sigma_y = param$sigma_y,
                                           param_x = param$param_x)
                               
                             },
                             get_design = function() {
                               return(c(design = private$design, overlap = private$overlap))
                             },
                             opt_weight = function(estimand = "ATE", augment = FALSE, solver = "mosek") {
                               if(estimand == "cATE") estimand <- "ATE"
                               estimand <- match.arg(estimand, choices = c("ATT","ATC","ATE"))
                               aug <- isTRUE(augment)
                               solver <- match.arg(solver, choices = c("mosek", "gurobi", "cplex"))
                               
                               # f = switch(as.character(aug),
                               #            "FALSE" = private$y,
                               #            "TRUE" = switch(as.character(private$design),
                               #                       "A" = private$x %*% private$param$beta_y,
                               #                       "1" = private$x %*% private$param$beta_y,
                               #                       "B" = (private$x[,c(1,2,5)] %*% private$param$beta_y[1:3])^2,
                               #                       "2" = (private$x[,c(1,2,5)] %*% private$param$beta_y[1:3])^2,
                               #                       stop("design must be one of 'A' or 'B'")))
                               m0 <- .lm.fit(private$x[private$z==0,,drop = FALSE], private$y[private$z==0,drop=FALSE])
                               m1 <- .lm.fit(private$x[private$z==1,,drop = FALSE], private$y[private$z==1,drop=FALSE])
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
                                                               "ATT" = mean(private$x[private$z==1,,drop= FALSE] %*% m0$coefficients),
                                                               "ATC" = mean(private$x[private$z==0,,drop= FALSE] %*% m1$coefficients),
                                                               "ATE" = c(mean(private$x[private$z==1,,drop= FALSE] %*% m0$coefficients),
                                                                         mean(private$x[private$z==0,,drop= FALSE] %*% m1$coefficients))))
                               
                               mu_1 = switch(paste0(c(private$design, private$overlap), collapse =", "),
                                             "A, high" = 2.0944332,
                                             "A, low" = 2.3886758,
                                             "B, high" = 9.003046,
                                             "B, low" = 9.497603 )
                               mu_0 = switch(paste0(c(private$design, private$overlap), collapse =", "),
                                             "A, high" = 0.9165373,
                                             "A, low" = 0.6069783,
                                             "B, high" = 7.011926,
                                             "B, low" = 6.528442)
                               
                               mu = switch(as.character(private$design),
                                           "A" = 1.5,
                                           "1" = 1.5,
                                           "B" = 8,
                                           "2" = 8,
                                           stop("design must be one of 'A' or 'B'"))
                               
                               mu_var <- switch(estimand,
                                      "ATE" = mu,
                                      "ATT" = mu_1,
                                      "ATC" = mu_0)
                               
                               # f <- switch(estimand,
                               #             "ATE" = f,
                               #             "ATT" = f[private$z == 0],
                               #             "ATC" = f[private$z == 1])
                               if(estimand != "ATE" ) {
                                 Q <-Matrix::Matrix(tcrossprod(f), sparse = TRUE)
                                 L <- -2 * (mu_var - const) * t(f)
                                 A <- Matrix::sparseMatrix(i = rep(1, length(f)),
                                                           j = 1:length(f),
                                                           x = 1)
                                 problem <- list(obj = list(Q = Q, L = L),
                                                 LC = list(dir = "E",
                                                           vals = 1,
                                                           A = A))
                                 
                                 weights <- switch(solver,
                                                   "mosek" =  mosek_solver(problem),
                                                   "gurobi" = gurobi_solver(problem),
                                                   "cplex"  = cplex_solver(problem)
                                 )
                                 weights <- renormalize(weights)
                               } else {
                                 f1 <- f[[2]]
                                 f0 <- f[[1]]
                                 
                                 Q0 <- Matrix::Matrix(tcrossprod(f0), sparse = TRUE)
                                 L0 <- -2 * (mu_var - const[1]) * t(f0)
                                 A0 <- Matrix::sparseMatrix(i = rep(1, length(f0)),
                                                           j = 1:length(f0),
                                                           x = 1)
                                 problem0 <- list(obj = list(Q = Q0, L = L0),
                                                 LC = list(dir = "E",
                                                           vals = 1,
                                                           A = A0))
                                 
                                 Q1 <- Matrix::Matrix(tcrossprod(f1), sparse = TRUE)
                                 L1 <- -2 * (mu_var - const[2]) * t(f1)
                                 A1 <- Matrix::sparseMatrix(i = rep(1, length(f1)),
                                                            j = 1:length(f1),
                                                            x = 1)
                                 problem1 <- list(obj = list(Q = Q1, L = L1),
                                                 LC = list(dir = "E",
                                                           vals = 1,
                                                           A = A1))
                                 
                                 weights <- list(w0 = switch(solver,
                                                   "mosek" =  mosek_solver(problem0),
                                                   "gurobi" = gurobi_solver(problem0),
                                                   "cplex"  = cplex_solver(problem0)
                                 ),
                                             w1 = switch(solver,
                                                    "mosek" =  mosek_solver(problem1),
                                                    "gurobi" = gurobi_solver(problem1),
                                                    "cplex"  = cplex_solver(problem1)
                                 ) )
                                 weights <- lapply(weights, renormalize)
                               }
                               return(weights)
                               
                             },
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
                           private = list(design = "character",
                                          overlap = "character",
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
                                                            B = c(1,1,1)),
                                              sigma_z= list(low = sqrt(30),
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
                                            if(is.null(private$design) ) {
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
