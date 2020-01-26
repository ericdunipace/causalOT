Hainmueller <- R6::R6Class("Hainmueller", 
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
        x13 <- matrix(private$param$param_x$x_13$mean, nrow = private$n,
                      ncol = 3, byrow = TRUE) + 
          matrix(rnorm(private$n * 3), 
                 nrow = private$n, 
                 ncol = 3) %*% chol(private$param$param_x$x_13$covar)
        x4 <- runif(private$n, private$param$param_x$x4$lower, private$param$param_x$x4$upper)
        x5 <- rchisq(private$n, df = private$param$param_x$x5$df)
        x6 <- rbinom(private$n, size = 1, prob =private$param$param_x$x6$p)
        
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
           private$check_data()
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
                    if(is.null(overlap) & (miss.null(sigma_z) )) {
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

# hainmueller$lock(c("n","p","param","design","overlap"))

# hainmueller_sim <- function() {
#   gen_x <- function(n, param_x) {
#     x13 <- param_x$x_13$mean + matrix(rnorm(n * 3), nrow = n, ncol = 3) %*% chol(param_x$x_13$covar)
#     x4 <- runif(n, param_x$x4$lower, param_x$x4$upper)
#     x5 <- rchisq(n, df = param_x$x5$df)
#     x6 <- rbinom(n, size = 1, prob = param_x$x6$p)
#     
#     x <- cbind(x13, x4, x5, x6)
#     colnames(x) <- 1:6
#     return(x)
#   }
#   
#   gen_y <- function(n, x, beta_y, sigma_y, design = "A") {
#     mean_y <- if(design =="A" | design == 1) {
#       x %*% beta_y
#     } else if (design =="B" | design == 2) {
#       (x[,c(1,2,5)] %*% beta_y[1:3])^2
#     } else {
#       stop("design must be one of 'A' or 'B'")
#     }
#     y <- mean_y + rnorm(n, mean = 0, sd = sigma_y)
#     return(y)
#   }
#   
#   gen_z <- function(n, x, beta_z, sigma_z) {
#     mean_z <- x %*% beta_z
#     latent_z <- mean_z + rnorm(n, mean=0, sd = sigma_z)
#     z <- ifelse(latent_z <0, 0, 1)
#     return(z)
#   }
#   
#   default_param <- list(
#     beta_z = c(1,2,-2,-1,-0.5,1),
#     beta_y = list(A = c(1,1,1,-1,1,1),
#                   B = c(1,1,1)),
#     sigma_z= list(low_overlap = sqrt(30),
#                   high_overlap = sqrt(100)),
#     sigma_y = 1,
#     param_x = list(x_13 = list(mean = rep(0, 3),
#                                 covar = matrix(c(2,1,-1,1,1,-0.5, -1, -0.5, 1), nrow = 3,ncol = 3)),
#                     x4 = list(lower = -3, upper = 3),
#                     x5 = list(df = 1),
#                     x6 = list(p = 0.5))
#   )
#   return(list(
#     rY = gen_y,
#     rX = gen_x,
#     rZ = gen_z,
#     default_param = default_param
#   ))
# }