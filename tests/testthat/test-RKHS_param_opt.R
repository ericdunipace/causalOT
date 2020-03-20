testthat::test_that("parameter optimization for RKHS", {
  
  set.seed(290384)
  library(causalOT)
  
  n <- 2^9
  p <- 6
  overlap <- "low"
  design <- "A"
  distance <- c("mahalanobis")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  x <- data$get_x()
  y <- data$get_y()
  z <- data$get_z()
  
  # debugonce("RKHS_param_opt")
  
  opt <- RKHS_param_opt(x, y, z, power = 2:3,
                 metric = c("mahalanobis", "Lp"), is.dose = FALSE, 
                         opt.method = c("bayesian.optimization"),
                 bounds = list(theta_0 = c(0,5),
                      theta_1 = c(.Machine$double.xmin,5),
                      gamma_0 = c(.Machine$double.xmin,5),
                      gamma_1 = c(.Machine$double.xmin,5),
                      sigma2 = c(.Machine$double.xmin,1)),
                 initPoints = 10,
                 iters.k = 1,
                 iters.n = 10
                 ) 
  
  # debugonce("RKHS_param_opt")
  
  opt <- RKHS_param_opt(x, y, z, power = 2:3,
                        metric = c("mahalanobis"), is.dose = FALSE, 
                        opt.method = c("optim"), control = list(maxit = 1000))
                        
  # debugonce("RKHS_param_opt")
  
  opt <- RKHS_param_opt(x, y, z, p = 2:3,
                        metric = c("mahalanobis"), is.dose = FALSE, 
                        opt.method = c("stan")) 
})
