testthat::test_that("cost rkhs works", {
  testthat::skip_on_cran()
  set.seed(008868)
  
  n <- 2^6
  p <- 6
  nsims <- 2
  overlap <- "high"
  design <- "A"
  distance <- c("RKHS")
  power <- c(1,2)
  ground_power <- 1:2
  std_mean_diff <- c(0.001, 0.01, 0.1)
  solver <- "mosek"
  estimand <- c("ATT", "ATC","ATE")
  
  #### get simulation functions ####
  simulator <- Hainmueller$new(n = n, p = p, 
                              design = design, overlap = overlap)
  simulator$gen_data()
  RKHS.opt <- lapply(estimand, function(est)
    RKHS_param_opt(x = simulator$get_x(),
                 z = simulator$get_z(),
                 y = simulator$get_y(),
                 power = 2:3,
                 metric = "mahalanobis",
                 is.dose = FALSE,
                 opt.method = "stan",
                 estimand = est))
  names(RKHS.opt) <- estimand
  x <- simulator$get_x()
  z <- simulator$get_z()
  
  
  testthat::expect_silent(
    costs <- lapply(estimand, function(est) 
    cost_fun(x=x, z=z,
           power = 2, 
           metric = "RKHS",
           rkhs.args = RKHS.opt[[est]],
           estimand = est)
  )
  )
  
})
