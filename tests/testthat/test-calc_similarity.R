testthat::test_that("similarity calcs work, non dose, ATE", {
  set.seed(989080)
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
  
  b <- scale(x, center = TRUE, scale=FALSE)
  a <- b %*% solve(chol(cov(x)))
  
  
  # debugonce("calc_similarity")
  sim <- calc_similarity(X = x, z = z, metric = distance, kernel = "polynomial",
                         is.dose = FALSE,
                         estimand = "ATE")
  testthat::expect_equal(sim, tcrossprod(a))
  
  testthat::expect_equal(tcrossprod(x), calc_similarity(X = x, z = z, metric = "Lp", 
                                                        kernel = "polynomial",
                                                        is.dose = FALSE,
                                                        estimand = "ATE"))
})

testthat::test_that("similarity calcs work, non dose, ATC", {
  set.seed(989080)
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
  
  mean <- colMeans(x)
  cov  <- cov(x)
  
  b <- scale(x, scale = FALSE, center = TRUE)
  a <- b %*% solve(chol(cov))
  
  
  # debugonce("calc_similarity")
  sim <- calc_similarity(X = x, z = z, metric = distance, 
                         kernel = "polynomial", 
                         is.dose = FALSE,
                         estimand = "ATC")
  testthat::expect_equal(sim, tcrossprod(a))
  
  testthat::expect_equal(tcrossprod(x), calc_similarity(X = x, z = z, metric = "Lp", 
                                                        kernel = "polynomial",
                                                        is.dose = FALSE,
                                                        estimand = "ATC"))
})

testthat::test_that("similarity calcs work, non dose, ATT", {
  set.seed(989080)
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
  
  mean <- colMeans(x)
  cov  <- cov(x)
  
  b <- x - matrix(mean, n,p, byrow = TRUE)
  a <- b %*% solve(chol(cov))
  
  
  # debugonce("calc_similarity")
  sim <- calc_similarity(X = x, z = z, metric = distance, 
                         kernel = "polynomial",
                         is.dose = FALSE,
                         estimand = "ATT")
  testthat::expect_equal(sim, tcrossprod(a))
  
  testthat::expect_equal(tcrossprod(x), calc_similarity(X = x, z = z, 
                                                        kernel = "polynomial",
                                                        metric = "Lp", is.dose = FALSE,
                                                        estimand = "ATT"))
})

testthat::test_that("similarity calcs work, dose", {
  set.seed(979797)
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
  
  b <- scale(x, scale=FALSE)
  a <- b %*% solve(chol(cov(x)))
  
  z_b <- scale(z, scale= FALSE)
  z_a <- scale(z)
  
  # debugonce("calc_similarity")
  sims <- calc_similarity(X = x, z = z, kernel = "polynomial",
                          metric = distance, is.dose = TRUE)
  testthat::expect_equal(sims[["X"]], tcrossprod(a))
  testthat::expect_equal(tcrossprod(z_a), sims[["Z"]])
  
  
  sims2 <- calc_similarity(X = x, z = z, metric = "Lp", is.dose = TRUE)
  testthat::expect_equal(tcrossprod(b), sims2[["X"]])
  testthat::expect_equal(tcrossprod(z_b), sims2[["Z"]])
  
  testthat::expect_warning(calc_similarity(X = x, z = z, 
                                           kernel = "polynomial",
                                           metric = "Lp", is.dose = TRUE, estimand = "ATT"))
})
