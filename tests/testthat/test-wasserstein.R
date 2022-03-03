testthat::test_that("runs for current causal weights", {
  testthat::skip_on_cran()
  set.seed(23483)
  n <- 2^7
  p <- 6
  power <- 2
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(4)
  solver <- "osqp"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  
  cost <- causalOT:::cost_mahalanobis(data$get_x0(), data$get_x1(), power)
  
  weights <- lapply(estimates, function(e) causalOT::calc_weight(data = data, 
                                                       constraint = 4, 
                                                       estimand = e, 
                                                       method = "Wasserstein",
                                                       solver = solver,
                                                       p = power,
                                                       cost = cost))
  
  b <- rep(1/n0,n0)
  
  testthat::expect_silent(dists <- sapply(weights, function(w) causalOT:::wasserstein_p(a = w, p = power, cost = cost, estimand = w$estimand)))
  for(i in seq_along(dists)) {
    testthat::expect_equivalent(dists[i], approxOT::wasserstein(a = weights[[i]]$w0,
                                                            b = weights[[i]]$w1,
                                                            p = power,
                                                            cost = cost),
                                tolerance = 1e-3)
  }
  
  
})
