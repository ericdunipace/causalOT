cost_fun_deprecated <- function(x, y, power = 2, metric = c("mahalanobis","Lp","RKHS"), 
                                rkhs.args = NULL, estimand = "ATE", ...) {
  direction <- "rowwise"
  metric <- match.arg(metric)
  
  dist <- switch(metric, 
                 "Lp" = causalOT:::cost_calc_lp(x,y,ground_p = power, direction = direction),
                 "mahalanobis" = causalOT:::cost_mahalanobis(x,y, ground_p = power, direction = direction),
                 "RKHS" = causalOT:::cost_RKHS(X=x, Y=y, rkhs.args = rkhs.args, estimand = estimand, ...))
  
  return(dist)
  
}

testthat::test_that("cost function gives correct values", {
  n0 <- 100
  n1 <- 55
  d <- 5
  
  x1 <- matrix(rnorm(n1*d), n1, d)
  x0 <- matrix(rnorm(n0*d), n0, d)
  
  x <- rbind(x0,x1)
  z <- c(rep(0,n0), rep(1,n1))
  
 
  power <- 2.0
  
  #ATT Lp
  estimand <- "ATT"
  metric <- "Lp"
  cost.old <- cost_fun_deprecated(x0, x1, power, metric, rkhs.args = NULL)
  cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  
  testthat::expect_equal(cost.old, cost.new)
  
  #ATT mahal
  metric <- "mahalanobis"
  cost.old <- cost_fun_deprecated(x0, x1, power, metric, rkhs.args = NULL)
  cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  
  testthat::expect_equal(cost.old, cost.new)
  
  estimand <- "ATC"
  metric <- "Lp"
  cost.old <- cost_fun_deprecated(x0, x1, power, metric, rkhs.args = NULL)
  cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  
  testthat::expect_equal(cost.old, cost.new)
  metric <- "mahalanobis"
  cost.old <- cost_fun_deprecated(x0, x1, power, metric, rkhs.args = NULL)
  cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  
  testthat::expect_equal(cost.old, cost.new)
  
  
  estimand <- "ATE"
  metric <- "Lp"
  cost.old <- list(cost_fun_deprecated(x0, x, power, metric, rkhs.args = NULL),
                   cost_fun_deprecated(x1, x, power, metric, rkhs.args = NULL))
  cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  
  testthat::expect_equal(cost.old, cost.new)
  
  metric <- "mahalanobis" # need to update test
  # cost.old <- list(cost_fun_deprecated(x0, x, power, metric, rkhs.args = NULL),
  #                  cost_fun_deprecated(x1, x, power, metric, rkhs.args = NULL))
  # cost.new <- cost_fun(x, z, power = power, metric = metric, estimand = estimand)
  # 
  # testthat::expect_equal(cost.old, cost.new)
})
