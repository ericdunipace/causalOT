testthat::test_that("barycenter projection pow2, lp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  power <- 2
  nsims <- 1
  overlap <- "low"
  design <- "A"
  distance <- c("Lp")
  power <- c(2)
  solver <- "osqp"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metric = "Lp"))
  
  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 3, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       solver = solver,
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE)))
  
  
  testthat::expect_equal(c(crossprod(weights[[1]]$gamma, y0) * 1/colSums(weights[[1]]$gamma)), 
                         barycenter_pow2(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(weights[[1]]$gamma %*% y1) * 1/rowSums(weights[[1]]$gamma) , 
  barycenter_pow2(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
  testthat::expect_equal(c(crossprod(weights[[2]]$gamma, y0) * 1/colSums(weights[[2]]$gamma)), 
                         barycenter_pow2(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(weights[[2]]$gamma %*% y1) * 1/rowSums(weights[[2]]$gamma), 
                         barycenter_pow2(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
  
})

testthat::test_that("barycenter projection pow2, mahalanobis", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  power <- 2
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp")
  power <- c(2)
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metic = "mahalanobis"))  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 1.7, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       solver = "osqp",
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE)))
  
  
  testthat::expect_equal(c(crossprod(weights[[1]]$gamma, y0)*n1), 
                         barycenter_pow2(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(weights[[1]]$gamma %*% y1 * 1/rowSums(weights[[1]]$gamma)), 
                         barycenter_pow2(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
  testthat::expect_equal(c(crossprod(weights[[2]]$gamma, y0)*1/colSums(weights[[2]]$gamma)), 
                         barycenter_pow2(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(weights[[2]]$gamma %*% y1)*1/rowSums(weights[[2]]$gamma), 
                         barycenter_pow2(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
})

testthat::test_that("barycenter projection pow1, mahalanobis", {
  testthat::skip_if_not_installed(pkg = "osqp")
  set.seed(23483)
  n <- 2^7
  p <- 6
  power <- 1
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp")
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  e <- eigen(0.5*(cov(cbind(x0,y0)) + cov(cbind(x1,y1))))
  U <- e$vectors %*% diag(sqrt(abs(e$values))) %*% t(e$vectors)
  U_inv <- solve(U)
  dt0 <- (cbind(x0,y0) %*% U_inv)
  dt1 <- (cbind(x1,y1) %*% U_inv)
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metic = "mahalanobis")) 
  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 3, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE)))
  
  
  testthat::expect_equal(c((t(U) %*% apply(weights[[1]]$gamma,2,function(w) matrixStats::colWeightedMedians(x=dt0, w=w)))[p+1,]), 
                         barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "mahalanobis")$y0)
  testthat::expect_equal(c((t(U) %*% apply(weights[[1]]$gamma,1,function(w) matrixStats::colWeightedMedians(x=dt1, w=w)))[p+1,]), 
                         barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "mahalanobis")$y1)
  
  testthat::expect_equal(c((t(U) %*% apply(weights[[2]]$gamma,2,function(w) matrixStats::colWeightedMedians(x=dt0, w=w)))[p+1,]), 
                         barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "mahalanobis")$y0)
  testthat::expect_equal(c((t(U) %*% apply(weights[[2]]$gamma,1,function(w) matrixStats::colWeightedMedians(x=dt1, w=w)))[p+1,]), 
                         barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "mahalanobis")$y1)
  
  # library(ggridges)
  # fit0 <-lm(y0 ~ x0)
  # df <- rbind(data.frame(y = y0, method = "control"),
  #             data.frame(y = y1, method = "counterfactual"),
  #             data.frame(y = barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "mahalanobis")$y0, method = "barycenter_mahal"),
  #             data.frame(y = barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0, method = "barycenter_lp"),
  #             data.frame(y = cbind(1,x1) %*% fit0$coef, method = "regression"))
  # 
  # ggplot(df, aes(x = y, y = method)) + geom_density_ridges()
  # mean((barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "mahalanobis")$y0-
  #        y1)^2)
  # mean((barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0-
  #         y1)^2)
  # 
  # mean(((cbind(1,x1) %*% fit0$coef)-
  #         y1)^2)
  # 
  # 
  # fit1 <-lm(y1 ~ x1)
  # df <- rbind(data.frame(y = y1, method = "treated"),
  #             data.frame(y = y0, method = "counterfactual"),
  #             data.frame(y = barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "mahalanobis")$y1, method = "barycenter_mahal"),
  #             data.frame(y = barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1, method = "barycenter_lp"),
  #             data.frame(y = cbind(1,x0) %*% fit1$coef, method = "regression"))
  # 
  # ggplot(df, aes(x = y, y = method)) + geom_density_ridges()
  # 
  # mean((barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "mahalanobis")$y1-
  #         y0)^2)
  # mean((barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1-
  #         y0)^2)
  # 
  # mean((cbind(1,x0) %*% fit1$coef-
  #         y0)^2)
})

testthat::test_that("barycenter projection pow1, lp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  power <- 1
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp")
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metic = "Lp")) 
  
  testthat::expect_warning(weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 3, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE)))
  
  
  testthat::expect_equal(c(apply(weights[[1]]$gamma,2,function(w) matrixStats::weightedMedian(x=y0, w=w))), 
                         barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(apply(weights[[1]]$gamma,1,function(w) matrixStats::weightedMedian(x=y1, w=w))), 
                         barycenter_pow1(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
  testthat::expect_equal(c(apply(weights[[2]]$gamma,2,function(w) matrixStats::weightedMedian(x=y0, w=w))), 
                         barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp")$y0)
  testthat::expect_equal(c(apply(weights[[2]]$gamma,1,function(w) matrixStats::weightedMedian(x=y1, w=w))), 
                         barycenter_pow1(weights[[2]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATC", metric = "Lp")$y1)
  
})

testthat::test_that("barycenter projection pow3, mahalanobis", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("mahalanobis")
  power <- 3
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metic = "mahalanobis")) 
  
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 1.7, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE))
  
  # debugonce(barycenter_estimation)
  testthat::expect_silent(barycenter_estimation(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = distance, power = power)$y0)
  
  
})

testthat::test_that("barycenter projection pow4, lp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp")
  power <- 4
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  
  cost <- lapply(estimates, function(e)
    cost_fun(data$get_x(), data$get_z(), power = power,
             estimand = e, metic = "Lp"))   
  weights <- lapply(estimates, function(e) calc_weight(data = data, 
                                                       constraint = 1.7, 
                                                       estimand = e, 
                                                       method = "NNM",
                                                       p = power,
                                                       cost = cost[[e]],
                                                       transport.matrix = TRUE))
  
  # debugonce(barycenter_estimation)
  testthat::expect_silent(barycenter_estimation(weights[[1]]$gamma, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp", power = power)$y0)

})

testthat::test_that("barycenter projection pow4, lp", {
  set.seed(23483)
  n <- 2^7
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("Lp")
  power <- 4
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  
  # cost <- causalOT::cost_calc_lp(data$get_x0(), data$get_x1(), power)
  
  gamma1 <- matrix(rep(1/n1/n0,n0),n0,n1)
  gamma0 <- matrix(rep(1/n0/n1,n1),n0,n1, byrow = TRUE)
  
  # debugonce(barycenter_estimation)
  testthat::expect_equal(16.118159135815048444,
                         barycenter_estimation(gamma1, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "Lp", power = power)$y0[1],
                         check.attributes = FALSE)
  
  
})

testthat::test_that("barycenter projection pow4, mahalanobis", {
  set.seed(23483)
  n <- 2^6
  p <- 6
  nsims <- 1
  overlap <- "low"
  design <- "B"
  distance <- c("mahalanobis")
  power <- 4
  solver <- "mosek"
  estimates <- c("ATT", "ATC", "cATE", "ATE")
  
  #### get simulation functions ####
  data <- causalOT::Hainmueller$new(n = n, p = p, 
                                    design = design, overlap = overlap)
  data$gen_data()
  ns <- data$get_n()
  n0 <- ns["n0"]
  n1 <- ns["n1"]
  x0 <- data$get_x0()
  x1 <- data$get_x1()
  z <- data$get_z()
  y <- data$get_y()
  y0 <- y[z==0]
  y1 <- y[z==1]
  
  
  
  gamma1 <- matrix(rep(1/n1/n0,n0),n0,n1)
  gamma0 <- matrix(rep(1/n0/n1,n1),n0,n1, byrow = TRUE)
  
  testthat::expect_equal(11.13143461, #12.623257847612290306,
                         barycenter_estimation(gamma1, x0 = x0, x1 = x1, y0 = y0, y1 = y1, estimand = "ATT", metric = "mahalanobis", power = power)$y0[1],
                         check.attributes = FALSE)
  
  
})


