set.seed(908709)

#### Load Packages ####
library(ROI)
library(ROI.plugin.ecos)
library(ROI.plugin.qpoases)
library(ROI.plugin.cplex)
library(ggplot2)
library(slam)
library(limbs)
library(causalOT)
library(doParallel)

#### Data Dim ####
n <- 2^9
p <- 6
nsims <- 100
power <- 1

#### data Gen Fun ####
gen_x <- function(n, param_x) {
  x13 <- param_x$x_13$mean + matrix(rnorm(n * 3), nrow = n, ncol = 3) %*% chol(param_x$x_13$covar)
  x4 <- runif(n, param_x$x4$lower, param_x$x4$upper)
  x5 <- rchisq(n, df = param_x$x5$df)
  x6 <- rbinom(n, size = 1, prob = param_x$x6$p)
  
  x <- cbind(x13, x4, x5, x6)
  colnames(x) <- 1:6
  return(x)
}

gen_y <- function(n, x, beta_y, sigma_y, design = 1) {
  mean_y <- if(design ==1) {
    x %*% beta_y
  } else if (design == 2) {
    (x[,c(1,2,5)] %*% beta_y[1:3])^2
  } else {
    stop("design must be in c(1,2)")
  }
  y <- mean_y + rnorm(n, mean = 0, sd = sigma_y)
  return(y)
}

gen_z <- function(n, x, beta_z, sigma_z) {
  mean_z <- x %*% beta_z
  latent_z <- mean_z + rnorm(n, mean=0, sd = sigma_z)
  z <- ifelse(latent_z <0, 0, 1)
  return(z)
}


#### Bal weight fun, std dif mean ####
qp_sbw <- function(x, z, K) {
  stopifnot(length(K) == ncol(x))
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  x1_var <- apply(x1,2,var)
  x0_var <- apply(x0,2,var)
  pool_sd <- sqrt(mean(x1_var + x0_var))
  # pool_sd <- 1
  n <- nrow(x0)
  d <- ncol(x0)
  
  x_constraint <- t(x0)
  x1m <- colMeans(x1)
  # f_funct <- function(w) {
  #   diff <- crossprod(x0[,1,drop=FALSE], w) - x1m
  #   
  #   return(abs(diff/pool_sd))
  # }
  
  Konst <- K * pool_sd
  
  Q0 <- 2 * (diag(1, n, n) - matrix(1/n, n,n)) # *2 since the Q part is 1/2 a^\top (x^\top x) a in ROI!!!
  L0 <- c(w = rep(0,n))
  var_bounds <- ROI::V_bound(li =1:n, ui = 1:n, lb = rep(0,n), ub = rep(1,n))
  var_bounds <- NULL
  op <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = L0, names = as.character(1:n)),
                bounds = var_bounds,
                maximum = FALSE)
  
  A1 <- rep(1,n) # beta portion, gamma portion, t portion is 0. multiplies all vars each time!!!
  A2 <- x_constraint
  
  A <- rbind(A1, A2, A2)
  vals <- c(1, Konst + x1m, - Konst + x1m)
  dir <- c(ROI::eq(1), ROI::leq(d), ROI::geq(d))
  LC <- ROI::L_constraint(A, dir, vals)
  
  # A <- A1
  # vals <- 1
  # dir <- ROI::eq(1)
  # LC <- ROI::L_constraint(A, dir, vals)
  
  # Q1 <- matrix(0,0)
  # L1 <- rep(1,n) # beta portion, gamma portion, t portion is 0. multiplies all vars each time!!!
  # 
  # Q2 <- lapply(1:d, function(i) 2*crossprod(x_constraint[i,,drop=FALSE]))
  # L2 <- lapply(1:d, function(i) 2*crossprod(x_constraint[i,,drop=FALSE]))
  # A <- rbind(A1,A2)
  # vals <- c(1,K + x1m)
  # dir <- c(ROI::eq(1), ROI::leq(d))
  # LC <- ROI::L_constraint(A, dir, vals)
  
  ROI::constraints(op) <- LC
  
  
  return(op)
}

#### Bal weight fun, w2 ####
qp_w2 <- function(x, z, K, p = 2, target = c("control", "treated", 
                                             "feasible control",
                                             "feasbile treated"),
                  dist=c("Lp", "mahalanobis")) {
  
  target <- match.arg(target)
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  x1_var <- cov(x1)
  x0_var <- cov(x0)
  
  dist <- match.arg(dist)
  dist.fun <- switch(dist, 
                     "Lp" = causalOT::cost_calc_lp,
                     "mahalanobis" = causalOT::cost_mahalanobis)
  # covar <- (x1_var + x0_var)/2 # add later
  cost <- dist.fun(x0, x1, ground_p = p)
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
  
  cost_vec <- simple_triplet_matrix(i = rep(1, n0*n1),
                                    j = 1:(n0*n1),
                                    v= c(cost)^p,
                                    nrow = 1, ncol = n0*n1)
  sum_const <- simple_triplet_matrix(i = rep(1, n0*n1),
                                     j = 1:(n0*n1),
                                     v = rep(1, n1 * n0),
                                     nrow = 1, ncol = n0*n1)
  q_mat <- marg_const <- NULL
  
  # Qmat
  if(target == "treated") {
    # q_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1,n0*n1,n0)
    q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
                                j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                x = 1,
                                dims = c(n0,n1*n0))
    # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
    marg_mass <- 1/n0
    marg_n <- n0
    marg_const_mat <- simple_triplet_matrix(i = rep(1:n1, each = n0),
                                            j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
                                            v = rep(1,n0*n1),
                                            nrow = n1, ncol = n0*n1)
    # marg_const_mat <- matrix(0, n1,n1*n0)
    # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
    marg_const <- rep( 1/n1, n1)
    marg_const_n <- n1
  } else if (target == "control") {
    stop("not yet correct")
    q_mat <- matrix(0, n1,n1*n0)
    for(i in 1:n1) q_mat[i, (i-1) * n1 + (1:n0)] <- 1
    # for(i in 1:n1) q_mat[(i-1) * n0 + (1:n0),i] <- 1
    marg_mass <- 1/n1
    marg_n <- n1
    marg_const_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1,n0*n1,n0)
    for(i in 1:n1) marg_const_mat[i, (i-1) + n1_idx] <- 1
    marg_const <- rep( 1/n0, n0 )
    marg_const_n <- n0
  } else if (target == "feasible treated") {
    stop("not yet implemented")
  } else if (target == "feasible control") {
    stop("not yet implemented")
  }
  # mean_mat <- as.simple_triplet_matrix(matrix(marg_mass, marg_n, n1*n0))
  # q_s <- Matrix::sparseMatrix(i = q_mat_list$i,
  #                             j = q_mat_list$j,
  #                             x = q_mat_list$val)  #Matrix::Matrix(q_mat, sparse = TRUE)
  qtq_x2 <-  as( Matrix::crossprod(q_s),"dgTMatrix")
  Q0 <- simple_triplet_matrix(i = qtq_x2@i+1, j = qtq_x2@j+1, v = qtq_x2@x)
  rm(q_s, qtq_x2)
  # Q0 <- simple_triplet_diag_matrix(1, nrow=n0*n1)
  
  L0 <- simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
                              nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  # var_bounds <- ROI::V_bound(li =1:(n0*n1), ui = 1:(n0*n1), lb = rep(0,n0*n1), ub = rep(1,n0*n1))
  var_bounds <- NULL
  op <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = L0),
                bounds = var_bounds,
                maximum = FALSE)
  
  # A1 <- rep(1,n) # beta portion, gamma portion, t portion is 0. multiplies all vars each time!!!
  # A2 <- x_constraint
  
  A <- rbind(cost_vec, sum_const, marg_const_mat)
  vals <- c(K^p, 1, marg_const)
  dir <- c(ROI::leq(1), ROI::eq(1), ROI::eq(marg_const_n))
  LC <- ROI::L_constraint(A, dir, vals)
  
  ROI::constraints(op) <- LC
  
  return(op)
}

update_wp_tol <- function(x, new_val = NULL, p=2) {
  stopifnot(inherits(x, "OP"))
  stopifnot(is.numeric(new_val))
  LC <- ROI::constraints(x)
  LC$rhs[1] <- new_val^p
  ROI::constraints(x) <- LC
  return(x)
}

#### Outcome fun ####
outcome_model_DRH <- function(x,y,z, weights, ...) {
  dots <- list(...)
  est <- dots$est
  if(is.null(est)) est <- "ATE"
  design <- x
  n <- nrow(x)
  # design0 <- design[z==0,,drop=FALSE]
  # design1 <- design[z==1,,drop=FALSE]
  # y0 <- y[z==0]
  # y1 <- y[z==1]
  w0 <- weights$z0
  w1 <- weights$z1
  w <- rep(NA, n)
  w[z==0] <- w0
  w[z==1] <- w1
  w <- w/sum(w)
  
  model0 <- lm(y ~ design, subset = z == 0)
  f0 <- predict(model0, data.frame(design))
  mu_0 <- mean(f0)
  e_0 <- model0$resid #c(y - f0)
  y_0 <- mu_0 + e_0 %*% w0 #(w * (1-z))/sum((w * (1-z)))
  
  if(est == "ATE") {
    model1 <- lm(y ~ design, subset = z == 1)
    f1 <- predict(model1, data.frame(design))
    mu_1 <- mean(f1)
    e_1 <- model1$resid #c(y - f1)
    y_1 <- mu_1 + e_1 %*% w1 #(w * z)/sum(w*z)
  } else if (est == "ATT") {
    y_1 <- mean(y[z==1])
  }
  
  tx_effect <- y_1 - y_0
  
  return(tx_effect)
}

outcome_model_H <- function(x,y,z, weights, ...) {
  dots <- list(...)
  est <- dots$est
  if(is.null(est)) est <- "ATE"
  design <- x
  n <- nrow(x)
  # design0 <- design[z==0,,drop=FALSE]
  # design1 <- design[z==1,,drop=FALSE]
  y0 <- y[z==0]
  y1 <- y[z==1]
  w0 <- weights$z0
  w1 <- weights$z1
  if(is.null(w1)) w1 <- rep(1/n1,n1)
  w <- rep(NA, n)
  w[z==0] <- w0
  w[z==1] <- w1
  w <- w/sum(w)
  w1 <- w1/sum(w1)
  w0 <- w0/sum(w0)
  
  y_0 <- y0 %*% w0
  
  if(est == "ATE") {
    y_1 <- y1 %*% w1
  } else if (est == "ATT") {
    y_1 <- mean(y[z==1])
  }
  
  tx_effect <- y_1 - y_0
  
  return(tx_effect)
}

#### Other functions ####
weighted.var <- function(x, weights) {
  weights <- weights/sum(weights)
  v2 <- sum(weights^2)
  mean <- c(crossprod(x, weights))
  wt.ss <- c(crossprod((x-mean)^2, weights))
  var <- wt.ss/(1-v2)
  return(var)
}
colVar <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,2,var))
  } else {
    return(apply(x, 2, weighted.var, weights = weights))
  }
}
rowVar <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,1,var))
  } else {
    return(apply(x, 1, weighted.var, weights = weights))
  }
}
colSD <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,2,sd))
  } else {
    return(sqrt(apply(x, 2, weighted.var, weights = weights)))
  }
}
rowSD <- function(x, ...) {
  dots <- list(...)
  weights <- dots$weights
  if(is.null(weights) & length(dots)>0) weights <- dots[[1]]
  if(is.null(weights)){
    return(apply(x,1,sd))
  } else {
    return(sqrt(apply(x, 1, weighted.var, weights = weights)))
  }
}

#### Data Param ####
beta_z <- c(1,2,-2,-1,-0.5,1)
beta_y1 <- c(1,1,1,-1,1,1)
beta_y2 <- c(1,1,1)
sigma_z1 <- sqrt(30)
sigma_z2 <- sqrt(100)
sigma_y <- 1
param_x <- list(x_13 = list(mean = rep(0, 3),
                            covar = matrix(c(2,1,-1,1,1,-0.5, -1, -0.5, 1), nrow = 3,ncol = 3)),
                x4 = list(lower = -3, upper = 3),
                x5 = list(df = 1),
                x6 = list(p = 0.5))



cl <- parallel::makeCluster(parallel::detectCores()-1)
registerDoParallel(cl)
# registerDoSEQ()

low_overlap <- foreach(sim = 1:nsims) %dopar% {
  library(ROI)
  library(ROI.plugin.ecos)
  library(ROI.plugin.qpoases)
  library(ROI.plugin.cplex)
  library(ggplot2)
  library(slam)
  library(limbs)
  library(causalOT)
  library(doParallel)
  #### Gen Data ####
  x <- gen_x(n, param_x)
  z <- gen_z(n, x, beta_z, sigma_z1)
  y <- gen_y(n, x, beta_y1, sigma_y, design=1)
  
  x0 <- x[z==0,,drop=FALSE]
  x1 <- x[z==1,, drop = FALSE]
  y1 <- y[z==1]
  y0 <- y[z==0]
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  #### dist fun ####
  dist <- "Lp"
  cost.calc <- switch(dist, "Lp" = causalOT::cost_calc_lp,
                      "mahalanobis" = causalOT::cost_mahalanobis)
  
  #### Original ####
  OP <- qp_sbw(x,z, rep(0.1,6))
  res <- ROI.plugin.cplex:::solve_OP(OP)
  # return(NULL)
  
  # res <- ROI::ROI_solve(CP, solver = "ecos")
  sol <- ROI::solution(res)
  # print(sol)
  weights <- sol[1:sum(1-z)]
  # print(c("sum of unadj wts" = sum(sol[1:sum(1-z)])))
  weights <- weights/sum(weights)
  # ESS <- 1/sum(weights^2)
  # print(c("sum of adj wts" = sum(weights)))
  # print(c("ESS" = ESS))
  # print(c("orig N0" = sum(1-z)))
  
  
  # pool_sd <- sqrt((apply(x[z==0,],2,var) + apply(x[z==1,],2,var))/2)
  
  # print("Std. Means")
  # data.frame(unmatch = abs(colMeans(x[z==0,]) - colMeans(x[z==1,])),
  #       match = abs(t(x[z==0,]) %*% weights - colMeans(x[z==1,])))/pool_sd
  # 
  # print("Var")
  # data.frame(unmatch = abs(colVar(x[z==0,]) - colVar(x[z==1,])),
  #            match = abs(colVar(x[z==0,], weights = weights) - colVar(x[z==1,])))
  
  mass0 <- weights[weights>0]
  mass1 <- rep(1/n1,n1)
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = 2))[1:n0, n0 + (1:n1)]
  cost <- cost.calc(x0,x1, 1)
  cost <- cost[weights > 0,]
  
  # w2_b4 <- transport::wasserstein(rep(1/n0,n0), mass1, p = 2, tplan = NULL, costm = cost)
  w2_match <- transport::wasserstein(mass0, mass1, p = power, tplan = NULL, costm = cost)
  
  # print(c("pre-match" = w2_b4, "After match" = w2_match))
  
  # full_weight <- rep(0, n)
  # full_weight[z==1] <- 1/n1
  # full_weight[z==0] <- weights
  # unmatch_weight <- rep(0,n)
  # unmatch_weight[z==1] <- 1/n1
  # unmatch_weight[z==0] <- 1/n0
  
  # df <- rbind(data.frame(x,z = factor(z), w = unmatch_weight, match = "unmatch"),
  #             data.frame(x,z = factor(z), w = full_weight, match = "match"))
  
  # ggplot(df, aes(x = X1, stat(density), fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1)
  # 
  # ggplot(df, aes(x = X2, stat(density), fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1)
  # 
  # ggplot(df, aes(x = X4, stat(density), fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1)
  # 
  # ggplot(df, aes(x = X5, stat(density), fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1)
  # 
  # ggplot(df, aes(x = X6, stat(density), fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), position = "dodge", alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1)
  
  #### W2 ####
  dist <- "Lp"
  cost.calc <- switch(dist, "Lp" = causalOT::cost_calc_lp,
                      "mahalanobis" = causalOT::cost_mahalanobis)
  OP_w2 <- qp_w2(x, z, 1.3, p = power, target = "treated",dist = dist)
  # OP_w2 <- update_wp_tol(OP_w2, new_val = 1.25)
  res_w2 <- ROI.plugin.cplex:::solve_OP(OP_w2)
  if(all(is.na(ROI::solution(res_w2)))){
    OP_w2 <- update_wp_tol(OP_w2, new_val = 0.9 * w2_match)
    res_w2 <- ROI.plugin.cplex:::solve_OP(OP_w2)
  }
  # print(res_w2)
  
  sol_w2 <- ROI::solution(res_w2)
  w0_w2 <- sol_w2
  w0_w2[w0_w2 < 0] <- 0
  w0_w2 <- w0_w2/sum(w0_w2)
  gamma <- matrix(w0_w2, n0,n1)
  w0_marg_w2 <- rowSums(gamma)#/sum(rowSums(gamma))
  # ESS_w2 <- 1/sum(w0_marg_w2^2)
  
  # print(c("sum of unadj wts" = sum(sol_w2)))
  # print(c("sum of adj wts" = sum(w0_w2)))
  # print(c("ESS" = ESS_w2))
  # print(c("orig N0" = sum(1-z)))
  
  
  
  
  
  # pool_sd <- sqrt((apply(x[z==0,],2,var) + apply(x[z==1,],2,var))/2)
  # 
  # print("Std. Mean")
  # data.frame(unmatch = abs(colMeans(x[z==0,]) - colMeans(x[z==1,])),
  #            sbw_match = abs(t(x[z==0,]) %*% weights - colMeans(x[z==1,])),
  #            match = abs(t(x[z==0,]) %*% w0_marg_w2 - colMeans(x[z==1,])))/pool_sd
  # 
  # print("Var")
  # data.frame(unmatch = abs(colVar(x[z==0,]) - colVar(x[z==1,])),
  #            sbw_match = abs(colVar(x[z==0,], weights = weights) - colVar(x[z==1,])),
  #            w2_match = abs(colVar(x[z==0,], weights = w0_marg_w2) - colVar(x[z==1,])))
  
  mass0 <- weights
  mass0_w2 <- w0_marg_w2
  mass1 <- rep(1/n1,n1)
  cost <- cost.calc(x0,x1, 2) #cost_calc(t(x0),t(x1), 2)
  mass0_rm <- mass0[mass0 > 0]
  mass0_w2_rm <- mass0_w2[mass0_w2>0]
  cost_rm <- cost[mass0_w2 >0,]
  cost_rm_sbw <- cost[mass0>0,]
  
  w2_b4_w2 <- transport::wasserstein(rep(1/n0,n0), mass1, p = 2, tplan = NULL, costm = cost)
  w2_match_sbw <- transport::wasserstein(mass0_rm, mass1, p = 2, tplan = NULL, costm = cost_rm_sbw)
  w2_match_w2 <- transport::wasserstein(mass0_w2_rm, mass1, p = 2, tplan = NULL, costm = cost_rm)
  
  # print(c("pre-match" = w2_b4_w2, "sbw" = w2_b4_sbw, "After match" = w2_match_w2))
  
  
  # full_weight_w2 <- rep(0, n)
  # full_weight_w2[z==1] <- 1/n1
  # full_weight_w2[z==0] <- w0_marg_w2
  # unmatch_weight_w2 <- rep(0,n)
  # unmatch_weight_w2[z==1] <- 1/n1
  # unmatch_weight_w2[z==0] <- 1/n0
  # 
  # df_w2 <- rbind(data.frame(x,z = factor(z), w = unmatch_weight_w2, match = "unmatch"),
  #             data.frame(x,z = factor(z), w = full_weight_w2, match = "match"))
  
  # ggplot(df_w2, aes(x = X1, fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X1")
  # 
  # ggplot(df_w2, aes(x = X2, fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X2")
  # 
  # ggplot(df_w2, aes(x = X4, fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X4")
  # 
  # ggplot(df_w2, aes(x = X5, fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X5")
  # 
  # ggplot(df_w2, aes(x = X6, fill = z, group = z, weight = w)) + 
  #   geom_histogram(aes(fill = z), position = "dodge", alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X6")
  # 
  # ggplot(df_w2, aes(x = X1, fill = z, group = z, weight = w)) + 
  #   geom_density(alpha = 0.5) +
  #   facet_wrap(~match, nrow = 2, ncol = 1) + ggtitle("X1")
  
  # pc <- data.frame(id = 1:n, prcomp(x)$x[,1:2], z= factor(z), unmatch = unmatch_weight, 
  #                  sbw = full_weight, 
  #                  w2 = full_weight_w2)
  # pc_long <- tidyr::gather(pc, "condition", "weight", unmatch:w2, factor_key = TRUE)
  # 
  # ggplot(pc_long, aes(x = PC1, y = PC2, color = z, size = weight)) + 
  #   geom_point(alpha = 0.7) + facet_wrap(~condition) + 
  #   ggthemes::scale_color_tableau() + theme_bw()
  
  
  
  
  #### ATT compare bias ####
  naive_att <- mean(y1) - mean(predict(lm(y0~x0), data.frame(x1)))
  sbw_drh <- outcome_model_DRH(x,y,z, weights = list(z0 = weights, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  sbw_h <- outcome_model_H(x,y,z, weights = list(z0 = weights, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  
  w2_drh <- outcome_model_DRH(x,y,z, weights = list(z0 = w0_marg_w2, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  w2_h <- outcome_model_H(x,y,z, weights = list(z0 = w0_marg_w2, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  
  # print(c("Naive" = naive,
  #         "SBW Hajek" = sbw_h,
  #         "W2 Hajek" =w2_h,
  #         "SBW DR Hajek" = sbw_drh,
  #         "W2 DR Hajek" = w2_drh))
  
  #### ATE calc ####
  OPate <- qp_sbw(x,1-z, rep(0.1,6))
  res_ate <- ROI.plugin.cplex:::solve_OP(OPate)
  
  # res <- ROI::ROI_solve(CP, solver = "ecos")
  sol_ate <- ROI::solution(res_ate)
  weights_ate <- sol_ate[1:sum(z)]
  
  # dist <- "Lp"
  # cost.calc <- switch(dist, "Lp" = cost_calc_lp,
  # "mahalanobis" = cost_mahalanobis)
  OP_w2_ctrl <- qp_w2(x,1-z, 1.3, p = power, target = "treated",dist = dist)
  # OP_w2_ctrl <- update_wp_tol(OP_w2_ctrl, new_val = 1.25)
  res_w2_ctrl <- ROI.plugin.cplex:::solve_OP(OP_w2_ctrl)
  if(all(is.na(ROI::solution(res_w2_ctrl)))){
    OP_w2_ctrl <- update_wp_tol(OP_w2_ctrl, new_val = 0.9 * w2_match)
    res_w2_ctrl <- ROI.plugin.cplex:::solve_OP(OP_w2_ctrl)
  }
  # print(res_w2_ctrl)
  
  sol_w2_ctrl <- ROI::solution(res_w2_ctrl)
  sol_w2_ctrl[sol_w2_ctrl < 0] <- 0
  sol_w2_ctrl <- sol_w2_ctrl/sum(sol_w2_ctrl)
  gamma_ate <- t(matrix(sol_w2_ctrl, n1,n0))
  w1_marg_w2_ate <- colSums(gamma_ate)#/sum(rowSums(gamma))
  
  # sbw_ate_weight <- full_weight
  # w2_ate_weight <- full_weight_w2
  # sbw_ate_weight[z==1] <- weights_ate
  # w2_ate_weight[z==1] <- w1_marg_w2_ate
  
  mass0_ate <- weights[weights>0]
  mass0_ate_w2 <- w0_marg_w2[w0_marg_w2>0]
  mass1_ate <- weights_ate[weights_ate>0]
  mass1_ate_w2 <- w1_marg_w2_ate[w1_marg_w2_ate>0]
  cost_rm_ate <- cost[weights>0, weights_ate>0]
  cost_rm_w2_ate <- cost[w0_marg_w2>0, w1_marg_w2_ate>0]
  
  w2_b4_sbw_ate <- transport::wasserstein(mass0_ate, mass1_ate, p = 2, tplan = NULL, costm = cost_rm_ate)
  w2_match_w2_ate <- transport::wasserstein(mass0_ate_w2, mass1_ate_w2, p = 2, tplan = NULL, costm = cost_rm_w2_ate)
  
  # pc_ate <- data.frame(id = 1:n, prcomp(x)$x[,1:2], z= factor(z), unmatch = unmatch_weight, 
  #                  sbw = full_weight, 
  #                  w2 = full_weight_w2,
  #                  sbw_ate = sbw_ate_weight,
  #                  w2_ate = w2_ate_weight)
  # pc_long_ate <- tidyr::gather(pc_ate, "condition", "weight", unmatch:w2_ate, factor_key = TRUE)
  # 
  # ggplot(pc_long_ate, aes(x = PC1, y = PC2, color = z, size = weight)) + 
  #   geom_point(alpha = 0.7) + facet_wrap(~condition) + 
  #   ggthemes::scale_color_tableau() + theme_bw()
  naive <- mean(y1) - mean(y0)
  sbw_drh_ate <- outcome_model_DRH(x,y,z, weights = list(z0 = weights, z1 = weights_ate), est = "ATE")
  sbw_h_ate <- outcome_model_H(x,y,z, weights = list(z0 = weights, z1 = weights_ate), est = "ATE")
  
  w2_drh_ate <- outcome_model_DRH(x,y,z, weights = list(z0 = w0_marg_w2, z1 = w1_marg_w2_ate), est = "ATE")
  w2_h_ate <- outcome_model_H(x,y,z, weights = list(z0 = w0_marg_w2, z1 = w1_marg_w2_ate), est = "ATE")
  
  # print(c("Naive" = naive,
  #         "SBW Hajek" = sbw_h_ate,
  #         "W2 Hajek" =w2_h_ate,
  #         "SBW DR Hajek" = sbw_drh_ate,
  #         "W2 DR Hajek" = w2_drh_ate))
  out <- list(ATT = list("Naive" = naive_att,
                         "SBW Hajek" = sbw_h,
                         "W2 Hajek" =w2_h,
                         "SBW DR Hajek" = sbw_drh,
                         "W2 DR Hajek" = w2_drh),
              ATE = list("Naive" = naive,
                         "SBW Hajek" = sbw_h_ate,
                         "W2 Hajek" =w2_h_ate,
                         "SBW DR Hajek" = sbw_drh_ate,
                         "W2 DR Hajek" = w2_drh_ate),
              W2 = list("pre-match" = w2_b4_w2, "ATT sbw" = w2_match_sbw, "ATT w2" = w2_match_w2,
                        "ATE sbw" = w2_b4_sbw_ate, "ATE w2" = w2_match_w2_ate)
  )
  return(out) 
}
stopCluster(cl)

ATT <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ATT)))
ATE <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ATE)))
W2 <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$W2)))

colMeans(ATT)
colMeans(ATE)
colMeans(W2)

colVar(ATT)#*(nsims-1)/nsims
colVar(ATE)#*(nsims-1)/nsims
colVar(W2)
