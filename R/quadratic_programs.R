quadprog.DataSim <- function(data, constraint,  estimand = c("ATT", "ATC", "ATE","cATE","feasible"), 
                             method = supported.methods(),
                             ...) {
  meth <- match.arg(method, supported.methods())
  est <- match.arg(estimand)
  
  qp <- if(meth == "SBW") {
    if(length(constraint) != data$get_p())  {
      K <- rep(constraint, data$get_p())[1:data$get_p()]
    } else {
      K <- constraint
    }
    if(est == "cATE") {
      list(qp_sbw(x=data$get_x(), z = data$get_z(), K = K, estimand = "ATC"),
           qp_sbw(x=data$get_x(), z = data$get_z(), K = K, estimand = "ATT"))
    } else {
      list(qp_sbw(x=data$get_x(), z = data$get_z(), K = K, estimand = est))
    } 
  } else if (meth == "Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    if(dots$dist == "RKHS" & is.null(dots$rkhs.args)) {
      if(is.null(dots$opt.method)) dots$opt.method <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
                         )
      dots$rkhs.args <- RKHS_param_opt(x=data$get_x(), 
                                       z=data$get_z(), 
                                       y = data$get_y(),
                                       metric = "mahalanobis",
                                       is.dose = dots$is.dose,
                                       opt.method = dots$opt.method,
                                       estimand = temp.est,
                                       ...)
    }
    if(est == "cATE") {
      list(qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
                   p=dots$p, estimand = "ATC", dist = dots$dist, cost = dots$cost,
                   rkhs.args = dots$rkhs.args),
           qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
                   p=dots$p, estimand = "ATT", dist = dots$dist, cost = dots$cost,
                   rkhs.args = dots$rkhs.args))
    } else if (est == "ATE") {
      qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
              p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost,
              rkhs.args = dots$rkhs.args)
    } else {
      list(qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
                   p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost,
                   rkhs.args = dots$rkhs.args))
    }
    
  } else if (meth == "Constrained Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    if(dots$dist == "RKHS" & is.null(dots$rkhs.args)) {
      if(is.null(dots$opt.method)) dots$opt.method <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
      )
      dots$rkhs.args <- RKHS_param_opt(x=data$get_x(), 
                                       z=data$get_z(), 
                                       y = data$get_y(),
                                       metric = "mahalanobis",
                                       is.dose = dots$is.dose,
                                       opt.method = dots$opt.method,
                                       estimand = temp.est,
                                       ...)
    }
    if (est == "cATE") {
      list(qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
                   p=dots$p, estimand = "ATC", dist = dots$dist, 
                   cost = dots$cost, 
                   rkhs.args = dots$rkhs.args),
           qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
                   p=dots$p, estimand = "ATT", 
                   dist = dots$dist, cost = dots$cost,
                   rkhs.args = dots$rkhs.args))
    } else if (est == "ATE") {
      qp_wass_const(x=data$get_x(), z=data$get_z(),K=constraint, 
              p=dots$p, estimand = est, 
              dist = dots$dist, cost = dots$cost,
              rkhs.args = dots$rkhs.args)
    } else {
      list(qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
                    p=dots$p, estimand = est, 
                    dist = dots$dist, cost = dots$cost,
                    rkhs.args = dots$rkhs.args))
    }
  } else if (meth == "RKHS" | meth == "RKHS.dose") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "mahalanobis"
    if(is.null(dots$theta)) dots$theta <- c(1,1)
    if(is.null(dots$gamma)) dots$gamma <- c(1,1)
    if(is.null(dots$lambda)) dots$lambda <- 0
    if(is.null(dots$sigma_2)) dots$sigma_2 <- 0
    
    if(est == "cATE") {
      list(qp_rkhs(x=data$get_x(), z=data$get_z(),
            p=dots$p, theta = dots$theta, gamma = dots$gamma,
            lambda = dots$lambda, sigma_2 = dots$sigma_2,
            dist = dots$dist, cost = dots$cost,
            is.dose = isTRUE(meth == "RKHS.dose"),
            estimand = "ATT"),
           qp_rkhs(x=data$get_x(), z=data$get_z(),
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$dist, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = "ATT"))
    } else {
      list(qp_rkhs(x=data$get_x(), z=data$get_z(),
              p=dots$p, theta = dots$theta, gamma = dots$gamma,
              lambda = dots$lambda, sigma_2 = dots$sigma_2,
              dist = dots$dist, cost = dots$cost,
              is.dose = isTRUE(meth == "RKHS.dose"),
              estimand = estimand))
    }
  }
  return(qp)
}

quadprog.data.frame <- function(data, constraint,  
                                estimand = c("ATT", "ATC", "ATE", "cATE", "feasible"), 
                                method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                                ...) {
  meth <- match.arg(method)
  est <- match.arg(estimand)
  
  df <- prep_data(data, ...)
  z <- as.numeric(df$z)
  x <- as.matrix(df$df[,!(colnames(df$df) == "y")])
  col_x <- ncol(x)
  
  qp <- if(meth == "SBW") {
    if(length(constraint) != col_x)  {
      K <- rep(constraint, col_x)[1:col_x]
    } else {
      K <- constraint
    }
    if(est == "cATE") {
      list(qp_sbw(x=x, z = z, K = K, estimand = "ATC"),
           qp_sbw(x=x, z = z, K = K, estimand = "ATT"))
    } else {
      list(qp_sbw(x=x, z = z, K = K, estimand = est))
    } 
  } else if (meth == "Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    if(dots$dist == "RKHS" & is.null(dots$rkhs.args)) {
      if(is.null(dots$opt.method)) dots$opt.method <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
      )
      dots$rkhs.args <- RKHS_param_opt(x=x, 
                                       z=z, 
                                       y = y,
                                       metric = "mahalanobis",
                                       is.dose = dots$is.dose,
                                       opt.method = dots$opt.method,
                                       estimand = temp.est,
                                       ...)
    }
    if(est == "cATE") {
      list(qp_wass(x=x, z=z, 
                   p=dots$p, estimand = "ATC", dist = dots$dist, cost = dots$cost),
           qp_wass(x=x, z=z,
                   p=dots$p, estimand = "ATT", dist = dots$dist, cost = dots$cost))
    } else if (est == "ATE") {
      qp_wass(x=x, z=z,
              p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost)
    } else {
      list(qp_wass(x=x, z=z,
                   p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost))
    }
    
  } else if (meth == "Constrained Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    if(dots$dist == "RKHS" & is.null(dots$rkhs.args)) {
      if(is.null(dots$opt.method)) dots$opt.method <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
      )
      dots$rkhs.args <- RKHS_param_opt(x=x, 
                                       z=z, 
                                       y = y,
                                       metric = "mahalanobis",
                                       is.dose = dots$is.dose,
                                       opt.method = dots$opt.method,
                                       estimand = temp.est,
                                       ...)
    }
    if (est == "cATE") {
      list(qp_wass_const(x=x, z=z,K=constraint, 
                   p=dots$p, estimand = "ATC", dist = dots$dist, cost = dots$cost),
           qp_wass_const(x=x, z=z,K=constraint, 
                   p=dots$p, estimand = "ATT", dist = dots$dist, cost = dots$cost))
    } else if (est == "ATE") {
      qp_wass_const(x=x, z=z,K=constraint, 
              p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost)
    } else {
      list(qp_wass_const(x=x, z=z, K=constraint, 
                         p=dots$p, estimand = est, dist = dots$dist, cost = dots$cost))
    }
  } else if (meth == "RKHS" | meth == "RKHS.dose") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "mahalanobis"
    if(is.null(dots$theta)) dots$theta <- c(1,1)
    if(is.null(dots$gamma)) dots$gamma <- c(1,1)
    if(is.null(dots$lambda)) dots$lambda <- 0
    if(is.null(dots$sigma_2)) dots$sigma_2 <- 0
    
    if(est == "cATE") {
      list(qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$dist, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = "ATC"),
           qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$dist, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = "ATT"))
    } else {
      list(qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$dist, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand - est))
    }
  }
  return(qp)
}

qp_sbw <- function(x, z, K, estimand = c("ATT", "ATC", 
                                         "ATE", "feasible")) {
  est <- match.arg(estimand)
  if (est == "ATC") {
    z <- 1 - z
  }
  
  stopifnot(length(K) == ncol(x))
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  x1_var <- colVar(x1)
  x0_var <- colVar(x0)
  pool_sd <- sqrt((x1_var + x0_var)/2)
  K_sd <- K * pool_sd
  K_upr <- K_sd
  K_lwr <- (-K_sd)
  
  if(est == "ATE") {
    n <- nrow(x)
    n1 <- nrow(x1)
    n0 <- nrow(x0)
    
    d <- ncol(x)
    Q0 <- Matrix::Diagonal(n, 1) 
    A1 <- Matrix::sparseMatrix(i = c(rep(1,n0),
                                     rep(2,n1)),
                               j = c(1:n0, n0+1:n1),
                               x = 1,
                               dims = c(2,n))
    
    x_constraint <- t(cbind(rbind(x0,matrix(0,nrow(x1),d)),rbind(matrix(0,nrow(x0),d),x1)))
    xm <- rep(colMeans(x), 2)
    K_upr <- rep(K_upr,2) + xm
    K_lwr <- rep(K_lwr,2) + xm
  } else if (est != "feasible") {
    n <- nrow(x0)
    d <- ncol(x0)
    x_constraint <- t(x0)
    x1m <- colMeans(x1)
    K_upr <- K_upr + x1m
    K_lwr <- K_lwr + x1m
    Q0 <- Matrix::Diagonal(n, 1) 
    A1 <- matrix(1,nrow = 1, ncol = n)
  } else if(est == "feasible") {
    n <- nrow(x)
    d <- ncol(x)
    x_diff <- x * matrix(ifelse(z == 0, 1, -1), nrow=n, ncol=d)
    x_constraint <- t(x_diff[order(z),])
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    Q0_c <- diag(1, n0, n0)#2 * (diag(1, n0, n0) - matrix(1/n0, n0, n0))
    Q0_t <- diag(1, n1, n1)#2 * (diag(1, n1, n1) - matrix(1/n1, n1, n1))
    Q0 <- Matrix::sparseMatrix(i = c(rep(1:n0, each = n0), rep(n0 + 1:n1, each = n1)),
                               j = c(rep.int(1:n0, n0), rep.int(n0 + 1:n1, n1)),
                               x = c(Q0_c, Q0_t),
                               dims = c(n,n), giveCsparse = FALSE
    )
    A1 <- Matrix::sparseMatrix(i = c(rep.int(1,n0), rep.int(2,n1)),
                               j = c(1:n0, n0 + 1:n1),
                               x = 1,
                               dims = c(2,n))
  }
  
  
  L0 <- c(w = rep(0,n))
  # var_bounds <- ROI::V_bound(li =1:n, ui = 1:n, lb = rep(0,n), ub = rep(1,n))
  # var_bounds <- NULL
  # op <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = L0, names = as.character(1:n)),
  #               bounds = var_bounds,
  #               maximum = FALSE)
  
  
  A2 <- x_constraint
  
  A <- rbind(A1, A2, A2)
  vals <- c(rep(1, nrow(A1)), K_lwr, K_upr)
  dir <- c(rep("E", nrow(A1)),rep("G",length(K_lwr)), rep("L",length(K_upr)))
  # dir <- c(ROI::eq(1), ROI::geq(d), ROI::leq(d))
  # LC <- ROI::L_constraint(A, dir, vals)
  
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
  
  # ROI::constraints(op) <- LC
  
  program <- list(obj = list(Q = Q0, L = L0),
                   LC = list(A = A, dir = dir,
                             vals = vals))
  return(program)
}

qp_wass_const <- function(x, z, K, p = 2, estimand = c("ATC", "ATT", 
                                                       "ATE", "feasible"),
                          dist=c("Lp", "mahalanobis","RKHS"), cost = NULL,
                          rkhs.args = NULL) {
  estimand <- match.arg(estimand)
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  # x1_var <- cov(x1)
  # x0_var <- cov(x0)
  
  dist <- match.arg(dist)
  # dist.fun <- switch(dist, 
  #                    "Lp" = causalOT::cost_calc_lp,
  #                    "mahalanobis" = causalOT::cost_mahalanobis,
  #                    "RKHS" = causalOT::cost_RKHS)
  # covar <- (x1_var + x0_var)/2 # add later
  if(is.null(cost)) { 
    cost <- cost_fun(x0, x1, ground_p = p, metric = dist,
                     rkhs.args = rkhs.args, estimand = estimand)
  } else {
    if(estimand != "ATE") {
      stopifnot(all(dim(cost) %in% c(nrow(x1),nrow(x0))))
    }
  }
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
  
  if(estimand != "ATE") {
    cost_vec <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                   j = 1:(n0*n1),
                                   x= c(cost)^p,
                                   dims = c(1, n0*n1))
  }
  sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                    j = 1:(n0*n1),
                                    x = rep(1, n1 * n0),
                                    dims = c(1, n0*n1))
  q_s <- marg_const <- marg_const_mat <- NULL
  marg_const_n <- 0
  
  # Qmat
  if(estimand == "ATT") {
    # q_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1,n0*n1,n0)
    q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
                                j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                x = 1,
                                dims = c(n0,n1*n0), giveCsparse = FALSE)
    # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
    marg_mass <- 1/n0
    marg_n <- n0
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n1, each = n0),
                                           j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
                                           x = rep(1,n0*n1),
                                           dims = c(n1,n0*n1), giveCsparse = FALSE)
    # marg_const_mat <- matrix(0, n1,n1*n0)
    # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
    marg_const <- rep( 1/n1, n1)
    marg_const_n <- n1
    # q_m <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # q_m <- matrix(-1/n0, n0,n0)
  } 
  else if (estimand == "ATC") {
    n0_idx <- 1:n0
    q_s <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
                                j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
                                x = 1,
                                dims = c(n1,n1*n0), giveCsparse = FALSE)
    
    marg_mass <- 1/n1
    marg_n <- n1
    n1_idx <- seq(1,n0*n1,n0)
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n0, each = n1),
                                           j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                           x= rep.int(1,n0*n1),
                                           dims = c(n0, n0*n1), giveCsparse = FALSE)
    
    marg_const <- rep.int( 1/n0, n0 )
    marg_const_n <- n0
    # q_m <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
  } 
  else if (estimand == "feasible") {
    n1_idx <- seq(1,n0*n1,n0)
    n0_idx <- 1:n0
    q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
                                j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                x = 1,
                                dims = c(n0,n1*n0), giveCsparse = FALSE)
    
    q_t <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
                                j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
                                x = 1,
                                dims = c(n1,n1*n0), giveCsparse = FALSE)
    q_s <- rbind(q_c,q_t)
    # Q0_c <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # Q0_t <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
    # Q0_c <- (matrix(-1/n0, n0,n0))
    # Q0_t <- (matrix(-1/n1, n1,n1))
    # q_m <- Matrix::sparseMatrix(i = c(rep(1:n0, each = n0), rep(n0 + 1:n1, each = n1)),
    #                             j = c(rep.int(1:n0, n0), rep.int(n0 + 1:n1, n1)),
    #                             x = c(Q0_c, Q0_t),
    #                             dims = c(n,n)
    # )
    rm(q_t, q_c)
  } else if (estimand == "ATE") {
    if(is.list(cost)) {
      cost_n0 <- cost[[1]]
      cost_n1 <- cost[[2]]
    } else {
      if(dist != "RKHS") {
        cost_n0 <- cost_fun(x0, x, ground_p = p, metric = dist)
        cost_n1 <- cost_fun(x1, x, ground_p = p, metric = dist)
        
      } else {
        cost    <- cost_fun(x0, x1, ground_p = p, metric = dist,
                            rkhs.args = rkhs.args, estimand = estimand)
        cost_n0 <- cost[[1]]
        cost_n1 <- cost[[2]]
      }
    }
    
    if(length(K) >=2) {
      K <- K[1:2]
    } else {
      K <- rep(K,2)
    }
    
    cost_vec_n1 <- Matrix::sparseMatrix(i = rep(1, n*n1),
                                        j = 1:(n*n1),
                                        x= c(cost_n1)^p,
                                        dims = c(1, n*n1))
    
    cost_vec_n0 <- Matrix::sparseMatrix(i = rep(1, n*n0),
                                        j = 1:(n*n0),
                                        x= c(cost_n0)^p,
                                        dims = c(1, n0*n))
    
    sum_const_n1 <- Matrix::sparseMatrix(i = rep(1, n*n1),
                                         j = 1:(n*n1),
                                         x = rep(1, n1 * n),
                                         dims = c(1, n*n1))
    sum_const_n0 <- Matrix::sparseMatrix(i = rep(1, n0*n),
                                         j = 1:(n0*n),
                                         x = rep(1, n * n0),
                                         dims = c(1, n0*n))
    
    marg_const <- rep(1/n,n)
    marg_n <- n
    
    n1_idx <- 1:n1
    n0_idx <- 1:n0
    n_idx <- 1:n
    
    marg_const_mat_n1 <- Matrix::sparseMatrix(i = rep(1:n, each = n1),
                                              j = c(sapply(0:(n-1), function(i) i * n1 + n1_idx)),
                                              x = rep(1,n*n1),
                                              dims = c(n,n*n1))
    
    marg_const_mat_n0 <- Matrix::sparseMatrix(i = rep(1:n, each = n0),
                                              j = c(sapply(0:(n-1), function(i) i* n0 + n0_idx)),
                                              x= rep.int(1,n0*n),
                                              dims = c(n, n0*n))
    
    q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n) , 
                                j = c(sapply(0:(n-1), function(i) i * n0 + n0_idx)),
                                x = 1,
                                dims = c(n0,n*n0), giveCsparse = FALSE)
    
    q_t <- Matrix::sparseMatrix(i = rep(1:n1, each = n), 
                                j = c(sapply(0:(n-1), function(i) i * n1 + n1_idx)),
                                x = 1,
                                dims = c(n1,n1*n), giveCsparse = FALSE)
    op <- vector("list",2)
    op[[1]] <- list(obj = list(Q =  Matrix::crossprod(q_c),
                               L = c(rep(0, n0*n))),
                    LC = list(A = rbind(cost_vec_n0, sum_const_n0, marg_const_mat_n0),
                              vals = c(K[1]^p, 1, marg_const),
                              dir = c("L", rep("E", 1 + marg_n))))
    rm(q_c)
    
    op[[2]] <- list(obj = list(Q =  Matrix::crossprod(q_t),
                               L = c(rep(0, n1*n))),
                    LC = list(A = rbind(cost_vec_n1, sum_const_n1, marg_const_mat_n1),
                              vals = c(K[2]^p, 1, marg_const),
                              dir = c("L", rep("E", 1 + marg_n))))
    rm(q_t)
    return(op)
  }
  # mean_mat <- as.simple_triplet_matrix(matrix(marg_mass, marg_n, n1*n0))
  # q_s <- Matrix::sparseMatrix(i = q_mat_list$i,
  #                             j = q_mat_list$j,
  #                             x = q_mat_list$val)  #Matrix::Matrix(q_mat, sparse = TRUE)
  # q_self  <- Matrix::crossprod(q_s)
  # q_mult  <- Matrix::crossprod(q_m_chol, q_s)
  # q_cross  <- Matrix::crossprod(q_s)
  Q0 <-  Matrix::crossprod(q_s)
  rm(q_s)
  
  L0 <- c(rep(0, n0*n1))#simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  obj <- list(Q = Q0,
              L = L0)
  
  LC <- list()
  LC$A <- rbind(cost_vec, sum_const, marg_const_mat)
  LC$vals <- c(K^p, 1, marg_const)
  LC$dir <- c("L", rep("E", 1 + marg_const_n))
  
  op <- list(obj=obj, LC=LC)
  return(op)
}


qp_wass <- function(x, z, K = NULL, p = 2, estimand = c("ATC", "ATT", "ATE",
                                            "feasible"),
                    dist=c("Lp", "mahalanobis"), cost = NULL,
                    rkhs.args = NULL) {
  sqrt_mat <- function(X) {
    p <- ncol(X)
    decomp <- eigen(X)
    return(decomp$vectors %*% diag(sqrt(decomp$values), p, p) %*% t(decomp$vectors))
  }
  estimand <- match.arg(estimand)
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  d_plus <- d+1
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  # x1_var <- cov(x1)
  # x0_var <- cov(x0)
  
  dist <- match.arg(dist)
  # dist.fun <- switch(dist, 
  #                    "Lp" = causalOT::cost_calc_lp,
  #                    "mahalanobis" = causalOT::cost_mahalanobis)
  # covar <- (x1_var + x0_var)/2 # add later
  if(is.null(cost)) { 
    if(estimand %in% c("ATT","ATC")) {
      cost <- lapply(1:d, function(i) cost_fun(x0[,i,drop=FALSE], x1[,i,drop=FALSE], 
                                               ground_p = p, metric = dist,
                                               rkhs.args = rkhs.args, estimand = estimand))
      cost[d+1] <- list(cost_fun(x0, x1, 
                              ground_p = p, metric = dist,
                              rkhs.args = rkhs.args, estimand = estimand))
    } else if (estimand == "ATE") {
      if(dist != "RKHS") {
        cost <- list(z0 = lapply(1:d, function(i) cost_fun(x0[,i,drop=FALSE], x[,i,drop=FALSE], 
                                                           ground_p = p, metric = dist,
                                                           rkhs.args = rkhs.args, estimand = "ATE")),
                     z1 = lapply(1:d, function(i) cost_fun(x1[,i,drop=FALSE], x[,i,drop=FALSE], 
                                                           ground_p = p, metric = dist,
                                                           rkhs.args = rkhs.args, estimand = "ATE")))
        cost$z0[d+1] <- list(cost_fun(x0, x, 
                                 ground_p = p, metric = dist,
                                 rkhs.args = rkhs.args, estimand = "ATE"))
        cost$z1[d+1] <- list(cost_fun(x1, x, 
                                 ground_p = p, metric = dist,
                                 rkhs.args = rkhs.args, estimand = "ATE"))
      } else {
        cost_temp    <- lapply(1:d, function(i) cost_fun(x0[,i,drop=FALSE], x1[,i,drop=FALSE], 
                                                         ground_p = p, metric = dist,
                                                         rkhs.args = rkhs.args, estimand = "ATE"))
        cost <- list(z0 = lapply(cost_temp, function(cc) cc[[1]]),
                     z1 = lapply(cost_temp, function(cc) cc[[2]]))
        cost$z0[d+1] <- list(cost_fun(x0, x, 
                                 ground_p = p, metric = dist,
                                 rkhs.args = rkhs.args, estimand = "ATE"))
        cost$z1[d+1] <- list(cost_fun(x1, x, 
                                 ground_p = p, metric = dist,
                                 rkhs.args = rkhs.args, estimand = "ATE"))
      }
      
    }
    
  } else {
    stopifnot(all(sapply(cost, dim) %in% c(nrow(x1),nrow(x0))))
  }
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
  
  if(estimand != "ATE") {
    weight.dim <- n0 * n1
    cost_vec <- Matrix::sparseMatrix(i = rep(1:d_plus, each = n0*n1),
                                     j = rep(1:(n0*n1), d_plus),
                                     x= c(sapply(cost, c)^p),
                                     dims = c(d_plus, n0*n1))
  } 
  
  sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                    j = 1:(n0*n1),
                                    x = rep(1, n1 * n0),
                                    dims = c(1, n0*n1))
  q_s <- marg_const_mat <- marg_const <- NULL
  marg_const_n <- 0
  
  # Qmat
  if(estimand == "ATT") {
    # q_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1, weight.dim, n0)
    # q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
    #                             j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
    #                             x = 1,
    #                             dims = c(n0,n1*n0))
    # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
    marg_mass <- 1/n0
    marg_n <- n0
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n1, each = n0),
                                           j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
                                           x = rep(1, weight.dim),
                                           dims = c(n1, weight.dim))
    # marg_const_mat <- matrix(0, n1,n1*n0)
    # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
    marg_const <- rep(1/n1, n1)
    marg_const_n <- n1
    # q_m <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # q_m <- matrix(-1/n0, n0,n0)
  } else if (estimand == "ATC") {
    n0_idx <- 1:n0
    # q_s <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
    #                             j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
    #                             x = 1,
    #                             dims = c(n1,n1*n0))
    
    marg_mass <- 1/n1
    marg_n <- n1
    n1_idx <- seq(1, weight.dim, n0)
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n0, each = n1),
                                           j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                           x= rep.int(1, weight.dim),
                                           dims = c(n0, weight.dim))
    
    marg_const <- rep.int( 1/n0, n0 )
    marg_const_n <- n0
    # q_m <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
  } else if (estimand == "feasible") {
    n1_idx <- seq(1,n0*n1,n0)
    n0_idx <- 1:n0
    # q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) ,
    #                             j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
    #                             x = 1,
    #                             dims = c(n0,n1*n0))
    # 
    # q_t <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) ,
    #                             j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
    #                             x = 1,
    #                             dims = c(n1,n1*n0))
    # q_s <- rbind(q_c,q_t)
    # q_s   <- Matrix::sparseMatrix(i = 1:(n1*n0),
    #                               j = 1:(n1*n0),
    #                               x = 1,
    #                               dims = c(n1*n0,n1*n0))
    # Q0_c <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # Q0_t <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
    # Q0_c <- (matrix(-1/n0, n0,n0))
    # Q0_t <- (matrix(-1/n1, n1,n1))
    # q_m <- Matrix::sparseMatrix(i = c(rep(1:n0, each = n0), rep(n0 + 1:n1, each = n1)),
    #                             j = c(rep.int(1:n0, n0), rep.int(n0 + 1:n1, n1)),
    #                             x = c(Q0_c, Q0_t),
    #                             dims = c(n,n)
    # )
    # rm(q_t, q_c)
  } else if (estimand == "ATE") {
    
    weight.dim1 <- n * n1
    weight.dim0 <- n * n0
    
    if(is.list(cost)) {
      cost_n0 <- cost[[1]]
      cost_n1 <- cost[[2]]
    } else {
      stop("Cost must be a list")
    } 
    
    cost_vec_n1 <- Matrix::sparseMatrix(i = rep(1:d_plus, each = weight.dim1),
                                        j = rep(1:weight.dim1, d_plus ),
                                        x= c(sapply(cost_n1,c))^p,
                                        dims = c( d_plus, weight.dim1))
    
    cost_vec_n0 <- Matrix::sparseMatrix(i = rep(1:d_plus, each = weight.dim0),
                                        j = rep(1:weight.dim0,d_plus),
                                        x= c(sapply(cost_n0,c))^p,
                                        dims = c(d_plus , weight.dim0))
    
    sum_const_n1 <- Matrix::sparseMatrix(i = rep(1, weight.dim1),
                                         j = 1:weight.dim1,
                                         x = rep(1, weight.dim1),
                                         dims = c(1, weight.dim1))
    sum_const_n0 <- Matrix::sparseMatrix(i = rep(1, weight.dim0),
                                         j = 1:weight.dim0,
                                         x = rep(1, weight.dim0),
                                         dims = c(1, weight.dim0))
    
    marg_const <- rep(1/n,n)
    marg_n <- n
    
    n1_idx <- 1:n1
    n0_idx <- 1:n0
    n_idx <- 1:n 
    
    marg_const_mat_n1 <- Matrix::sparseMatrix(i = rep(1:n, each = n1),
                                              j = c(sapply(0:(n-1), function(i) i * n1 + n1_idx)),
                                              x = rep(1,n*n1),
                                              dims = c(n,n*n1))
    
    marg_const_mat_n0 <- Matrix::sparseMatrix(i = rep(1:n, each = n0),
                                              j = c(sapply(0:(n-1), function(i) i* n0 + n0_idx)),
                                              x= rep.int(1,n0*n),
                                              dims = c(n, n0*n))
    
    if(is.null(K)) {
      sigma1 <- cov(x1)
      sigma0 <- cov(x0)
      sigma  <- cov(x)
      v1 <- diag(sigma1)
      v0 <- diag(sigma0)
      v <- diag(sigma)
      shared_var0 <- 0.5 * v0 + 0.5 * v
      shared_var1 <- 0.5 * v1 + 0.5 * v
      
      max0 <- sapply(1:d, function(i) transport::wasserstein1d(a = x0[,i], b = x[,i],p = p))
      max1 <- sapply(1:d, function(i) transport::wasserstein1d(a = x1[,i], b = x[,i],p = p))
      max0 <- c(max0, transport::wasserstein(a = rep(1/n0,n0), b = rep(1/n,n), 
                                             p = p, costm = cost$z0[[d_plus]]))^p
      max1 <- c(max1, transport::wasserstein(a = rep(1/n1,n1), b = rep(1/n,n),
                                             p = p, costm = cost$z1[[d_plus]]))^p
      
      min0 <- sapply(1:d_plus, function(i) mean(apply(cost$z0[[i]]^p, 2, min)))
      min1 <- sapply(1:d_plus, function(i) mean(apply(cost$z1[[i]]^p, 2, min)))
      
      # w1 <-  tabulate(apply(cost$z1[[d_plus]]^p, 2, which.min), nbins = n1)/n
      # w0 <-  tabulate(apply(cost$z0[[d_plus]]^p, 2, which.min), nbins = n0)/n
      
      # min1 <- sapply(1:6, function(i) transport::wasserstein(a = w1, b = rep(1/n,n),
      #                                p = p, costm = cost$z1[[i]]))
      # min0 <- sapply(1:6, function(i) transport::wasserstein(a = w0, b = rep(1/n,n),
      #                                                        p = p, costm = cost$z0[[i]]))
      
      K_vals_0 <- c(
             sqrt(0.2^2 * shared_var0  + v + v0 - 2 * sqrt(sqrt(v) *v0 * sqrt(v))),
             sqrt(sum(0.2^2 * shared_var0) + sum(v) + sum(v0) - 
               2 * sum(diag(sqrt_mat(sqrt_mat(sigma) %*% sqrt_mat(sigma0)))))
      )^p
      K_vals_1 <- c(
        sqrt(abs(0.2^2 * shared_var1   + v + v0 - 2 * sqrt(sqrt(v) *v1 * sqrt(v)))),
        sqrt(abs(sum(0.2^2 * shared_var1 ) + sum(v) + sum(v1) -
          2 * sum(diag(sqrt_mat(sqrt_mat(sigma) %*% sigma1 %*% sqrt_mat(sigma))))))
      )^p
      increase.factor.0 <- ifelse(max0/min0 < 2, max0/min0, 2)
      increase.factor.1 <- ifelse(max1/min1 < 2, max1/min1, 2)
      K_vals_0 <- ifelse(K_vals_0 < min0, min0 * increase.factor.0, K_vals_0)
      K_vals_1 <- ifelse(K_vals_1 < min1, min1 * increase.factor.1, K_vals_1)
    } else if(length(K) == 2 ) {
      if(is.numeric(K)) {
        K <- K^p
        K_vals_0 <- rep(K[1], d_plus)
        K_vals_1 <- rep(K[2], d_plus)
        
        K_vals_0[d_plus] <- K[1] * d
        K_vals_1[d_plus] <- K[2] * d
      } else if (is.list(K)) {
        K0 <- K[[1]]
        K1 <- K[[2]]
        
        if(length(K1) != d_plus ){
          K1 <- K1^p
          K_vals_1 <- rep(K1, d_plus)[1:d_plus]
        } else if (length(K1) == d_plus) {
          K_vals_1 <- K1
        }
        if(length(K0) == 1 ){
          K0 <- K1^p
          K_vals_2 <- rep(K0, d_plus)[1:d_plus]
        } else if (length(K0) == d_plus) {
          K_vals_0 <- K0
        }
      }
      
      
    } else if(length(K) == 1 & is.numeric(K)) {
      K <- K^p
      K_vals_0 <- K_vals_1 <-  rep(K, d_plus)[1:d_plus]
      
      K_vals_0[d_plus] <- K * d
      K_vals_1[d_plus] <- K * d
    } else {
      K <- K^p
      K_vals_0 <- K_vals_1 <- K
    }
    
    # q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n) , 
    #                             j = c(sapply(0:(n-1), function(i) i * n0 + n0_idx)),
    #                             x = 1,
    #                             dims = c(n0,n*n0))
    # 
    # q_t <- Matrix::sparseMatrix(i = rep.int(1:n1, each = n) , 
    #                             j = c(sapply(0:(n1-1), function(i) i * n1 + n1_idx)),
    #                             x = 1,
    #                             dims = c(n1,n1*n))
    op <- vector("list",2)
    op[[1]] <- list(obj = list(Q =  Matrix::Diagonal(weight.dim0, 1),
                               L = rep(0, weight.dim0)),
                    LC = list(A = rbind(sum_const_n0, marg_const_mat_n0,
                                        cost_vec_n0),
                              vals = c(1, marg_const,
                                       K_vals_0),
                              dir = c(rep("E", 1 + marg_n),
                                      rep("L", d_plus)
                                      )))
    # rm(q_c)
    
    op[[2]] <- list(obj = list(Q =  Matrix::Diagonal(weight.dim1, 1),
                               L = rep(0, weight.dim1)),
                    LC = list(A = rbind(sum_const_n1, marg_const_mat_n1,
                                        cost_vec_n1),
                              vals = c(1, marg_const,
                                       K_vals_1),
                              dir = c(rep("E", 1 + marg_n),
                                      rep("L", d_plus))))
    # rm(q_t)
    return(op)
  }
  # mean_mat <- as.simple_triplet_matrix(matrix(marg_mass, marg_n, n1*n0))
  # q_s <- Matrix::sparseMatrix(i = q_mat_list$i,
  #                             j = q_mat_list$j,
  #                             x = q_mat_list$val)  #Matrix::Matrix(q_mat, sparse = TRUE)
  # q_self  <- Matrix::crossprod(q_s)
  # q_mult  <- Matrix::crossprod(q_m_chol, q_s)
  # q_cross  <- Matrix::crossprod(q_s)
  # if(!Matrix::isDiagonal(q_s)) {
  #   Q0 <-  2*Matrix::crossprod(q_s)
  # } else {
  #   Q0 <-  2*q_s
  # }
  # rm(q_s)
  if(is.null(K)) {
    sigma1 <- cov(x1)
    sigma0 <- cov(x0)
    v1 <- diag(sigma1)
    v0 <- diag(sigma0)
    shared_var <- 0.5 * v1 + 0.5 * v0
    
    K <- c(sqrt(0.2^2 * shared_var * 2  + v1 + v0 - 2 * sqrt(sqrt(v1) *v0 * sqrt(v1))),
           sum(0.2^2 * shared_var * 2) + sum(v1) + sum(v0) - 
             2 * sum(diag(sqrt_mat(sqrt_mat(sigma1) %*% sigma0 %*% sqrt_mat(sigma1))))
           )
  }
  if(length(K) == 1) {
    K_vals <- rep(K, d_plus)
  } else if (length(K) == d_plus) {
    K_vals <- K^p
  } else {
    stop("K must be length 1 or ncol covariates +1 for ATC/ATT")
  }
  
  # K_vals <- rep(K, d_plus)
  L0 <- rep(0, n0*n1) #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  obj <- list(Q = Matrix::Diagonal(weight.dim, 1),#Q0,
              L = L0)
  
  LC <- list()
  LC$A <- rbind(sum_const, marg_const_mat, cost_vec)
  LC$vals <- c(1, marg_const, K_vals)
  LC$dir <- c(rep("E", 1 + marg_const_n),
              rep("L", d_plus))
  
  op <- list(obj=obj, LC=LC)
  return(op)
}

qp_rkhs <- function(x, z, p = 1, estimand = c("ATC", "ATT", "ATE"),
                    theta = c(1,1),
                    gamma = c(1,1), lambda = 0, 
                    sigma_2 = c(1,1),
                    dist = c("mahalanobis","Lp"), cost = NULL, 
                    is.dose = FALSE, ...) {
  # est <- match.arg(estimand)
  # if (est == "ATC") {
  #   z <- 1 - z
  # }
  n <- nrow(x)
  
  dist <- match.arg(dist)
  # cost.fun <- switch(isTRUE(is.dose),
  #                    "TRUE" = "kernel_calculation",
  #                    "FALSE" = 
  
  if(is.null(cost)) { 
    kernels <- kernel_calculation(x, z, p = p, theta = theta, 
                               gamma = gamma, sigma_2 = sigma_2, 
                               metric = dist, is.dose = is.dose,
                               estimand = estimand)
    cost <- list(kernels$cov_kernel, kernels$mean_kernel)
  } else {
    stopifnot(all(dim(cost) %in% n))
  }
  
  cost[[1]] <- Matrix::Matrix(data = cost[[1]], sparse = TRUE)
  
  if(is.dose) {
    Q0 <- (1/(n^2)) * (cost[[1]] + Matrix::Diagonal(n, lambda/(n^2)))
    
    L0 <- c(w = -2/(n^2) * c(cost[[2]]))
    
    A <- t(rep(1.0/n, n))
    vals <- as.double(n)
    dir <- "E"
    
  } else {
    Q0 <- cost[[1]]
    
    L0 <- c(w = -2 * c(cost[[2]]))
    
    A <- rbind(t(1/n * z),
               t(1/n * (1-z)  ))
    vals <- as.double(c(n,n))
    dir <- c("E", "E")
  }
  
  Q0 <- as(Q0, "dsTMatrix")
  
  
  quick_op <- list(obj = list(Q = Q0, L = L0),
                   LC = list(A = A, dir = dir,
                             vals = vals))
  return(quick_op)
}

check_wass_const <- function(opt_problem) {
  cost <- c(opt_problem$LC$A[1,])
  const <- c(opt_problem$LC$vals[1])
  n <- length(cost)
  mass <- rep.int(1/n,n)
  wass <- sum(cost * mass)
  
  output <- if(wass < const) {
    list(res = mass, skip_cplex = TRUE)
  } else {
    list(res = NULL, skip_cplex = FALSE)
  }
  return(output)
}


# setOldClass("DataSim")
setGeneric("quadprog", function(data, ...) UseMethod("quadprog"))
setMethod("quadprog", "DataSim", quadprog.DataSim)
setMethod("quadprog", "data.frame", quadprog.data.frame)

