quadprog.DataSim <- function(data, constraint,  estimate = c("ATT", "ATC","feasible"), 
                             method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                             ...) {
  meth <- match.arg(method)
  est <- match.arg(estimate)
  qp <- if(meth == "SBW") {
    if(length(constraint) != data$get_p())  {
      K <- rep(constraint, data$get_p())[1:data$get_p()]
    } else {
      K <- constraint
    }
    qp_sbw(x=data$get_x(), z = data$get_z(), K = K, estimate = est)
  } else if (meth == "Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    qp_wass(x=data$get_x(), z=data$get_z(),
            p=dots$p, target = est, dist = dots$dist, cost = dots$cost)
  } else if (meth == "Constrained Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
            p=dots$p, target = est, dist = dots$dist, cost = dots$cost)
  }
  return(qp)
}

quadprog.data.frame <- function(data, constraint,  estimate = c("ATT", "ATC","feasible"), 
                             method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                             ...) {
  meth <- match.arg(method)
  est <- match.arg(estimate)
  
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
    qp_sbw(x=x, z = z, K = K, estimate = est)
  } else if (meth == "Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    qp_wass(x=x, z=z,
            p=dots$p, target = est, dist = dots$dist, cost = dots$cost)
  } else if (meth == "Constrained Wasserstein") {
    dots <- list(...)
    if(is.null(dots$p)) dots$p <- 2
    if(is.null(dots$dist)) dots$dist <- "Lp"
    qp_wass_const(x=x, z=z, K=constraint, 
                  p=dots$p, target = est, dist = dots$dist, cost = dots$cost)
  }
  return(qp)
}
qp_sbw <- function(x, z, K, estimate = c("ATT", "ATC", "feasible")) {
  est <- match.arg(estimate)
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
  
  # pool_sd <- 1
  if (est != "feasible") {
    n <- nrow(x0)
    d <- ncol(x0)
    x_constraint <- t(x0)
    x1m <- colMeans(x1)
    K_upr <- K_upr + x1m
    K_lwr <- K_lwr + x1m
    Q0 <- diag(2, n, n) # times 2 since the Q part is 1/2 a^\top (x^\top x) a in ROI!!!
    A1 <- matrix(1,nrow = 1, ncol = n)
  } else {
    n <- nrow(x)
    d <- ncol(x)
    x_diff <- x * matrix(ifelse(z == 0, 1, -1), nrow=n, ncol=d)
    x_constraint <- t(x_diff[order(z),])
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    Q0_c <- diag(2, n0, n0)#2 * (diag(1, n0, n0) - matrix(1/n0, n0, n0))
    Q0_t <- diag(2, n1, n1)#2 * (diag(1, n1, n1) - matrix(1/n1, n1, n1))
    Q0 <- Matrix::sparseMatrix(i = c(rep(1:n0, each = n0), rep(n0 + 1:n1, each = n1)),
                               j = c(rep.int(1:n0, n0), rep.int(n0 + 1:n1, n1)),
                               x = c(Q0_c, Q0_t),
                               dims = c(n,n)
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
  dir <- c(rep("E", nrow(A1)),rep("G",d), rep("L",d))
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
  
  quick_op <- list(obj = list(Q = Q0, L = L0),
                   LC = list(A = A, dir = dir,
                             vals = vals))
  return(quick_op)
}

qp_wass_const <- function(x, z, K, p = 2, target = c("ATC", "ATT", 
                                               "feasible"),
                    dist=c("Lp", "mahalanobis"), cost = NULL) {
  target <- match.arg(target)
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  # x1_var <- cov(x1)
  # x0_var <- cov(x0)
  
  dist <- match.arg(dist)
  dist.fun <- switch(dist, 
                     "Lp" = causalOT::cost_calc_lp,
                     "mahalanobis" = causalOT::cost_mahalanobis)
  # covar <- (x1_var + x0_var)/2 # add later
  if(is.null(cost)) { 
    cost <- dist.fun(x0, x1, ground_p = p)
  } else {
    stopifnot(all(dim(cost) %in% c(nrow(x1),nrow(x0))))
  }
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
  
  cost_vec <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                   j = 1:(n0*n1),
                                   x= c(cost)^p,
                                   dims = c(1, n0*n1))
  sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                    j = 1:(n0*n1),
                                    x = rep(1, n1 * n0),
                                    dims = c(1, n0*n1))
  q_s <- marg_const <- marg_const_mat <- NULL
  marg_const_n <- 0
  
  # Qmat
  if(target == "ATT") {
    # q_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1,n0*n1,n0)
    q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
                                j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                x = 1,
                                dims = c(n0,n1*n0))
    # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
    marg_mass <- 1/n0
    marg_n <- n0
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n1, each = n0),
                                           j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
                                           x = rep(1,n0*n1),
                                           dims = c(n1,n0*n1))
    # marg_const_mat <- matrix(0, n1,n1*n0)
    # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
    marg_const <- rep( 1/n1, n1)
    marg_const_n <- n1
    # q_m <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # q_m <- matrix(-1/n0, n0,n0)
  } else if (target == "ATC") {
    n0_idx <- 1:n0
    q_s <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
                                j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
                                x = 1,
                                dims = c(n1,n1*n0))
    
    marg_mass <- 1/n1
    marg_n <- n1
    n1_idx <- seq(1,n0*n1,n0)
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n0, each = n1),
                                           j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                           x= rep.int(1,n0*n1),
                                           dims = c(n0, n0*n1))
    
    marg_const <- rep.int( 1/n0, n0 )
    marg_const_n <- n0
    # q_m <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
  } else if (target == "feasible") {
    n1_idx <- seq(1,n0*n1,n0)
    n0_idx <- 1:n0
    q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
                                j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                x = 1,
                                dims = c(n0,n1*n0))
    
    q_t <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
                                j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
                                x = 1,
                                dims = c(n1,n1*n0))
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
  }
  # mean_mat <- as.simple_triplet_matrix(matrix(marg_mass, marg_n, n1*n0))
  # q_s <- Matrix::sparseMatrix(i = q_mat_list$i,
  #                             j = q_mat_list$j,
  #                             x = q_mat_list$val)  #Matrix::Matrix(q_mat, sparse = TRUE)
  # q_self  <- Matrix::crossprod(q_s)
  # q_mult  <- Matrix::crossprod(q_m_chol, q_s)
  # q_cross  <- Matrix::crossprod(q_s)
  Q0 <-  2*Matrix::crossprod(q_s)
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

# qp_wass <- function(x, z, p = 2, target = c("ATC", "ATT", 
#                                             "feasible"),
#                     dist=c("Lp", "mahalanobis"), cost = NULL) {
#   target <- match.arg(target)
#   x1 <- x[z==1,,drop=FALSE]
#   x0 <- x[z==0,,drop=FALSE]
#   
#   n <- nrow(x)
#   d <- ncol(x)
#   
#   n1 <- nrow(x1)
#   n0 <- nrow(x0)
#   # x1_var <- cov(x1)
#   # x0_var <- cov(x0)
#   
#   dist <- match.arg(dist)
#   dist.fun <- switch(dist, 
#                      "Lp" = causalOT::cost_calc_lp,
#                      "mahalanobis" = causalOT::cost_mahalanobis)
#   # covar <- (x1_var + x0_var)/2 # add later
#   if(is.null(cost)) { 
#     cost <- dist.fun(x0, x1, ground_p = p)
#   } else {
#     stopifnot(all(dim(cost) %in% c(nrow(x1),nrow(x0))))
#   }
#   # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
#   
#   cost_vec <- c(cost)^p
#   sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
#                                     j = 1:(n0*n1),
#                                     x = rep(1, n1 * n0),
#                                     dims = c(1, n0*n1))
#   q_s <- marg_const_mat <- marg_const <- NULL
#   marg_const_n <- 0
#   
#   # Qmat
#   if(target == "ATT") {
#     # q_mat <- matrix(0, n0,n1*n0)
#     n1_idx <- seq(1,n0*n1,n0)
#     q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
#                                 j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
#                                 x = 1,
#                                 dims = c(n0,n1*n0))
#     # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
#     marg_mass <- 1/n0
#     marg_n <- n0
#     marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n1, each = n0),
#                                            j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
#                                            x = rep(1,n0*n1),
#                                            dims = c(n1,n0*n1))
#     # marg_const_mat <- matrix(0, n1,n1*n0)
#     # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
#     marg_const <- rep( 1/n1, n1)
#     marg_const_n <- n1
#     # q_m <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
#     # q_m <- matrix(-1/n0, n0,n0)
#   } else if (target == "ATC") {
#     n0_idx <- 1:n0
#     q_s <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
#                                 j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
#                                 x = 1,
#                                 dims = c(n1,n1*n0))
#     
#     marg_mass <- 1/n1
#     marg_n <- n1
#     n1_idx <- seq(1,n0*n1,n0)
#     marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n0, each = n1),
#                                            j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
#                                            x= rep.int(1,n0*n1),
#                                            dims = c(n0, n0*n1))
#     
#     marg_const <- rep.int( 1/n0, n0 )
#     marg_const_n <- n0
#     # q_m <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
#   } else if (target == "feasible") {
#     n1_idx <- seq(1,n0*n1,n0)
#     n0_idx <- 1:n0
#     q_c <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) ,
#                                 j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
#                                 x = 1,
#                                 dims = c(n0,n1*n0))
#     
#     q_t <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) ,
#                                 j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
#                                 x = 1,
#                                 dims = c(n1,n1*n0))
#     q_s <- rbind(q_c,q_t)
#     # q_s   <- Matrix::sparseMatrix(i = 1:(n1*n0),
#     #                               j = 1:(n1*n0),
#     #                               x = 1,
#     #                               dims = c(n1*n0,n1*n0))
#     # Q0_c <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
#     # Q0_t <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
#     # Q0_c <- (matrix(-1/n0, n0,n0))
#     # Q0_t <- (matrix(-1/n1, n1,n1))
#     # q_m <- Matrix::sparseMatrix(i = c(rep(1:n0, each = n0), rep(n0 + 1:n1, each = n1)),
#     #                             j = c(rep.int(1:n0, n0), rep.int(n0 + 1:n1, n1)),
#     #                             x = c(Q0_c, Q0_t),
#     #                             dims = c(n,n)
#     # )
#     rm(q_t, q_c)
#   }
#   # mean_mat <- as.simple_triplet_matrix(matrix(marg_mass, marg_n, n1*n0))
#   # q_s <- Matrix::sparseMatrix(i = q_mat_list$i,
#   #                             j = q_mat_list$j,
#   #                             x = q_mat_list$val)  #Matrix::Matrix(q_mat, sparse = TRUE)
#   # q_self  <- Matrix::crossprod(q_s)
#   # q_mult  <- Matrix::crossprod(q_m_chol, q_s)
#   # q_cross  <- Matrix::crossprod(q_s)
#   if(!Matrix::isDiagonal(q_s)) {
#     Q0 <-  2*Matrix::crossprod(q_s)
#   } else {
#     Q0 <-  2*q_s
#   }
#   rm(q_s)
#   
#   L0 <- cost_vec #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
#   #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
#   obj <- list(Q = Q0,
#               L = L0)
#   
#   LC <- list()
#   LC$A <- rbind(sum_const, marg_const_mat)
#   LC$vals <- c(1, marg_const)
#   LC$dir <- c(rep("E", 1 + marg_const_n))
#   
#   op <- list(obj=obj, LC=LC)
#   return(op)
# }

qp_wass <- function(x, z, p = 2, target = c("ATC", "ATT", 
                                                     "feasible"),
                          dist=c("Lp", "mahalanobis"), cost = NULL) {
  target <- match.arg(target)
  x1 <- x[z==1,,drop=FALSE]
  x0 <- x[z==0,,drop=FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  # x1_var <- cov(x1)
  # x0_var <- cov(x0)
  
  dist <- match.arg(dist)
  dist.fun <- switch(dist, 
                     "Lp" = causalOT::cost_calc_lp,
                     "mahalanobis" = causalOT::cost_mahalanobis)
  # covar <- (x1_var + x0_var)/2 # add later
  if(is.null(cost)) { 
    cost <- dist.fun(x0, x1, ground_p = p)
  } else {
    stopifnot(all(dim(cost) %in% c(nrow(x1),nrow(x0))))
  }
  # cost <- as.matrix(dist(rbind(x0,x1), method = "minkowski", p = p))[1:n0, (n0+1):n]
  
  cost_vec <- c(cost)^p
  sum_const <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                    j = 1:(n0*n1),
                                    x = rep(1, n1 * n0),
                                    dims = c(1, n0*n1))
  q_s <- marg_const_mat <- marg_const <- NULL
  marg_const_n <- 0
  
  # Qmat
  if(target == "ATT") {
    # q_mat <- matrix(0, n0,n1*n0)
    n1_idx <- seq(1,n0*n1,n0)
    # q_s <- Matrix::sparseMatrix(i = rep(1:n0, each = n1) , 
    #                             j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
    #                             x = 1,
    #                             dims = c(n0,n1*n0))
    # for(i in 1:n0) q_mat[i, (i-1) + n1_idx] <- 1
    marg_mass <- 1/n0
    marg_n <- n0
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n1, each = n0),
                                           j = c(sapply(0:(n1-1), function(i) i*n0 + 1:n0)),
                                           x = rep(1,n0*n1),
                                           dims = c(n1,n0*n1))
    # marg_const_mat <- matrix(0, n1,n1*n0)
    # for(i in 1:n1) marg_const_mat[i, (i-1) * n0 + (1:n0)] <- 1
    marg_const <- rep( 1/n1, n1)
    marg_const_n <- n1
    # q_m <- (diag(1, n0, n0) - matrix(1/n0, n0,n0))
    # q_m <- matrix(-1/n0, n0,n0)
  } else if (target == "ATC") {
    n0_idx <- 1:n0
    # q_s <- Matrix::sparseMatrix(i = rep.int(1:n1, n0) , 
    #                             j = c(sapply(0:(n1-1), function(i) i * n0 + n0_idx)),
    #                             x = 1,
    #                             dims = c(n1,n1*n0))
    
    marg_mass <- 1/n1
    marg_n <- n1
    n1_idx <- seq(1,n0*n1,n0)
    marg_const_mat <- Matrix::sparseMatrix(i = rep(1:n0, each = n1),
                                           j = c(sapply(0:(n0-1), function(i) i + n1_idx)),
                                           x= rep.int(1,n0*n1),
                                           dims = c(n0, n0*n1))
    
    marg_const <- rep.int( 1/n0, n0 )
    marg_const_n <- n0
    # q_m <- (diag(1, n1, n1) - matrix(1/n1, n1,n1))
  } else if (target == "feasible") {
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
  
  L0 <- cost_vec #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  obj <- list(Q = NULL,#Q0,
              L = L0)
  
  LC <- list()
  LC$A <- rbind(sum_const, marg_const_mat)
  LC$vals <- c(1, marg_const)
  LC$dir <- c(rep("E", 1 + marg_const_n))
  
  op <- list(obj=obj, LC=LC)
  return(op)
}

qp_rkhs <- function(x, z, estimate = c("ATT", "ATC", "feasible"), d = 1, theta = c(1,1),
                    gamma = c(1,1), lambda = 1, cost = NULL, ...) {
  # est <- match.arg(estimate)
  # if (est == "ATC") {
  #   z <- 1 - z
  # }
  n <- nrow(x)
  
  dist <- match.arg(dist)
  
  if(is.null(cost)) { 
    cost <- kernel_calculation(x, z, d = d, theta = theta, gamma = gamma)
  } else {
    stopifnot(all(dim(cost) %in% c(nrow(x1),nrow(x0))))
  }
  
  Q0 <- (1/(n^2)) * (cost + diag(lambda/(n^2), n ,n))
  
  L0 <- c(w = -2/(n^2) * c(colMeans(cost)))
  
  A <- rep(1.0/n, n)
  vals <- n
  dir <- "E"
  
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

