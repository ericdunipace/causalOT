sbw_grid_search <- function(data, grid = NULL, 
                            estimand = c("ATT", "ATC","cATE","ATE","feasible"),
                            n.boot = 100,
                            ...) 
{
  # if(is.null(grid) & !is.null(list(...)$constraint)) grid <- constraint
  if(all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 1/sqrt(get_p(data)), length.out = 10)
  estimand <- match.arg(estimand)
  
  args <- list(data = data, constraint = grid[1],  estimand = estimand, 
               method = "SBW",
               ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(setNames(c(as.name("calc_weight_bal"), argn), c("", names(args))))
  
  pd <- prep_data(data, ...)
  x1 <- as.matrix(pd$df[pd$z == 1,-1])
  x0 <- as.matrix(pd$df[pd$z == 0,-1])
  x <- rbind(x0,x1)
  
  
  n <- nrow(pd$df)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  mean.bal.dat <- cbind(z = c(rep(0,n0),rep(1,n1)), x)
  
  weight.list <- lapply(grid, function(delta) {
    args$constraint <- delta
    out <- tryCatch(eval(f.call, envir = args),
                    error = function(e) {return(list(w0 = rep(NA_real_, n0),
                                                     w1 = rep(NA_real_, n1)))})
    # out <- do.call("calc_weight_bal", args)
    if(list(...)$solver == "gurobi") Sys.sleep(0.1)
    return(out)
  })
  pd <- prep_data(data, ...)
  x1 <- as.matrix(pd$df[pd$z == 1,-1])
  x0 <- as.matrix(pd$df[pd$z == 0,-1])
  x <- rbind(x0,x1)
  mean.bal.dat <- cbind(z = pd$z, x)
  
  n <- nrow(pd$df)
  
  bootIdx <- lapply(1:n.boot, function(ii) {sample.int(n,n, replace = TRUE)})
  output <- rep(NA, length(grid))
  names(output) <- as.character(grid)
  for(g in seq_along(grid)) {
    output[g] <- mean(sapply(bootIdx, mean_bal_grid, weight = weight.list[[g]], 
                             data = mean.bal.dat, tx_ind = "z", balance.covariates = colnames(x)))
  }
  if(all(is.na(output))) stop("sbw_grid_search: All grid values generated errors")
  
  min.idx <- which(output == min(output, na.rm=TRUE))
  weight.list[[min.idx]]$args$standardized.mean.difference <- grid[min.idx]
  return(weight.list[[min.idx]])
}

RKHS_grid_search <- function(data, grid = NULL, 
                             method = c("RKHS.dose"),
                             estimand = c("ATT", "ATC","ATE"),
                             n.boot = 100, opt.hyperparam = TRUE,
                             ...) 
{
  meth <- match.arg(method)
  estimand <- match.arg(estimand)
  pd <- prep_data(data, ...)
  if(opt.hyperparam) {
    opt_args <- list(x=pd$df[,-1], y = pd$df$y, z = pd$z, power = 2:3, estimand = estimand, ...)
    opt_args <- opt_args[!duplicated(names(opt_args))]
    param <- do.call("RKHS_param_opt", opt_args)
    rm(opt_args)
    
    args <- c(list(data = data,
                   lambda = 0,
                   method = meth,
                   estimand = estimand),
              param,
              opt.hyperparam = FALSE,
              list(...))
  } else {
    args <- list(data = data,
                 lambda = 0,
                 method = meth,
                 estimand = estimand,
                 opt.hyperparam = FALSE,
                 ...)
  }
  args <- args[!duplicated(names(args))]
  if(all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 100, length.out = 11)
  weight.list <- lapply(grid, function(delta) {
    args$lambda <- delta
    out <- do.call("calc_weight_RKHS", args)
    return(out)
  })
  
  x1 <- as.matrix(pd$df[pd$z == 1,-1])
  x0 <- as.matrix(pd$df[pd$z == 0,-1])
  x <- rbind(x0,x1)
  mean.bal.dat <- cbind(z = pd$z, x)
  
  n <- nrow(pd$df)
  
  bootIdx <- lapply(1:n.boot, function(ii) {sample.int(n,n, replace = TRUE)})
  output <- rep(NA, length(grid))
  names(output) <- as.character(grid)
  for(g in seq_along(grid)) {
    output[g] <- mean(sapply(bootIdx, mean_bal_grid, weight = weight.list[[g]], 
                             data = mean.bal.dat, tx_ind = "z", balance.covariates = colnames(x)))
  }
  min.idx <- which(output == min(output, na.rm=TRUE))[1]
  weight.list[[min.idx]]$estimand <- "ATE"
  weight.list[[min.idx]]$lambda <- grid[min.idx]
  weight.list[[min.idx]]$args$grid.search <- TRUE
  return(weight.list[[min.idx]])
}


wass_grid_search <- function(data, grid = NULL, 
                             estimand = c("ATT", "ATC","cATE","ATE"),
                             n.boot = 100,
                             method = c("Wasserstein","Constrained Wasserstein"),
                             wass.method = "shortsimplex", wass.iter = 0,
                             verbose = FALSE,
                             ...) 
{
  # if(is.null(grid) & !is.null(list(...)$constraint)) grid <- constraint
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  dots <- list(...)
  solver <- dots$solver
  p      <- dots$p
  metric <- dots$metric
  if(is.null(solver)) solver <- "gurobi"
  if(is.null(p)) p <- 2
  if(is.null(metric)) metric <- "mahalanobis"
  
  pd <- prep_data(data, ...)
  z  <- pd$z
  
  if(!is.null(pd$df$y)) {
    pd$df$y <- NULL
  }
  x  <- as.matrix(pd$df)
  x1 <- x[z == 1,, drop = FALSE]
  x0 <- x[z == 0,, drop = FALSE]
  
  
  n <- nrow(x)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  wass.dat <- cbind(z = z, x)
  
  if(is.null(dots$cost)) {
    # if(is.null(dots$p)) dots$p <- 2
    # if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    
    # if(estimand == "ATE" & metric != "RKHS") {
    #   cost <- list(cost_fun(x0, x, ground_p = p, metric = metric),
    #                cost_fun(x1, x, ground_p = p, metric = metric))
    #     
    # } else {
      cost.fun.args <- list(x=x, z=z, power = p, metric = metric,
                            rkhs.args = dots$rkhs.args,
                            estimand = estimand, ...)
      cost.fun.args <- cost.fun.args[!duplicated(names(cost.fun.args))]
      cf.name <- lapply(names(cost.fun.args), as.name)
      
      cf.call <- as.call(setNames(c(as.name("cost_fun"), cf.name), c("",names(cost.fun.args))))
      cost <- eval(cf.call, envir = cost.fun.args)
    # }
    
  } else {
    cost <- dots$cost
  }
  
  if(estimand == "ATE") {
    if(all(is.null(grid)) | all(is.na(grid))) {
      
      nnm <- calc_weight(data, estimand = estimand, method = "NNM", 
                         cost = cost,
                         ...)
      w0 <- w1 <- nnm
      w1$w0 <- w1$w1
      w1$w1 <- w0$w1 <- rep(1/n,n)
      wass_nnm <- list(
        wass_dist_helper(a= w0, cost =cost[[1]], p = p, method = "networkflow", niter = wass.iter, ...),
        wass_dist_helper(a= w1, cost =cost[[2]], p = p, method = "networkflow", niter = wass.iter, ...)
      )
      wass_full <- list(
        wass_dist_helper(a=rep(1/n0,n0), 
                      b = rep(1/n,n),
                      cost = cost[[1]], 
                      p = p, method = "networkflow", niter = wass.iter, ...),
        wass_dist_helper(a=rep(1/n1,n1), 
                      b = rep(1/n,n),
                      cost = cost[[2]], 
                      p = p, method = "networkflow", niter = wass.iter, ...)
      )
      
      grid <- rbind(seq(wass_nnm[[1]], wass_full[[1]], length.out = 10),
                    seq(wass_nnm[[2]], wass_full[[2]], length.out = 10))
      grid <- lapply(1:10, function(i) grid[,i])
    }
    
    # args <- list(data = data, constraint = grid[[1]],  estimand = estimand, 
    #              method = method,
    #              ...)
    # args <- args[!duplicated(names(args))]
    # argn <- lapply(names(args), as.name)
    # names(argn) <- names(args)
    # 
    # f.call <- as.call(c(as.name("calc_weight_bal"), argn))
    # 
    # weight.list <- lapply(grid, function(delta) {
    #   args$constraint <- delta
    #   out <- tryCatch(eval(f.call, envir = args),
    #                   error = function(e) {return(list(w0 = rep(NA_real_, n0),
    #                                                    w1 = rep(NA_real_, n1),
    #                                                    gamma = NULL,
    #                                                    estimand = estimand,
    #                                                    method = method, args = list()))})
    #   # out <- do.call("calc_weight_bal", args)
    #   if(solver == "gurobi") Sys.sleep(0.1)
    #   return(out)
    # })
  } else {
    if(all(is.null(grid)) | all(is.na(grid))) {
      
      nnm <- calc_weight(data, estimand = estimand, method = "NNM", 
                         cost = cost,
                         ...)
      wass_nnm <- wass_dist_helper(a = nnm, cost = cost, p = p, method = "networkflow", niter = wass.iter, ...)
      wass_full <- wass_dist_helper(a = rep(1/n0,n0), 
                                 b = rep(1/n1,n1),
                                 cost = cost, 
                                 p = p, method = "networkflow", niter = wass.iter, ...
      )
      
      grid <- seq(wass_nnm, wass_full, length.out = 10)
    }
    
    
    
  }
  args <- list(data = data, constraint = grid[[1]],  estimand = estimand, 
               method = method, solver = solver, metric = metric,
               p = p, cost = cost,
               ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(setNames(c(as.name("calc_weight_bal"), argn), c("", names(args))))
  if(verbose) message("\nEstimating Wasserstein values for each constraint")
  weight.list <- lapply(grid, function(delta) {
    args$constraint <- delta
    out <- tryCatch(eval(f.call, envir = args),
                    error = function(e) {return(list(w0 = rep(NA_real_, n0),
                                                     w1 = rep(NA_real_, n1),
                                                     gamma = NULL,
                                                     estimand = estimand,
                                                     method = method, args = list(constraint = delta)))})
    # out <- do.call("calc_weight_bal", args)
    if(solver == "gurobi") Sys.sleep(0.1)
    class(out) <- "causalWeights"
    return(out)
  })
  
  bootIdx <- lapply(1:n.boot, function(ii) {sample.int(n,n, replace = TRUE)})
  output <- rep(NA, length(grid))
  names(output) <- as.character(grid)
  tx_idx <- which(pd$z == 1)
  boot.args <- list(X = bootIdx, FUN = wass_grid, weight = weight.list[[1]], 
                           data = wass.dat, tx_idx = tx_idx, cost = cost,
                           estimand = estimand,
                           wass.method = wass.method, wass.iter = wass.iter,
                           ...
  )
  boot.args <- boot.args[!duplicated(boot.args)]
  boot.n <- lapply(names(boot.args), as.name)
  names(boot.n) <- names(boot.args)
  f.call <- as.call(c(list(quote(sapply)), boot.n))
  if (verbose ) {
    message("\nCalculating out of sample balance")
    pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
    
  }
  
  for(g in seq_along(grid)) {
    boot.args$weight <- weight.list[[g]]
    output[g] <- mean(eval(f.call, envir = boot.args))
    if (verbose) setTxtProgressBar(pb, g)
  }
  if(all(is.na(output))) stop("wass_grid_search: All grid values generated errors")
  
  min.idx <- which(output == min(output, na.rm=TRUE))
  # weight.list[[min.idx]]$args$standardized.mean.difference <- grid[min.idx]
  return(weight.list[[min.idx]])
}


mean_bal_grid <- function(bootIdx, weight, data, tx_ind,...) {
  wvec <- c(weight$w0, weight$w1)
  if(all(is.na(wvec))) return(NA)
  dataResamp <- data[bootIdx,]
  weightResamp <- wvec[bootIdx]
  z <- dataResamp[,tx_ind]
  wl <- list(w0 = renormalize(weightResamp[z == 0]), w1 = renormalize(weightResamp[z == 1]))
  bals <- mean_bal(dataResamp, weights = wl, treatment.indicator = tx_ind, ...)
  return(mean(bals))
}

wass_dist_helper <- function(...) {
  args <- list(...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  f.call <- as.call(setNames(c(list(quote(wasserstein_p)), argn), c("",names(args))))
  return(eval(f.call, envir = args))
}

wass_grid <- function(bootIdx, weight, cost, tx_idx, estimand, wass.method, wass.iter, ...) {
  if(all(is.na(weight$w1) | all(is.na(weight$w0)))) return(NA)
  n1 <- length(weight$w1)
  n0 <- length(weight$w0)
  ord.idx <- rep(NA_integer_, n1+n0)
  ord.idx[tx_idx] <- 1:n1
  ord.idx[-tx_idx] <- 1:n0
  tx_boot <- ord.idx[bootIdx[bootIdx %in% tx_idx]]
  cn_boot <- ord.idx[bootIdx[!(bootIdx %in% tx_idx)]]
  n1 <- length(tx_boot)
  n0 <- length(cn_boot)
  n <- n0 + n1
  weight$w0 <- renormalize(weight$w0[cn_boot])
  weight$w1 <- renormalize(weight$w1[tx_boot])
  p <- weight$args$power
  
  if(estimand == "ATE") {
    w0 <- w1 <- weight
    w0$gamma <- w1$gamma <- NULL
    w1$w0 <- w1$w1
    w1$w1 <- w0$w1 <- rep(1/n, n)
    wass0 <- wass_dist_helper(a = w0, b = NULL,
                           cost = cost[[1]][cn_boot, bootIdx],
                           p = p, method = wass.method, niter = wass.iter, ...)
    wass1 <- wass_dist_helper(a = w1, b = NULL,
                           cost = cost[[2]][tx_boot, bootIdx],
                           p = p, method = wass.method, niter = wass.iter, ...)
    return(((n1 * wass0^p + wass1^p * n0)/n)^(1/p))
  } else {
    # if(!is.null(weight$gamma)) weight$gamma <- weight$gamma[cn_boot, tx_boot]
    return(wass_dist_helper(a = weight, b = NULL,
                         cost = cost[cn_boot, tx_boot],
                         p = p, method = wass.method, niter = wass.iter, ...))
  }
}
