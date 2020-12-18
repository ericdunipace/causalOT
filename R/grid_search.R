sbw_grid_search <- function(data, grid = NULL, 
                            estimand = c("ATT", "ATC","cATE","ATE","feasible"),
                            n.boot = 100,
                            ...) 
{
  # if(is.null(grid) & !is.null(list(...)$constraint)) grid <- constraint
  if (all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 1/sqrt(get_p(data)), length.out = 10)
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
  for (g in seq_along(grid)) {
    output[g] <- mean(sapply(bootIdx, mean_bal_grid, weight = weight.list[[g]], 
                             data = mean.bal.dat, estimand = estimand,
                             tx_ind = "z", balance.covariates = colnames(x)))
  }
  if (all(is.na(output))) stop("sbw_grid_search: All grid values generated errors")
  
  min.idx <- which(output == min(output, na.rm = TRUE))
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
  if (opt.hyperparam) {
    opt_args <- list(x = pd$df[,-1], y = pd$df$y, z = pd$z, power = 2:3, estimand = estimand, ...)
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
  if (all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 100, length.out = 11)
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
  for (g in seq_along(grid)) {
    output[g] <- mean(sapply(bootIdx, mean_bal_grid, weight = weight.list[[g]], 
                             data = mean.bal.dat, estimand = estimand,
                             tx_ind = "z", balance.covariates = colnames(x)))
  }
  min.idx <- which(output == min(output, na.rm = TRUE))[1]
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
                             add.joint = FALSE,
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
  if (is.null(solver)) solver <- "gurobi"
  if (is.null(p)) p <- 2
  if (is.null(metric)) metric <- "mahalanobis"
  
  pd <- prep_data(data, ...)
  z  <- pd$z
  
  if (!is.null(pd$df$y)) {
    pd$df$y <- NULL
  }
  x  <- as.matrix(pd$df)
  x1 <- x[z == 1,, drop = FALSE]
  x0 <- x[z == 0,, drop = FALSE]
  
  
  n <- nrow(x)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  wass.dat <- cbind(z = z, x)
  
  cost <- dots$cost
  rerun <- FALSE
  if (method == "Wasserstein" ) {
    if (estimand == "ATE")  {
      if ( length(cost[[1]]) != (ncol(x) + 1) | 
                              length(cost[[2]]) != (ncol(x) + 1)) {
        rerun <- TRUE 
      }
    } else {
      if ( length(cost) != (ncol(x) + 1) ) {
        rerun <- TRUE
      }
    }
  }
  if (is.null(cost) | rerun) {
    # if(is.null(dots$p)) dots$p <- 2
    # if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    
    # if(estimand == "ATE" & metric != "RKHS") {
    #   cost <- list(cost_fun(x0, x, ground_p = p, metric = metric),
    #                cost_fun(x1, x, ground_p = p, metric = metric))
    #     
    # } else {
      
    # }
    cost <- wass_cost_default(x = x, z = z, estimand = estimand, 
                              metric = metric, method = method, p = p, 
                              rkhs.args = dots$rkhs.args)
  } else {
    cost <- dots$cost
  }
  
  if (is.null(grid) | is.na(grid)) {
    gargs <- list(x = x, z = z, p = p, data = data, cost = cost, estimand = estimand,
                  method = method, metric = metric, wass.iter = wass.iter, 
                  add.joint = add.joint, ...)
    gargs <- gargs[!duplicated(names(gargs))]
    n.gargs <- lapply(names(gargs), as.name)
    names(n.gargs) <- names(gargs)
    g.call <- as.call(c(list(as.name("wass_grid_default")), n.gargs))
    grid <- eval(g.call, envir = gargs)
  }
  
  args <- list(data = data, constraint = grid[[1]],  estimand = estimand, 
               method = method, solver = solver, metric = metric,
               p = p, cost = cost, add.joint = add.joint,
               ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(setNames(c(as.name("calc_weight_bal"), argn), c("", names(args))))
  if (verbose) message("\nEstimating Wasserstein values for each constraint")
  weight.list <- lapply(grid, function(delta) {
    args$constraint <- delta
    out <- tryCatch(eval(f.call, envir = args),
                    error = function(e) {return(list(w0 = rep(NA_real_, n0),
                                                     w1 = rep(NA_real_, n1),
                                                     gamma = NULL,
                                                     estimand = estimand,
                                                     method = method, args = list(constraint = delta)))})
    # out <- do.call("calc_weight_bal", args)
    if (solver == "gurobi") Sys.sleep(0.1)
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
                           add.joint = add.joint,
                           ...
  )
  boot.args <- boot.args[!duplicated(names(boot.args))]
  boot.n <- lapply(names(boot.args), as.name)
  names(boot.n) <- names(boot.args)
  b.call <- as.call(c(list(quote(sapply)), boot.n))
  if (verbose ) {
    message("\nCalculating out of sample balance")
    pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
    
  }
  
  for (g in seq_along(grid)) {
    boot.args$weight <- weight.list[[g]]
    output[g] <- mean(eval(b.call, envir = boot.args))
    if (verbose) setTxtProgressBar(pb, g)
  }
  if (all(is.na(output))) stop("wass_grid_search: All grid values generated errors")
  
  min.idx <- which(output == min(output, na.rm = TRUE))
  # weight.list[[min.idx]]$args$standardized.mean.difference <- grid[min.idx]
  return(weight.list[[min.idx]])
}


mean_bal_grid <- function(bootIdx, weight, data, tx_ind, estimand, ...) {
  wvec <- c(weight$w0, weight$w1)
  if (all(is.na(wvec))) return(NA)
  dataResamp <- data[bootIdx,]
  weightResamp <- wvec[bootIdx]
  z <- dataResamp[,tx_ind]
  if (estimand == "ATE" | estimand == "cATE") {
    n <- nrow(data)
    w0 <- list(w0 = renormalize(weightResamp[z == 0]), w1 = rep(1/n,n))
    z0 <- c(rep(0, sum(z == 0)), rep(1,n))
    d0 <- rbind(dataResamp[z == 0, ], dataResamp)
    d0[, tx_ind] <- z0

    w1 <- list(w0 = renormalize(weightResamp[z == 1]), w1 = rep(1/n,n))
    z1 <- c(rep(0, sum(z == 1)), rep(1,n))
    d1 <- rbind(dataResamp[z == 1, ], dataResamp)
    d1[, tx_ind] <- z1
    
    bals <- c(sum(1 - z) / n * mean_bal(d0, weights = w0, treatment.indicator = tx_ind, ...),
              sum(z) / n * mean_bal(d1, weights = w1, treatment.indicator = tx_ind, ...))
    
  } else if (estimand == "ATT" | estimand == "ATC" | estimand == "feasible") {
    wl <- list(w0 = renormalize(weightResamp[z == 0]), w1 = renormalize(weightResamp[z == 1]))
    bals <- mean_bal(dataResamp, weights = wl, treatment.indicator = tx_ind, ...)
  }
  return(mean(bals))
}

wass_dist_helper <- function(...) {
  args <- list(...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  f.call <- as.call(setNames(c(list(quote(wasserstein_p)), argn), c("",names(args))))
  return(eval(f.call, envir = args))
}

wass_cost_default <- function(x, z, estimand, metric, method, p, rkhs.args) {
  cost.fun.args <- list(x = x, z = z, power = p, metric = metric,
                        rkhs.args = rkhs.args,
                        estimand = estimand)
  cost.fun.args <- cost.fun.args[!duplicated(names(cost.fun.args))]
  cf.name <- lapply(names(cost.fun.args), as.name)
  
  cf.call <- as.call(setNames(c(as.name("cost_fun"), cf.name), c("",names(cost.fun.args))))
  
  if (method == "Wasserstein") {
    D <- ncol(x)
    cost.fun.args2 <- cost.fun.args
    cost.temp <- c(lapply(1:D, function(d) {
      cost.fun.args2$x <- x[,d, drop = FALSE]
      return(eval(cf.call, envir = cost.fun.args2))
    }),
    list(eval(cf.call, envir = cost.fun.args)))
    
    if (estimand == "ATE") {
      cost <- list(lapply(cost.temp, function(cc) cc[[1]]),
                 lapply(cost.temp, function(cc) cc[[2]]))
    } else {
      cost <- cost.temp
    }
    
  } else if (method == "Constrained Wasserstein") {
    
    cost <- eval(cf.call, envir = cost.fun.args)
  }
  return(cost)
}


wass.fun.grid <- function(x, z, p, data, cost, estimand, metric, wass.iter, add.joint, ...) {
  D <- ncol(x)
  
  # x0  <- x[z]
  # x <- rbind(x0, x1)
  
  n <- nrow(x)
  
  # z <- c(rep(0, n0), 
  #        rep(1, n1))
  # if (metric == "sdLp") {
  #   div <- 1/matrixStats::colSds(rbind(x))
  # } else {
  #   div <- rep(1, D)
  # }
  # if (metric == "sdLp") {
  #   x <- scale(x, center = FALSE, scale = TRUE)
  # } else if (metric == "mahalanobis") {
  #   x <- mahal_transform(x)
  # } 
  
  x0 <- x[z == 0, , drop = FALSE]
  x1 <- x[z == 1, , drop = FALSE]
  
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  if (estimand == "ATE") {
    
    nnm <- lapply(1:(D + 1), function(d) calc_weight(data, estimand = estimand, method = "NNM", 
                                                     cost = list(cost[[1]][[d]],
                                                                 cost[[2]][[d]]),
                                                     ...))
    
    w0 <- w1 <- nnm
    for (d in 1:(D + 1)) {
      w1[[d]]$w0 <- w1[[d]]$w1
      w1[[d]]$w1 <- w0[[d]]$w1 <- rep(1/n,n)
    }
    
    wass_nnm_0 <- 
      # c(
      sapply(1:(D + 1), 
             function(d) 
               # approxOT::transport_plan(X = t(x0[,d]), 
               #                                         Y = t(x[,d]),
               #                                         a = w0[[d]]$w0,
               #                                         b = w0[[d]]$w1,
               #                                         method = "univariate",
               #                                         observation.orientation = "colwise",
               #                                         ground_p = p, p = p)$cost),
      wass_dist_helper(a = w0[[d]], cost = cost[[1]][[d]], p = p, method = "networkflow", niter = wass.iter, ...) 
      )
    
    wass_nnm_1 <- 
      # c(
      sapply(1:(D + 1), function(d) 
        # approxOT::transport_plan(X =  t(x1[,d]), 
        #                                                Y = t(x[,d]),
        #                                                a = w1[[d]]$w0,
        #                                                b = w1[[d]]$w1,
        #                                                method = "univariate",
        #                                                observation.orientation = "colwise",
        #                                                ground_p = p, p = p)$cost),
      wass_dist_helper(a = w1[[ d ]], cost = cost[[2]][[d]], p = p, method = "networkflow", niter = wass.iter, ...) 
      )
    wass_full_0 <- 
      # c(
      sapply(1:(D + 1), function(d) 
        # approxOT::transport_plan(X = t(x0[,d]), 
        #                                                Y = t(x[,d]),
        #                                                method = "univariate",
        #                                                observation.orientation = "colwise",
        #                                                ground_p = p, p = p)$cost),
      wass_dist_helper(a = rep(1/n0,n0), 
                       b = rep(1/n,n),
                       cost = cost[[1]][[d]], 
                       p = p, method = "networkflow", niter = wass.iter, ...
      ))
    wass_full_1 <- 
      # c(
      sapply(1:(D + 1), function(d)
             # approxOT::transport_plan(X = t(x1[,d]), 
             #                                           Y = t(x[,d]),
             #                                           method = "univariate",
             #                                           observation.orientation = "colwise",
             #                                           ground_p = p, p = p)$cost),
      wass_dist_helper(a = rep(1/n1,n1), 
                       b = rep(1/n, n),
                       cost = cost[[2]][[d]], 
                       p = p, method = "networkflow", niter = wass.iter, ...
      ))
    if (metric == "Lp") {
      scales <- matrixStats::colSds(x)^p
      
    } else {
      scales <- rep(1,D)
    }
    
    #for control
    nnm.adjusted0 <- (wass_nnm_0[1:D]^p * 1/scales)^(1/p)
    full.adjusted0 <- (wass_full_0[1:D]^p * 1/scales)^(1/p)
    min.grid0 <- (min(nnm.adjusted0)^p * scales)^(1/p)
    max.grid0 <- (max(full.adjusted0)^p * scales)^(1/p)
    
    #for treated
    nnm.adjusted1 <- (wass_nnm_1[1:D]^p * 1/scales)^(1/p)
    full.adjusted1 <- (wass_full_1[1:D]^p * 1/scales)^(1/p)
    min.grid1 <- (min(nnm.adjusted1)^p * scales)^(1/p)
    max.grid1 <- (max(full.adjusted1)^p * scales)^(1/p)
    
    if (add.joint) {
      jl.grid0 <- sum(min.grid0^p)^(1/p)
      ju.grid0 <- sum(max.grid0^p)^(1/p)
      jl.grid1 <- sum(min.grid1^p)^(1/p)
      ju.grid1 <- sum(max.grid1^p)^(1/p)
      if (ju.grid0 < wass_full_0[[D + 1]]) ju.grid0 <- wass_full_0[[D + 1]]
      if (jl.grid0 < wass_nnm_0[[D + 1]]) jl.grid0 <- wass_nnm_0[[D + 1]]
      if (ju.grid1 < wass_full_1[[D + 1]]) ju.grid1 <- wass_full_1[[D + 1]]
      if (jl.grid1 > wass_nnm_1[[D + 1]]) jl.grid1 <- wass_nnm_1[[D + 1]]
      
      jgrid0 <- seq(jl.grid0, ju.grid0, length.out = 7)
      jgrid1 <- seq(jl.grid1, ju.grid1, length.out = 7)
      
      grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) seq(x[1], x[2], length.out = 7)))
      grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) seq(x[1], x[2], length.out = 7)))
      
      grid_0 <- do.call("cbind", lapply(jgrid0, function(joint) rbind(grid_1, joint)))
      grid_1 <- do.call("cbind", lapply(jgrid1, function(joint) rbind(grid_1, joint)))
      
    } else {
      grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) seq(x[1], x[2], length.out = 50)))
      grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) seq(x[1], x[2], length.out = 50)))
      
      
    }
    
    grid <- lapply(1:ncol(grid_0), function(i) list(grid_0[,i],
                                          grid_1[,i]))
    
    
    
  } else {
    nnm <- lapply(cost, function(cc) calc_weight(data, estimand = estimand, method = "NNM", 
                                                 cost = cc,
                                                 ...))
    wass_nnm <- sapply(1:(D + 1), function(d) wass_dist_helper(a = nnm[[d]], cost = cost[[d]], p = p, method = "networkflow", niter = wass.iter, ...) )
    wass_full <- 
      # c(
      sapply(1:(D + 1), function(d) 
        # approxOT::transport_plan(X = t(x0[,d]), 
        #                                                Y = t(x1[,d]),
        #                                                method = "univariate",
        #                                                observation.orientation = "colwise",
        #                                                ground_p = p, p = p)$cost),
      wass_dist_helper(a = rep(1/n0,n0), 
                       b = rep(1/n1,n1),
                       cost = cost[[d]], 
                       p = p, method = "networkflow", niter = wass.iter, ...
      ))
    if (metric == "Lp") {
      scales <- matrixStats::colSds(x)^p
      
    } else {
      scales <- rep(1,D)
    }
    nnm.adjusted <- (wass_nnm[1:D]^p * 1/scales)^(1/p)
    full.adjusted <- (wass_full[1:D]^p * 1/scales)^(1/p)
    min.grid <- (min(nnm.adjusted)^p * scales)^(1/p)
    max.grid <- (max(full.adjusted)^p * scales)^(1/p)
    
    if (add.joint) {
      jl.grid <- sum(min.grid^p)^(1/p)
      ju.grid <- sum(max.grid^p)^(1/p)
      if (jl.grid < wass_nnm[[D + 1]]) jl.grid <- wass_nnm[[D + 1]]
      if (ju.grid < wass_full[[D + 1]]) ju.grid <- wass_full[[D + 1]]
      
      jgrid <- seq(jl.grid, ju.grid, length.out = 7)

      grid_temp <- t(apply(cbind(min.grid, max.grid), 1, function(x) seq(x[1], x[2], length.out = 7)))
      
      grid_temp <- do.call("cbind", lapply(jgrid, function(joint) rbind(grid_temp, joint)))
      
    } else {
        grid_temp <- t(apply(cbind(min.grid, max.grid), 1, function(x) seq(x[1], x[2], length.out = 50)))
        
    }
    grid <- lapply(1:ncol(grid_temp), function(i) grid_temp[,i] )
  }
  return(grid)
}

cwass.fun.grid <- function(cost, p, estimand, wass.iter, ...) {
  
  if (estimand == "ATE") {
    n0  <- nrow(cost[[1]])
    n1  <- nrow(cost[[2]])
    n   <- n0 + n1
    nnm <- calc_weight(data, estimand = estimand, method = "NNM", 
                       cost = cost, p = p,
                       ...)
    w0 <- w1 <- nnm
    w1$w0 <- w1$w1
    w1$w1 <- w0$w1 <- rep(1/n,n)
    wass_nnm <- list(
      wass_dist_helper(a = w0, cost = cost[[1]], p = p, method = "networkflow", niter = wass.iter, ...),
      wass_dist_helper(a = w1, cost = cost[[2]], p = p, method = "networkflow", niter = wass.iter, ...)
    )
    wass_full <- list(
      wass_dist_helper(a = rep(1/n0,n0), 
                       b = rep(1/n,n),
                       cost = cost[[1]], 
                       p = p, method = "networkflow", niter = wass.iter, ...),
      wass_dist_helper(a = rep(1/n1,n1), 
                       b = rep(1/n,n),
                       cost = cost[[2]], 
                       p = p, method = "networkflow", niter = wass.iter, ...)
    )
    
    grid <- rbind(seq(wass_nnm[[1]], wass_full[[1]], length.out = 10),
                  seq(wass_nnm[[2]], wass_full[[2]], length.out = 10))
    grid <- lapply(1:10, function(i) grid[,i])
    
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
    n0  <- nrow(cost)
    n1  <- ncol(cost)
    n   <- n0 + n1
    nnm <- calc_weight(data, estimand = estimand, method = "NNM", 
                       cost = cost, p = p,
                       ...)
    wass_nnm <- wass_dist_helper(a = nnm, cost = cost, p = p, method = "networkflow", niter = wass.iter, ...)
    wass_full <- wass_dist_helper(a = rep(1/n0,n0), 
                                  b = rep(1/n1,n1),
                                  cost = cost, 
                                  p = p, method = "networkflow", niter = wass.iter, ...
    )
    
    grid <- seq(wass_nnm, wass_full, length.out = 10)
    
  }
  return(grid)
}

wass_grid_default <- function(x, z, p, data, cost, estimand, method, metric, wass.iter, ...) {
  
    
  get_defaults <- switch(method,
                         "Constrained Wasserstein" = as.name("cwass.fun.grid"),
                         "Wasserstein" = as.name("wass.fun.grid"),
                         "Sliced Wasserstein" = NULL)
  args <- list(x = x, z = z, p = p, data = data, 
               cost = cost, estimand = estimand, metric = metric, 
               wass.iter = wass.iter, ...)
  args <- args[!duplicated(names(args))]
  n.args <- lapply(names(args), as.name)
  names(n.args) <- names(args)
  f.call <- as.call(c(list(get_defaults), args))
  
  return(eval(f.call, envir = args))
}

wass_grid <- function(bootIdx, weight, cost, tx_idx, estimand, wass.method, wass.iter, ...) {
  if (all(is.na(weight$w1) | all(is.na(weight$w0)))) return(NA)
  n1 <- length(weight$w1)
  n0 <- length(weight$w0)
  ord.idx <- rep(NA_integer_, n1 + n0)
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
  
  if (estimand == "ATE") {
    w0 <- w1 <- weight
    w0$gamma <- w1$gamma <- NULL
    w1$w0 <- w1$w1
    w1$w1 <- w0$w1 <- rep(1/n, n)
    if (is.list(cost[[1]]) & is.list(cost[[2]])) {
      c0 <- cost[[1]][[length(cost[[1]])]][cn_boot, bootIdx]
      c1 <- cost[[2]][[length(cost[[2]])]][tx_boot, bootIdx]
    } else {
      c0 <- cost[[1]][cn_boot, bootIdx]
      c1 <- cost[[2]][tx_boot, bootIdx]
    }
    wass0 <- wass_dist_helper(a = w0, b = NULL,
                           cost = c0,
                           p = p, method = wass.method, niter = wass.iter, ...)
    wass1 <- wass_dist_helper(a = w1, b = NULL,
                           cost = c1,
                           p = p, method = wass.method, niter = wass.iter, ...)
    return(((n1 * wass0^p + wass1^p * n0)/n)^(1/p))
  } else {
    # if(!is.null(weight$gamma)) weight$gamma <- weight$gamma[cn_boot, tx_boot]
    if ( is.list(cost) ) {
      cc <- cost[[length(cost)]][cn_boot, tx_boot]
    } else {
      cc <- cost[cn_boot, tx_boot]
    }
    return(wass_dist_helper(a = weight, b = NULL,
                         cost = cc,
                         p = p, method = wass.method, niter = wass.iter, ...))
  }
}
