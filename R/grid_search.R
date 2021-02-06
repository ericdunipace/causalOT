sbw_grid_search <- function(data, grid = NULL, 
                            estimand = c("ATT", "ATC","cATE","ATE","feasible"),
                            n.boot = 1000, grid.length = 10,
                            ...) 
{
  # if(is.null(grid) & !is.null(list(...)$constraint)) grid <- constraint
  if (all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 1/sqrt(get_p(data, ...)), length.out = grid.length)
  estimand <- match.arg(estimand)
  solver <- match.arg(list(...)$solver, c("mosek","gurobi","cplex"))
  
  args <- list(data = data, constraint = grid[1],  estimand = estimand, 
               method = "SBW", solver = solver,
               ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(setNames(c(as.name("calc_weight_bal"), argn), c("", names(args))))
  
  pd <- extract_x(data, ...)
  x0 <- as.matrix(pd$x0)
  x1 <- as.matrix(pd$x1)
  x <- rbind(x0,x1)

  n <- nrow(x)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  weight.list <- lapply(grid, function(delta) {
    args$constraint <- delta
    out <- tryCatch(eval(f.call, envir = args),
                    error = function(e) {
                      warning(e$message)
                      return(list(w0 = rep(NA_real_, n0),
                                                     w1 = rep(NA_real_, n1)))})
    if (solver == "gurobi") Sys.sleep(0.1)
    return(out)
  })
  
  
  if (estimand != "ATE") {
    mean.bal.dat <- cbind(z = c(rep(0,n0), rep(1,n1)), 
                          x)
    bootIdx       <- lapply(1:n.boot, function(ii) 
    {sample.int(n,n, replace = TRUE)})
    
    output        <- rep(NA, length(grid))
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
    
  } else if (estimand == "ATE") {
    
    output_0 <- output_1 <- rep(NA, length(grid))
    names(output_0) <- names(output_1) <- as.character(grid)
    mean.bal.dat.0  <- cbind(z = c(rep(0, n0), rep(1,n)),
                            rbind(x0,x))
    mean.bal.dat.1  <- cbind(z = c(rep(0, n1), rep(1,n)),
                           rbind(x1,x))
    # bootIdx       <- lapply(1:n.boot, function(ii) 
    #     {sample.int(n,n, replace = TRUE)})
    bootIdx.0       <- lapply(1:n.boot, function(ii) 
        {c(sample.int(n0,n0, replace = TRUE), n0 + sample.int(n,n, replace = TRUE))})
    bootIdx.1       <- lapply(1:n.boot, function(ii) 
        {c(sample.int(n1,n1, replace = TRUE), n1 + sample.int(n,n, replace = TRUE))})
    
    full.sample.wt <- rep(1/n,n)
    
    for (g in seq_along(grid)) {
      w0 <- list(w0 = weight.list[[g]]$w0, w1 = full.sample.wt)
      w1 <- list(w0 = weight.list[[g]]$w1, w1 = full.sample.wt)
      output_0[g] <- mean( sapply(bootIdx.0, mean_bal_grid, weight = w0, 
                               data = mean.bal.dat.0, estimand = "ATT",
                               tx_ind = "z", balance.covariates = colnames(x)))
      output_1[g] <- mean( sapply(bootIdx.1, mean_bal_grid, weight = w1, 
                               data = mean.bal.dat.1, estimand = "ATT",
                               tx_ind = "z", balance.covariates = colnames(x)))
    }
    if (all(is.na(output_0)) | all(is.na(output_1))) stop("All grid values generated errors")
    
    min.idx.0 <- which(output_0 == min(output_0, na.rm = TRUE))
    min.idx.1 <- which(output_1 == min(output_1, na.rm = TRUE))
    output.weight <- weight.list[[min.idx.0]]
    output.weight$w1 <- weight.list[[min.idx.1]]$w1
    output.weight$args$constraint <- output.weight$args$standardized.mean.difference <- c(grid[min.idx.0], grid[min.idx.1])
    return(output.weight)
    
  } else {
    stop("Estimand not recognized")
  }
  
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
                             grid.length = 7,
                             estimand = c("ATT", "ATC","cATE","ATE"),
                             n.boot = 1000,
                             method = c("Wasserstein","Constrained Wasserstein"),
                             sample_weight = NULL,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE,
                             add.margins = FALSE,
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
  add.margins <- isTRUE(add.margins)
  add.joint  <- isTRUE(add.joint)
  
  if (!add.margins & !add.joint & method == "Constrained Wasserstein") {
    stop("Must have marginal or joint constraints or both for Constrained Wasserstein")
  }
  
  pd <- prep_data(data, ...)
  z  <- pd$z
  sample_weight <- get_sample_weight(sample_weight, z)
  
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
  if (add.margins) {
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
                              rkhs.args = dots$rkhs.args, 
                              add.margins = add.margins)
  } else {
    cost <- dots$cost
  }
  
  if (all(is.null(grid)) | all(is.na(grid))) {
    gargs <- list(x = x, z = z, 
                  grid.length = grid.length,
                  p = p, data = data, cost = cost, estimand = estimand,
                  method = method, metric = metric, wass.iter = wass.iter, 
                  add.joint = add.joint, add.margins = add.margins,
                  ...)
    gargs <- gargs[!duplicated(names(gargs))]
    n.gargs <- lapply(names(gargs), as.name)
    names(n.gargs) <- names(gargs)
    g.call <- as.call(c(list(as.name("wass_grid_default")), n.gargs))
    grid <- eval(g.call, envir = gargs)
  }
  
  if (!all(is.na(grid)) && isTRUE(all.equal(grid[1], grid[length(grid)], check.attributes = FALSE)) ) {
    return(calc_weight_NNM(data = data, estimand = estimand,
                           transport.matrix = TRUE,
                           ...))
  }
  
  args <- list(data = data, constraint = grid[[1]],  estimand = estimand, 
               method = method, solver = solver, metric = metric,
               p = p, cost = cost, add.joint = add.joint,
               add.margins = add.margins, 
               # save.solution = TRUE,
               # sol = NULL,
               ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(setNames(c(as.name("calc_weight_bal"), argn), c("", names(args))))
  if (verbose) {
    message("\nEstimating Optimal Transport weights for each constraint")
    pbapply::pboptions(type = "timer", style = 3, char = "=")
  } else {
    pbapply::pboptions(type = "none")
  }
  
  # qp.constructor <- switch(args$method,
  #                          "qp_wass",
  #                          "qp_const_wass")
  # 
  # f.call <- as.call(setNames(c(as.name(qp.constructor), argn), c("", names(args))))
  # qp  <- eval(f.call, envir = args)
  # 
  # for(g in grid) {
  #   if (args$method == "Wasserstein") {
  #     qp$Q <- qp$Q * g[[1]][length(g[[1]])]
  #     if (args$add.margins) {
  #       
  #     } 
  #   }
  # }
  
  # out <- vector("list", length(grid))
  # 
  # for (g in seq_along(grid)) {
  #   args$constraint <- unlist(grid[[g]])
  #   out[[g]] <- tryCatch(eval(f.call, envir = args),
  #                               error = function(e) {
  #                                 warning(e$message)
  #                                 return(list(w0 = rep(NA_real_, n0),
  #                                             w1 = rep(NA_real_, n1),
  #                                             gamma = NULL,
  #                                             estimand = estimand,
  #                                             method = method, args = list(constraint = delta)))})
  #   if (!is.null(out[[g]]$args$sol)) {
  #     args$sol <- out[[g]]$args$sol
  #   }
  # }
  
  weight.list <- pbapply::pblapply(grid, function(delta) {
    args$constraint <- delta
    out <- tryCatch(eval(f.call, envir = args),
                    error = function(e) {
                      warning(e$message)
                      return(list(w0 = rep(NA_real_, n0),
                                                     w1 = rep(NA_real_, n1),
                                                     gamma = NULL,
                                                     estimand = estimand,
                                                     method = method, args = list(constraint = delta)))})
    # out <- do.call("calc_weight_bal", args)
    if (solver == "gurobi") Sys.sleep(0.1)
    class(out) <- "causalWeights"
    return(out)
  })
  
  if (estimand != "ATE") {
    # bootIdx.rows <- lapply(1:n.boot, function(ii) {sample.int(n0,n0, replace = TRUE)})
    # bootIdx.cols <- lapply(1:n.boot, function(ii) {sample.int(n1,n1, replace = TRUE)})
    if (wass.method == "networkflow" | wass.method == "exact" | wass.method == "hilbert") {
      rowCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n0), prob = sample_weight$a)), simplify = FALSE)
      colCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n1), prob = sample_weight$b)), simplify = FALSE)
    } else {
      rowCount <- replicate(n.boot, c(rmultinom(1, n0, prob = sample_weight$a)), simplify = FALSE)
      colCount <- replicate(n.boot, c(rmultinom(1, n1, prob = sample_weight$b)), simplify = FALSE)
    }
    
    if (is.null(dots$cost_a)) {
      cost_a <- cost_fun(rbind(x0,x0), z = c(rep(1,n0),
                                   rep(0,n0)),
               power = p,
               metric = metric, estimand = "ATT")
    } else {
      cost_a <- dots$cost_a
      if (is.list(cost_a)) cost_a <- cost_a[[1]]
    }
    
    if (is.null(dots$cost_b)) {
      cost_b <- cost_fun(rbind(x1,x1), z = c(rep(1, n1),
                                             rep(0, n1)),
                         power = p,
                         metric = metric, estimand = "ATT")
    } else {
      cost_b <- dots$cost_b
      if (is.list(cost_b)) cost_b <- cost_b[[1]]
    }
    
    tx_idx <- which(pd$z == 1)
    boot.args <- list(FUN = wass_grid, 
                      # bootIdx.row = bootIdx.rows, 
                      # bootIdx.col = bootIdx.cols,
                      rowCount = rowCount, 
                      colCount = colCount,
                      MoreArgs    = list(weight = weight.list[[1]], 
                          data = wass.dat, 
                          # tx_idx = tx_idx, 
                          cost = cost,
                          p = p,
                          metric = metric,
                      # estimand = estimand,
                          wass.method = wass.method, 
                          wass.iter = wass.iter,
                          add.joint = add.joint,
                          x0 = x0,
                          x1 = x1,
                          cost_a = cost_a,
                          cost_b = cost_b,
                          ...)
    )
    boot.args <- boot.args[!duplicated(names(boot.args))]
    boot.args$MoreArgs <- boot.args$MoreArgs[!duplicated(names(boot.args$MoreArgs))]
    boot.n <- lapply(names(boot.args), as.name)
    names(boot.n) <- names(boot.args)
    b.call <- as.call(c(list(quote(mapply)), boot.n))
    
    output <- rep(NA, length(grid))
    names(output) <- as.character(grid)
    
    if (verbose ) {
      message("\nCalculating out of sample balance:")
      # pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
      pbapply::pboptions(type = "timer", style = 3, char = "=")
    } else {
      pbapply::pboptions(type = "none")
    }
    output <- pbapply::pbsapply(weight.list, function(ww) {
      boot.args$MoreArgs$weight <- ww
      return(mean(eval(b.call, envir = boot.args)^p))
    })
    # for (g in seq_along(grid)) {
    #   boot.args$MoreArgs$weight <- weight.list[[g]]
    #   output[g] <- mean(eval(b.call, envir = boot.args))
    #   if (verbose) setTxtProgressBar(pb, g)
    # }
    pbapply::pboptions(type = "none")
    if (all(is.na(output))) stop("wass_grid_search: All grid values generated errors")
    
    min.idx <- which(output == min(output, na.rm = TRUE))[1]

    return(weight.list[[min.idx]])
  } else {
    # bootIdx.cols    <- lapply(1:n.boot, function(ii) {sample.int(n,n, replace = TRUE)})
    # bootIdx.rows.0  <- lapply(1:n.boot, function(ii) {sample.int(n0,n0, replace = TRUE)})
    # bootIdx.rows.1  <- lapply(1:n.boot, function(ii) {sample.int(n1,n1, replace = TRUE)})
    if (wass.method == "networkflow" | wass.method == "exact" | wass.method == "hilbert") {
      rowCount.0 <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n0), prob = sample_weight$a)),     simplify = FALSE)
      rowCount.1 <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n1), prob = sample_weight$b)),     simplify = FALSE)
      colCount   <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n), prob = sample_weight$total)), simplify = FALSE)
    } else {
      rowCount.0 <- replicate(n.boot, c(rmultinom(1, n0, prob = sample_weight$a)),     simplify = FALSE)
      rowCount.1 <- replicate(n.boot, c(rmultinom(1, n1, prob = sample_weight$b)),     simplify = FALSE)
      colCount   <- replicate(n.boot, c(rmultinom(1, n, prob = sample_weight$total)), simplify = FALSE)
    }
    
    
    
    if (is.null(dots$cost_a)) {
      cost_a <- list(cost_fun(rbind(x0,x0), z = c(rep(1,n0),
                                             rep(0,n0)),
                         power = p,
                         metric = metric, estimand = "ATT"),
                     cost_fun(rbind(x1,x1), z = c(rep(1,n1),
                                                  rep(0,n1)),
                              power = p,
                              metric = metric, estimand = "ATT")
      )
    } else {
      cost_a <- dots$cost_a
      stopifnot(length(cost_a) == 2)
    }
    
    if (is.null(dots$cost_b)) {
      cost_b <- list(cost_fun(rbind(x,x), z = c(rep(1,n),
                                                rep(0,n)),
                              power = p,
                              metric = metric, estimand = "ATT"))
    } else {
      cost_b <- dots$cost_b
      if (is.list(cost_b)) cost_b <- cost_b[[1]]
    }
    
    output_0        <- output_1 <- rep(NA, length(grid))
    names(output_0) <- names(output_1) <- as.character(grid)
    tx_idx          <- which(pd$z == 1)
    full.sample.weight <- rep(1/n, n)
    
    boot.args <- list(FUN = wass_grid,
                      # bootIdx.row  = bootIdx.rows.0,  
                      # bootIdx.col = bootIdx.cols,
                      rowCount  = rowCount.0,  
                      colCount = colCount,
                      MoreArgs = list(
                        weight = weight.list[[1]], 
                        data = wass.dat, 
                        # tx_idx = tx_idx, 
                        cost = cost[[1]],
                        p = p,
                        metric = metric,
                        x0 = x0,
                        x1 = x,
                        cost_a = cost_a[[1]],
                        cost_b = cost_b,
                        # estimand = "ATT",
                        wass.method = wass.method, wass.iter = wass.iter,
                        add.joint = add.joint,
                        ...)
    )
    boot.args <- boot.args[!duplicated(names(boot.args))]
    boot.args$MoreArgs <- boot.args$MoreArgs[!duplicated(names(boot.args$MoreArgs))]
    boot.n <- lapply(names(boot.args), as.name)
    names(boot.n) <- names(boot.args)
    b.call <- as.call(c(list(quote(mapply)), boot.n))
    if (verbose ) {
      message("\nCalculating out of sample balance")
      pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
      
    }
    boot.args1 <- boot.args0 <- boot.args
    boot.args0$MoreArgs$cost   <- cost[[1]]
    boot.args1$MoreArgs$cost   <- cost[[2]]
    
    # boot.args0$bootIdx.row     <- bootIdx.rows.0
    # boot.args1$bootIdx.row     <- bootIdx.rows.1
    boot.args0$rowCount     <- rowCount.0
    boot.args1$rowCount     <- rowCount.1
    boot.args1$x0           <- x1
    boot.args1$cost_a       <- cost_a[[2]]
    
    w0 <- w1 <- weight.list[[length(weight.list)]]
    w0$w1 <- w1$w1 <- full.sample.weight
    w1$w0 <- weight.list[[length(weight.list)]]$w1
    
    w0$args$power <- w1$args$power <- weight.list[[length(weight.list)]]$args$power
    boot.args0$MoreArgs$weight <- w0
    boot.args1$MoreArgs$weight <- w1
    
    for (g in seq_along(grid)) {
 
      boot.args0$MoreArgs$weight$w0 <- weight.list[[g]]$w0
      # boot.args$MoreArgs$cost   <- cost[[1]]
      output_0[g]               <- mean(eval(b.call, envir = boot.args0)^p)
      
     
      boot.args1$MoreArgs$weight$w0 <- weight.list[[g]]$w1
      # boot.args$MoreArgs$cost   <- cost[[2]]
      output_1[g]               <- mean(eval(b.call, envir = boot.args1)^p)
      if (verbose) setTxtProgressBar(pb, g)
    }
    if (all(is.na(output_0)) | all(is.na(output_1))) stop("wass_grid_search: All grid values generated errors")
    
    min.idx.0 <- which(output_0 == min(output_0, na.rm = TRUE))
    min.idx.1 <- which(output_1 == min(output_1, na.rm = TRUE))
    
    output.weight <- weight.list[[min.idx.0]]
    output.weight$w1 <- weight.list[[min.idx.1]]$w1
    output.weight$args$constraint <- list(weight.list[[min.idx.0]]$args$constraint[1],
                                       weight.list[[min.idx.1]]$args$constraint[2])
    # weight.list[[min.idx]]$args$standardized.mean.difference <- grid[min.idx]
    return(output.weight)
  }
  
}


mean_bal_grid <- function(bootIdx, weight, data, tx_ind, ...) {
  wvec <- c(weight$w0, weight$w1 )
  if (all(is.na(wvec))) return(NA)
  dataResamp <- data[bootIdx,]
  weightResamp <- wvec[bootIdx]
  z <- dataResamp[,tx_ind]
  # if (estimand == "ATE" | estimand == "cATE") {
  #   n <- nrow(data)
  #   w0 <- list(w0 = renormalize(weightResamp[z == 0]), w1 = rep(1/n,n))
  #   z0 <- c(rep(0, sum(z == 0)), rep(1,n))
  #   d0 <- rbind(dataResamp[z == 0, ], dataResamp)
  #   d0[, tx_ind] <- z0
  # 
  #   w1 <- list(w0 = renormalize(weightResamp[z == 1]), w1 = rep(1/n,n))
  #   z1 <- c(rep(0, sum(z == 1)), rep(1,n))
  #   d1 <- rbind(dataResamp[z == 1, ], dataResamp)
  #   d1[, tx_ind] <- z1
  #   
  #   bals <- c(sum(1 - z) / n * mean_bal(d0, weights = w0, treatment.indicator = tx_ind, ...),
  #             sum(z) / n * mean_bal(d1, weights = w1, treatment.indicator = tx_ind, ...))
  #   
  # } else if (estimand == "ATT" | estimand == "ATC" | estimand == "feasible") {
    wl <- list(w0 = renormalize(weightResamp[z == 0]), w1 = renormalize(weightResamp[z == 1]))
    bals <- mean_bal(dataResamp, weights = wl, treatment.indicator = tx_ind, ...)
  # }
  return(mean(bals))
}

wass_dist_helper <- function(...) {
  args <- list(...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  # names(argn) <- names(args)
  f.call <- as.call(setNames(c(list(quote(wasserstein_p)), argn), c("",names(args))))
  return(eval(f.call, envir = args))
}

wass_cost_default <- function(x, z, estimand, metric, method, p, rkhs.args,
                              add.margins) {
  cost.fun.args <- list(x = x, z = z, power = p, metric = metric,
                        rkhs.args = rkhs.args,
                        estimand = estimand)
  cost.fun.args <- cost.fun.args[!duplicated(names(cost.fun.args))]
  cf.name <- lapply(names(cost.fun.args), as.name)
  
  cf.call <- as.call(setNames(c(as.name("cost_fun"), cf.name), c("",names(cost.fun.args))))
  
  if (add.margins) {
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
    
  } else {
    
    cost <- eval(cf.call, envir = cost.fun.args)
  } 
  return(cost)
}


marg.cwass.fun.grid <- function(x, z, grid.length, p, data, cost, estimand, metric, wass.iter, add.joint, ...) {
  D <- ncol(x)
  
  D_plus <- if (add.joint) {
    D + 1
  } else {
    D
  }
  if (wass.iter < 1000 & wass.iter != 0) wass.iter <- wass.iter * 50
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
    
    
    nnm <- lapply(1:D_plus, function(d) calc_weight(data, estimand = estimand, method = "NNM", 
                                                     cost = list(cost[[1]][[d]],
                                                                 cost[[2]][[d]]),
                                                     ...))
    
    w0 <- w1 <- nnm
    for (d in 1:D_plus) {
      w1[[d]]$w0 <- w1[[d]]$w1
      w1[[d]]$w1 <- w0[[d]]$w1 <- rep(1/n,n)
    }
    
    wass_nnm_0 <- 
      # c(
      sapply(1:D_plus, 
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
      sapply(1:D_plus, function(d) 
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
      sapply(1:D_plus, function(d) 
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
      sapply(1:D_plus, function(d)
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
    keep <- scales != 0
    
    #for control
    nnm.adjusted0 <- (wass_nnm_0[1:D]^p * 1/scales)^(1/p)
    full.adjusted0 <- (wass_full_0[1:D]^p * 1/scales)^(1/p)
    min.grid0 <- pmax((max(nnm.adjusted0[keep])^p * scales)^(1/p), 1e-4)
    max.grid0 <- (max(full.adjusted0[keep])^p * scales)^(1/p)
    
    # min.grid0 <- wass_nnm_0[1:D]
    # max.grid0 <- wass_full_0[1:D]
    
    #for treated
    nnm.adjusted1 <- (wass_nnm_1[1:D]^p * 1/scales)^(1/p)
    full.adjusted1 <- (wass_full_1[1:D]^p * 1/scales)^(1/p)
    min.grid1 <- pmax((max(nnm.adjusted1[keep])^p * scales[keep])^(1/p), 1e-4)
    max.grid1 <- (max(full.adjusted1[keep])^p * scales[keep])^(1/p)

    # min.grid1 <- wass_nnm_1[1:D]
    # max.grid1 <- wass_full_1[1:D]
    
    if (add.joint) {
      # jl.grid0 <- sum(min.grid0^p)^(1/p)
      # ju.grid0 <- sum(max.grid0^p)^(1/p)
      # jl.grid1 <- sum(min.grid1^p)^(1/p)
      # ju.grid1 <- sum(max.grid1^p)^(1/p)
      # if (ju.grid0 < wass_full_0[[D + 1]]) 
        ju.grid0 <- wass_full_0[[D + 1]]
      # if (jl.grid0 < wass_nnm_0[[D + 1]]) 
        jl.grid0 <- wass_nnm_0[[D + 1]]
      # if (ju.grid1 < wass_full_1[[D + 1]]) 
        ju.grid1 <- wass_full_1[[D + 1]]
      # if (jl.grid1 > wass_nnm_1[[D + 1]]) 
        jl.grid1 <- wass_nnm_1[[D + 1]]
      
      jgrid0 <- exp(seq(log(jl.grid0), log(ju.grid0), length.out = grid.length))
      jgrid1 <- exp(seq(log(jl.grid1), log(ju.grid1), length.out = grid.length))
      
      # grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      # grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) seq(x[1], x[2], length.out = grid.length)))
      grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) seq(x[1], x[2], length.out = grid.length)))
      
      grid_0 <- do.call("cbind", lapply(jgrid0, function(joint) rbind(grid_1, joint)))
      grid_1 <- do.call("cbind", lapply(jgrid1, function(joint) rbind(grid_1, joint)))
      
    } else {
      grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      
      
    }
    
    grid <- lapply(1:ncol(grid_0), function(i) list(grid_0[,i],
                                          grid_1[,i]))
    
    
    
  } else {
    nnm <- lapply(cost, function(cc) calc_weight(data, estimand = estimand, method = "NNM", 
                                                 cost = cc, p = p,
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
    keep <- scales != 0
    
    nnm.adjusted  <- (wass_nnm[1:D]^p / scales)^(1/p)
    full.adjusted <- (wass_full[1:D]^p / scales)^(1/p)
    min.grid      <- pmax((max(nnm.adjusted[keep])^p * scales)^(1/p), 1e-4)
    max.grid      <- (max(full.adjusted[keep])^p * scales)^(1/p)
    
    # min.grid <- wass_nnm[1:D]
    # max.grid <- wass_full[1:D]
    
    if (add.joint) {
      # jl.grid <- sum(min.grid^p)^(1/p)
      # ju.grid <- sum(max.grid^p)^(1/p)
      # if (jl.grid < wass_nnm[[D + 1]]) 
      jl.grid <- wass_nnm[[D + 1]]
      # if (ju.grid > wass_full[[D + 1]]) 
      ju.grid <- wass_full[[D + 1]]
      
      jgrid <- exp(seq(log(jl.grid), log(ju.grid), length.out = grid.length))

      grid_temp <- t(apply(cbind(min.grid, max.grid), 1, function(x) seq(x[1], x[2], length.out = grid.length)))
      
      grid_temp <- do.call("cbind", lapply(jgrid, function(joint) rbind(grid_temp, joint)))
      
    } else {
        grid_temp <- t(apply(cbind(min.grid, max.grid), 1, function(x) seq(x[1], x[2], length.out = grid.length)))
        
    }
    grid <- lapply(1:ncol(grid_temp), function(i) grid_temp[,i] )
  }
  return(grid)
}

joint.cwass.fun.grid <- function(data, cost, grid.length, p, estimand, wass.iter, ...) {
  
  if (wass.iter < 1000) wass.iter <- wass.iter * 50
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
    
    grid <- rbind(exp(seq(log(wass_nnm[[1]]), log(wass_full[[1]]), length.out = grid.length)),
                  exp(seq(log(wass_nnm[[2]]), log(wass_full[[2]]), length.out = grid.length)))
    grid <- lapply(1:grid.length, function(i) grid[,i])
    
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
    
    grid <- exp(seq(log(wass_nnm), log(wass_full), length.out = grid.length))
    
  }
  return(grid)
}

cwass.fun.grid <- function(x, z, 
                           grid.length,
                           p, data, cost, estimand, metric, wass.iter, add.margins, add.joint, ...) {
  # D <- ncol(x)
  # 
  # n <- nrow(x)
  # 
  # x0 <- x[z == 0, , drop = FALSE]
  # x1 <- x[z == 1, , drop = FALSE]
  # 
  # n0 <- nrow(x0)
  # n1 <- nrow(x1)
  
  
  grid <- if (add.margins) {
    marg.cwass.fun.grid(x, z, grid.length, p, data, cost, estimand, metric, wass.iter, add.joint, ...)
  } else if (add.joint & !add.margins) {
    joint.cwass.fun.grid(data, cost, grid.length, p, estimand, wass.iter, ...)
  }
  
  return(grid)
  
}

wass.fun.grid <- function(x, z, 
                          grid.length,
                          p, 
                          data, cost, estimand, metric, wass.iter, add.margins, ...) {
  D <- ncol(x)
  
  n <- nrow(x)
  
  x0 <- x[z == 0, , drop = FALSE]
  x1 <- x[z == 1, , drop = FALSE]
  
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  marg.grid <- NULL
  if (add.margins) {
    args <- list(x = x, z = z, 
                 grid.length = grid.length,
                 p = p, data = data, 
                 cost = cost, 
                 estimand = estimand, metric = metric, 
                 wass.iter = wass.iter, add.joint = FALSE, ...)
    args <- args[!duplicated(names(args))]
    argn <- lapply(names(args), as.name)
    names(argn)  <- names(args)
    f.call <- as.call(c(list(as.name("marg.cwass.fun.grid")), argn) )
    marg.grid <- eval(f.call, envir = args)
    grid <- lapply(c(0, exp(seq(log(1e-3), log(1e6), length.out = length(marg.grid) - 1))),
                   function(nn) nn)
    # grid <- lapply(exp(seq(log(1e-3), log(1e6), length.out = length(marg.grid) )),
    #                function(nn) nn)
    keep <- round(seq.int(1L,length(marg.grid), length.out = grid.length))
    marg.grid <- marg.grid[keep]
    grid <- grid[keep]
    if (estimand == "ATE") {
      grid <- unlist(lapply(marg.grid, function(m) 
        lapply(grid, function(g) list(c(m[[1]],g),
                                          c(m[[2]],g)))
               ), recursive = FALSE)
    } else {
      grid <- unlist(lapply(marg.grid, function(m) 
        lapply(grid, function(g) c(m,g))),
        recursive = FALSE)
    }
    
  } else {
    # grid <- lapply(c(0, exp(seq(log(1e-3), log(1e6), length.out = grid.length - 1))),
    #                function(nn) nn)
    grid <- lapply(exp(seq(log(1e-3), log(1e6), length.out = grid.length )),
                   function(nn) nn)
    if (estimand == "ATE") {
      grid <- lapply(grid, function(gg) list(gg, gg))
    }
  }
  
  return(grid)
  
}

wass_grid_default <- function(x, z, grid.length,
                              p, data, cost, estimand, method, metric, wass.iter, 
                              add.joint, add.margins, ...) {
  
    
  get_defaults <- switch(method,
                         "Constrained Wasserstein" = as.name("cwass.fun.grid"),
                         "Wasserstein" = as.name("wass.fun.grid"),
                         "Sliced Wasserstein" = NULL)
  args <- list(x = x, z = z, 
               grid.length = grid.length,
               p = p, data = data, 
               cost = cost, estimand = estimand, metric = metric, 
               wass.iter = wass.iter, 
               add.joint = add.joint,
               add.margins = add.margins,
               ...)
  args <- args[!duplicated(names(args))]
  n.args <- lapply(names(args), as.name)
  names(n.args) <- names(args)
  f.call <- as.call(c(list(get_defaults), n.args))
  
  return(eval(f.call, envir = args))
}

wass_grid <- function(rowCount, colCount, weight, cost, x0, x1, wass.method, wass.iter, ...) {
  if (all(is.na(weight$w1) | all(is.na(weight$w0)))) return(NA)
  n1 <- length(weight$w1)
  n0 <- length(weight$w0)
  # ord.idx <- rep(NA_integer_, n1 + n0)
  # ord.idx[tx_idx] <- 1:n1
  # ord.idx[-tx_idx] <- 1:n0
  # tx_boot <- ord.idx[bootIdx[bootIdx %in% tx_idx]]
  # cn_boot <- ord.idx[bootIdx[!(bootIdx %in% tx_idx)]]
  # n1 <- length(tx_boot)
  # n0 <- length(cn_boot)
  # n <- n0 + n1
  
  # tab.row <- tabulate(bootIdx.row, nbins = n0)
  # tab.col <- tabulate(bootIdx.col, nbins = n0)
  # row.idx <- which(tab.row != 0)
  # col.idx <- which(tab.col != 0)
  # weight$w0 <- renormalize(weight$w0[row.idx] * tab.row[row.idx])
  # weight$w1 <- renormalize(weight$w1[col.idx] * tab.col[col.idx])
  
  weight$w0 <- renormalize(weight$w0 * rowCount)
  weight$w1 <- renormalize(weight$w1 * colCount)
  
  if (all(weight$w0 == 0) | all(weight$w1 == 0)) return(NULL)
  
  p.temp <- weight$args$power
  if (is.null(p.temp)) {
    p <- list(...)$p
  } else {
    p <- p.temp
  }
  
  # if (estimand == "ATE") {
  #   w0 <- w1 <- weight
  #   w0$gamma <- w1$gamma <- NULL
  #   w1$w0 <- w1$w1
  #   w1$w1 <- w0$w1 <- rep(1/n, n)
  #   if (is.list(cost[[1]]) & is.list(cost[[2]])) {
  #     c0 <- cost[[1]][[length(cost[[1]])]][cn_boot, bootIdx]
  #     c1 <- cost[[2]][[length(cost[[2]])]][tx_boot, bootIdx]
  #   } else {
  #     c0 <- cost[[1]][cn_boot, bootIdx]
  #     c1 <- cost[[2]][tx_boot, bootIdx]
  #   }
  #   wass0 <- wass_dist_helper(a = w0, b = NULL,
  #                          cost = c0,
  #                          p = p, method = wass.method, niter = wass.iter, ...)
  #   wass1 <- wass_dist_helper(a = w1, b = NULL,
  #                          cost = c1,
  #                          p = p, method = wass.method, niter = wass.iter, ...)
  #   return(((n1 * wass0^p + wass1^p * n0)/n)^(1/p))
  # } else {
    # if(!is.null(weight$gamma)) weight$gamma <- weight$gamma[cn_boot, tx_boot]
    if ( is.list(cost) ) {
      cc <- cost[[length(cost)]]#[rowCount > 0, colCount > 0]
    } else {
      cc <- cost#[rowCount > 0, colCount > 0]
    }
    return(wass_dist_helper(a = weight, b = NULL,
                         cost = cc, X = x0, Y = x1,
                         p = p, method = wass.method, niter = wass.iter, ...))
  # }
}
