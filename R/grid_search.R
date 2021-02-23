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
    out <- eval(f.call, envir = args)
    # out <- tryCatch(eval(f.call, envir = args),
    #                 error = function(e) {
    #                   warning(e$message)
    #                   return(list(w0 = rep(NA_real_, n0),
    #                                                  w1 = rep(NA_real_, n1)))})
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
    output.weight <- weight.list[[min.idx.0[1]]]
    output.weight$w1 <- weight.list[[min.idx.1[1]]]$w1
    output.weight$args$constraint <- output.weight$args$standardized.mean.difference <- c(grid[min.idx.0[1]], grid[min.idx.1[1]])
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
                             K = 10, R = 10,
                             n.boot = 1000,
                             eval.method = c("cross.validation", "bootstrap"),
                             method = c("Wasserstein","Constrained Wasserstein", "SCM"),
                             sample_weight = NULL,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE,
                             add.margins = FALSE,
                             joint.mapping = FALSE,
                             verbose = FALSE,
                             neg.weights = FALSE,
                             cgd = FALSE,
                             ...) 
{
  
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  eval.method <- match.arg(eval.method, c("cross.validation", "bootstrap"))
  
  add.margins <- isTRUE(add.margins)
  add.joint  <- isTRUE(add.joint)
  joint.mapping <- isTRUE(joint.mapping)
  neg.weights <- isTRUE(neg.weights)
  cgd <- isTRUE(cgd)
  
  # set optional arguments
  dots <- list(...)
  solver <- dots[["solver"]]
  p      <- dots[["p"]]
  metric <- dots[["metric"]]
  penalty <- dots[["penalty"]]
  if (is.null(solver)) solver <- "gurobi"
  if (is.null(p)) p <- 2
  if (is.null(metric)) metric <- "mahalanobis"
  if (is.null(penalty)) penalty <- "L2"
  
  
  if (method == "SCM") {
    add.margins <- add.joint <- joint.mapping <- FALSE
  }
  
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
  
  cost <- wg_cost_setup(cost = dots$cost, p = p,
                        estimand = estimand, x = x, z = z,
                        metric = metric, method = method,
                        rkhs.args = dots$rkhs.args,
                        add.margins = add.margins)
  
  
  if (all(is.null(grid)) || all(is.na(grid)) ) {
    grid <- f.call.list(fun = "wass_grid_default", 
                list.args = list(x = x, z = z, 
                  grid.length = grid.length,
                  p = p, data = data, cost = cost, estimand = estimand,
                  method = method, metric = metric, wass.iter = wass.iter, 
                  add.joint = add.joint, add.margins = add.margins,
                  joint.mapping = joint.mapping, penalty = penalty,
                  ...))
  }
  
  # if grid values return error, give out the nearest neighbor estimate
  if (!all(is.na(grid)) && isTRUE(all.equal(grid[1], grid[length(grid)], check.attributes = FALSE)) ) {
    return(calc_weight_NNM(data = data, estimand = estimand,
                           transport.matrix = TRUE,
                           ...))
  }
  
  args <- list(data = data, 
               x0 = x0,
               x1 = x1,
               x = x,
               z = z,
               grid = grid,  
               n.boot = n.boot,
               K = K, 
               R = R,
               eval.method = eval.method,
               wass.method = wass.method,
               wass.iter = wass.iter,
               sample_weight = sample_weight,
               estimand = estimand, 
               method = method, 
               solver = solver, 
               metric = metric,
               p = p, 
               cost = cost, 
               add.joint = add.joint,
               add.margins = add.margins, 
               joint.mapping = joint.mapping,
               neg.weights = neg.weights,
               cgd = cgd, 
               verbose = verbose,
               ...)
  
  weight.list <- weight_est_fun(args)
  
  output.weight <- eval_weights(weight.list, args)
  
  return(output.weight)
  
  
  # if (estimand != "ATE") {
  #   # bootIdx.rows <- lapply(1:n.boot, function(ii) {sample.int(n0,n0, replace = TRUE)})
  #   # bootIdx.cols <- lapply(1:n.boot, function(ii) {sample.int(n1,n1, replace = TRUE)})
  #   if (wass.method == "networkflow" | wass.method == "exact" | wass.method == "hilbert") {
  #     rowCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n0), prob = sample_weight$a)), simplify = FALSE)
  #     colCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n1), prob = sample_weight$b)), simplify = FALSE)
  #   } else {
  #     rowCount <- replicate(n.boot, c(rmultinom(1, n0, prob = sample_weight$a)), simplify = FALSE)
  #     colCount <- replicate(n.boot, c(rmultinom(1, n1, prob = sample_weight$b)), simplify = FALSE)
  #   }
  #   
  #   if (is.null(dots$cost_a)) {
  #     cost_a <- cost_fun(rbind(x0,x0), z = c(rep(1,n0),
  #                                  rep(0,n0)),
  #              power = p,
  #              metric = metric, estimand = "ATT")
  #   } else {
  #     cost_a <- dots$cost_a
  #     if (is.list(cost_a)) cost_a <- cost_a[[1]]
  #   }
  #   
  #   if (is.null(dots$cost_b)) {
  #     cost_b <- cost_fun(rbind(x1,x1), z = c(rep(1, n1),
  #                                            rep(0, n1)),
  #                        power = p,
  #                        metric = metric, estimand = "ATT")
  #   } else {
  #     cost_b <- dots$cost_b
  #     if (is.list(cost_b)) cost_b <- cost_b[[1]]
  #   }
  #   
  #   tx_idx <- which(pd$z == 1)
  #   boot.args <- list(FUN = wass_grid, 
  #                     # bootIdx.row = bootIdx.rows, 
  #                     # bootIdx.col = bootIdx.cols,
  #                     rowCount = rowCount, 
  #                     colCount = colCount,
  #                     MoreArgs    = list(weight = weight.list[[1]], 
  #                         data = wass.dat, 
  #                         # tx_idx = tx_idx, 
  #                         cost = cost,
  #                         p = p,
  #                         metric = metric,
  #                     # estimand = estimand,
  #                         wass.method = wass.method, 
  #                         wass.iter = wass.iter,
  #                         add.joint = add.joint,
  #                         x0 = x0,
  #                         x1 = x1,
  #                         cost_a = cost_a,
  #                         cost_b = cost_b,
  #                         ...)
  #   )
  #   boot.args <- boot.args[!duplicated(names(boot.args))]
  #   boot.args$MoreArgs <- boot.args$MoreArgs[!duplicated(names(boot.args$MoreArgs))]
  #   boot.n <- lapply(names(boot.args), as.name)
  #   names(boot.n) <- names(boot.args)
  #   b.call <- as.call(c(list(quote(mapply)), boot.n))
  #   
  #   output <- rep(NA, length(grid))
  #   names(output) <- as.character(grid)
  #   
  #   if (verbose ) {
  #     message("\nCalculating out of sample balance:")
  #     # pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
  #     pbapply::pboptions(type = "timer", style = 3, char = "=")
  #   } else {
  #     pbapply::pboptions(type = "none")
  #   }
  #   output <- pbapply::pbsapply(weight.list, function(ww) {
  #     boot.args$MoreArgs$weight <- ww
  #     return(mean(eval(b.call, envir = boot.args)^p))
  #   })
  #   # for (g in seq_along(grid)) {
  #   #   boot.args$MoreArgs$weight <- weight.list[[g]]
  #   #   output[g] <- mean(eval(b.call, envir = boot.args))
  #   #   if (verbose) setTxtProgressBar(pb, g)
  #   # }
  #   pbapply::pboptions(type = "none")
  #   if (all(is.na(output))) stop("wass_grid_search: All grid values generated errors")
  #   
  #   min.idx <- which(output == min(output, na.rm = TRUE))[1]
  #   if (length(min.idx) > 1) warning("Multiple penalties had minimum wasserstein value! Try increasing bootstrap number with `n.boot` parameter")
  #   
  #   return(weight.list[[min.idx[1]]])
  # } else {
  #   # bootIdx.cols    <- lapply(1:n.boot, function(ii) {sample.int(n,n, replace = TRUE)})
  #   # bootIdx.rows.0  <- lapply(1:n.boot, function(ii) {sample.int(n0,n0, replace = TRUE)})
  #   # bootIdx.rows.1  <- lapply(1:n.boot, function(ii) {sample.int(n1,n1, replace = TRUE)})
  #   if (wass.method == "networkflow" | wass.method == "exact" | wass.method == "hilbert") {
  #     rowCount.0 <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n0), prob = sample_weight$a)),     simplify = FALSE)
  #     rowCount.1 <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n1), prob = sample_weight$b)),     simplify = FALSE)
  #     colCount   <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n), prob = sample_weight$total)), simplify = FALSE)
  #   } else {
  #     rowCount.0 <- replicate(n.boot, c(rmultinom(1, n0, prob = sample_weight$a)),     simplify = FALSE)
  #     rowCount.1 <- replicate(n.boot, c(rmultinom(1, n1, prob = sample_weight$b)),     simplify = FALSE)
  #     colCount   <- replicate(n.boot, c(rmultinom(1, n, prob = sample_weight$total)), simplify = FALSE)
  #   }
  #   
  #   
  #   
  #   if (is.null(dots$cost_a)) {
  #     cost_a <- list(cost_fun(rbind(x0,x0), z = c(rep(1,n0),
  #                                            rep(0,n0)),
  #                        power = p,
  #                        metric = metric, estimand = "ATT"),
  #                    cost_fun(rbind(x1,x1), z = c(rep(1,n1),
  #                                                 rep(0,n1)),
  #                             power = p,
  #                             metric = metric, estimand = "ATT")
  #     )
  #   } else {
  #     cost_a <- dots$cost_a
  #     stopifnot(length(cost_a) == 2)
  #   }
  #   
  #   if (is.null(dots$cost_b)) {
  #     cost_b <- list(cost_fun(rbind(x,x), z = c(rep(1,n),
  #                                               rep(0,n)),
  #                             power = p,
  #                             metric = metric, estimand = "ATT"))
  #   } else {
  #     cost_b <- dots$cost_b
  #     if (is.list(cost_b)) cost_b <- cost_b[[1]]
  #   }
  #   
  #   output_0        <- output_1 <- rep(NA, length(grid))
  #   names(output_0) <- names(output_1) <- as.character(grid)
  #   tx_idx          <- which(pd$z == 1)
  #   full.sample.weight <- rep(1/n, n)
  #   
  #   boot.args <- list(FUN = wass_grid,
  #                     # bootIdx.row  = bootIdx.rows.0,  
  #                     # bootIdx.col = bootIdx.cols,
  #                     rowCount  = rowCount.0,  
  #                     colCount = colCount,
  #                     MoreArgs = list(
  #                       weight = weight.list[[1]], 
  #                       data = wass.dat, 
  #                       # tx_idx = tx_idx, 
  #                       cost = cost[[1]],
  #                       p = p,
  #                       metric = metric,
  #                       x0 = x0,
  #                       x1 = x,
  #                       cost_a = cost_a[[1]],
  #                       cost_b = cost_b,
  #                       # estimand = "ATT",
  #                       wass.method = wass.method, wass.iter = wass.iter,
  #                       add.joint = add.joint,
  #                       ...)
  #   )
  #   boot.args <- boot.args[!duplicated(names(boot.args))]
  #   boot.args$MoreArgs <- boot.args$MoreArgs[!duplicated(names(boot.args$MoreArgs))]
  #   boot.n <- lapply(names(boot.args), as.name)
  #   names(boot.n) <- names(boot.args)
  #   b.call <- as.call(c(list(quote(mapply)), boot.n))
  #   if (verbose ) {
  #     message("\nCalculating out of sample balance")
  #     pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
  #     
  #   }
  #   
  #   boot.args1 <- boot.args0 <- boot.args
  #   boot.args0$MoreArgs$cost   <- cost[[1]]
  #   boot.args1$MoreArgs$cost   <- cost[[2]]
  #   
  #   # boot.args0$bootIdx.row     <- bootIdx.rows.0
  #   # boot.args1$bootIdx.row     <- bootIdx.rows.1
  #   boot.args0$rowCount     <- rowCount.0
  #   boot.args1$rowCount     <- rowCount.1
  #   boot.args1$x0           <- x1
  #   boot.args1$cost_a       <- cost_a[[2]]
  #   
  #   w0 <- w1 <- weight.list[[length(weight.list)]]
  #   w0$w1 <- w1$w1 <- full.sample.weight
  #   w1$w0 <- weight.list[[length(weight.list)]]$w1
  #   
  #   w0$args$power <- w1$args$power <- weight.list[[length(weight.list)]]$args$power
  #   boot.args0$MoreArgs$weight <- w0
  #   boot.args1$MoreArgs$weight <- w1
  #   
  #   for (g in seq_along(grid)) {
  # 
  #     boot.args0$MoreArgs$weight$w0 <- weight.list[[g]]$w0
  #     # boot.args$MoreArgs$cost   <- cost[[1]]
  #     output_0[g]               <- mean(eval(b.call, envir = boot.args0)^p)
  #     
  #    
  #     boot.args1$MoreArgs$weight$w0 <- weight.list[[g]]$w1
  #     # boot.args$MoreArgs$cost   <- cost[[2]]
  #     output_1[g]               <- mean(eval(b.call, envir = boot.args1)^p)
  #     if (verbose) setTxtProgressBar(pb, g)
  #   }
  #   if (verbose) close(pb)
  #   if (all(is.na(output_0)) | all(is.na(output_1))) stop("wass_grid_search: All grid values generated errors")
  #   
    min.idx.0 <- which(output_0 == min(output_0, na.rm = TRUE))
    min.idx.1 <- which(output_1 == min(output_1, na.rm = TRUE))
    
    if(length(min.idx.0) > 1) warning("Multiple penalties for control had minimum wasserstein value! Try increasing bootstrap number with `n.boot` parameter")
    if(length(min.idx.1) > 1) warning("Multiple penalties for treated had minimum wasserstein value! Try increasing bootstrap number with `n.boot` parameter")
    
    output.weight <- weight.list[[min.idx.0[1]]]
    output.weight$w1 <- weight.list[[min.idx.1[1]]]$w1
    output.weight$args$constraint <- list(weight.list[[min.idx.0[1]]]$args$constraint[1],
                                       weight.list[[min.idx.1[1]]]$args$constraint[2])
    # weight.list[[min.idx]]$args$standardized.mean.difference <- grid[min.idx]
    return(output.weight)
  # }
  
}

wg_cost_setup <- function(cost, p, estimand, x, z,
                          metric, method, rkhs.args,
                          add.margins) {
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
  if (is.null(cost) || rerun) {
    # if(is.null(dots[["p"]])) dots[["p"]] <- 2
    # if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    
    # if(estimand == "ATE" & metric != "RKHS") {
    #   cost <- list(cost_fun(x0, x, ground_p = p, metric = metric),
    #                cost_fun(x1, x, ground_p = p, metric = metric))
    #     
    # } else {
    
    # }
    cost <- wass_cost_default(x = x, z = z, estimand = estimand, 
                              metric = metric, method = method, p = p, 
                              rkhs.args = rkhs.args, 
                              add.margins = add.margins)
  } 
  return(cost)
}

weight_est_fun <- function(args) {
  # args = list(data = data, grid = grid,  
  # estimand = estimand, 
  # method = method, solver = solver, 
  # metric = metric,
  # p = p, cost = cost, 
  # add.joint = add.joint,
  # add.margins = add.margins, 
  # joint.mapping = joint.mapping,
  # neg.weights = neg.weights,
  # cgd = cgd, verbose = verbose,
  # ...)
  
  if (args$verbose) {
    message("\nEstimating Optimal Transport weights for each constraint")
    pbapply::pboptions(type = "timer", style = 3, char = "=")
  } else {
    pbapply::pboptions(type = "none")
  }
  
  # if (args[["estimand"]] == "ATE") {
  #   args0 <- args1 <- args
  #   
  #   
  # } else {
    n0 <- nrow(args[["x0"]])
    n1 <- nrow(args[["x1"]])
    
    
    
    if (!args[["cgd"]]) {
      args[["x"]] <- args$z <- args$y <- args[["x1"]] <- args[["x0"]] <- NULL
      args$constraint <- args[["grid"]][1]
      weight.est.call <- f.call.list.no.eval("calc_weight",
                                             list.args = args)
      weights <-  pbapply::pblapply(args[["grid"]], function(delta) {
        weight.est.call$envir$constraint <- delta
        out <- eval(weight.est.call$expr, envir = weight.est.call$envir)
        # out <- tryCatch(eval(weight.est.call$expr, envir = weight.est.call$envir),
        #                 error = function(e) {
        #                   warning(e$message)
        #                   return(list(w0 = rep(NA_real_, n0),
        #                               w1 = rep(NA_real_, n1),
        #                               gamma = NULL,
        #                               estimand = args[["estimand"]],
        #                               method = args[["method"]], args = list(constraint = delta)))})
        if (args[["solver"]] == "gurobi") Sys.sleep(0.1)
        class(out) <- "causalWeights"
        return(out)
      })
    } else {
      weights <- cgd_grid_fun(args)
    }
    
  # }
  
  pbapply::pboptions(type = "none")
  
  return(weights)
}

cgd_grid_fun <- function(args) {
  # args = list(data = data, grid = grid,  
  # estimand = estimand, 
  # method = method, solver = solver, 
  # metric = metric,
  # p = p, cost = cost, 
  # add.joint = add.joint,
  # add.margins = add.margins, 
  # joint.mapping = joint.mapping,
  # neg.weights = neg.weights,
  # cgd = cgd, verbose = verbose,
  # ...)
  x0 <- args[["x0"]]
  x1 <- args[["x1"]]
  z  <- args[["z"]]
  
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  pen <- args$penalty
  if (length(pen) == 0) pen <- "none"
  qp_pen <- switch(pen,
                   "L2" = "L2",
                   "variance" = "variance",
                   "entropy" = "L2",
                   "none")
  
  if (args[["estimand"]] == "ATE" || args[["estimand"]] == "cATE") {
    x1_0 <- x0
    x1_1 <- x1
    x2 <- rbind(x0,x1)
    z_0 <- c(rep(0, nrow(x0)), rep(1,nrow(x)))
    z_1 <- c(rep(0, nrow(x1)), rep(1,nrow(x)))
  } else if (args[["estimand"]] == "ATT") {
    x1 <- x0
    x2 <- x1
    z <- z
  } else if (args[["estimand"]] == "ATT") {
    x1 <- x1
    x2 <- x0
    z <- 1 - z
  }
  
  if (args[["estimand"]] == "ATE" || args[["estimand"]] == "cATE") {
    optimizer_0 <- wassCGD$new(X1 = x1_0,
                             X2 = x2,
                             z = z_0,
                             cost = args[["cost"]][[1]],
                             qp_constraint = list(joint = args[["grid"]][[1]]$joint,
                                                  margins = args[["grid"]][[1]]$margins,
                                                  penalty = args[["grid"]][[1]]$penalty),
                             qp_solver = args[["solver"]],
                             qp_penalty = qp_pen,
                             lambda = list(joint = args[["grid"]][[1]]$joint,
                                           penalty = args[["grid"]][[1]]$penalty
                             ),
                             add.mapping = args[["joint.mapping"]],
                             add.margins = args[["add.margins"]],
                             penalty = pen,
                             metric = args[["metric"]],
                             power = args[["p"]],
                             niter = args[["niter"]],
                             tol = args[["tol"]],
                             sample_weight = args[["sample_weight"]])
    
    optimizer_1 <- wassCGD$new(X1 = x1_1,
                               X2 = x2,
                               z = z_1,
                               cost = args[["cost"]][[2]],
                               qp_constraint = list(joint = args[["grid"]][[2]]$joint,
                                                    margins = args[["grid"]][[2]]$margins,
                                                    penalty = args[["grid"]][[2]]$penalty),
                               qp_solver = args[["solver"]],
                               qp_penalty = qp_pen,
                               lambda = list(joint = args[["grid"]][[2]]$joint,
                                             penalty = args[["grid"]][[2]]$penalty
                               ),
                               add.mapping = args[["joint.mapping"]],
                               add.margins = args[["add.margins"]],
                               penalty = pen,
                               metric = args[["metric"]],
                               power = args[["p"]],
                               niter = args[["niter"]],
                               tol = args[["tol"]],
                               sample_weight = args[["sample_weight"]])
    
    weight0 <- pbapply::pblapply(args[["grid"]], function(delta) {
      optimizer_0$cold_start(list(joint = delta$joint, parameters =  0,
                                penalty = delta$penalty,
                                margins = delta$margins))
      out <- tryCatch(cg(optimizer_0, verbose = FALSE)$return_cw(),
                      error = function(e) {
                        warning(e$message)
                        return(list(w0 = rep(NA_real_, n0),
                                    w1 = rep(NA_real_, n1),
                                    gamma = NULL,
                                    estimand = estimand,
                                    method = method, args = list(constraint = delta)))})
      # out <- do.call("calc_weight_bal", args)
      if (args[["solver"]] == "gurobi") Sys.sleep(0.1)
      class(out) <- "causalWeights"
      return(out)
    })
    weight1 <- pbapply::pblapply(args[["grid"]], function(delta) {
      optimizer_0$cold_start(list(joint = delta$joint, parameters =  0,
                                  penalty = delta$penalty,
                                  margins = delta$margins))
      out <- tryCatch(cg(optimizer_0, verbose = FALSE)$return_cw(),
                      error = function(e) {
                        warning(e$message)
                        return(list(w0 = rep(NA_real_, n0),
                                    w1 = rep(NA_real_, n1),
                                    gamma = NULL,
                                    estimand = estimand,
                                    method = method, args = list(constraint = delta)))})
      # out <- do.call("calc_weight_bal", args)
      if (args[["solver"]] == "gurobi") Sys.sleep(0.1)
      class(out) <- "causalWeights"
      return(out)
    })
    
    weight <- mapply(FUN = combine_weight_ATE, weight0 = weight0, weight1 = weight1, SIMPLIFY = FALSE)
    
  } else {
    optimizer <- wassCGD$new(X1 = x1,
                             X2 = x2,
                             z = z,
                             cost = args[["cost"]],
                             qp_constraint = list(joint = args[["grid"]]$joint,
                                                  margins = args[["grid"]]$margins,
                                                  penalty = args[["grid"]]$penalty),
                             qp_solver = args[["solver"]],
                             qp_penalty = qp_pen,
                             lambda = list(joint = args[["grid"]]$joint,
                                           penalty = args[["grid"]]$penalty
                             ),
                             add.mapping = args[["joint.mapping"]],
                             add.margins = args[["add.margins"]],
                             penalty = pen,
                             metric = args[["metric"]],
                             power = args[["p"]],
                             niter = args[["niter"]],
                             tol = args[["tol"]],
                             sample_weight = args[["sample_weight"]])
    
    weight <- pbapply::pblapply(args[["grid"]], function(delta) {
      optimizer$cold_start(list(joint = delta$joint, parameters =  0,
                                penalty = delta$penalty,
                                margins = delta$margins))
      out <- tryCatch(cg(optimizer, verbose = FALSE)$return_cw(),
                      error = function(e) {
                        warning(e$message)
                        return(list(w0 = rep(NA_real_, n0),
                                    w1 = rep(NA_real_, n1),
                                    gamma = NULL,
                                    estimand = estimand,
                                    method = method, args = list(constraint = delta)))})
      # out <- do.call("calc_weight_bal", args)
      if (args[["solver"]] == "gurobi") Sys.sleep(0.1)
      class(out) <- "causalWeights"
      return(out)
    })
    weight.out <- lapply(weight, function(w) {
      if (all(is.na(w$w1))) {
        return(w)
      } else {
        convert_sol(sol = list(c(w$gamma)), 
                    estimand = args[["estimand"]],
                    method = args[["method"]], 
                    n0 = n0, n1 = n1, sample_weight = args[["sample_weight"]])
      }
    }
    )
    for (i in seq_along(weight)) {
      weight.out[[i]]$args <- weight[[i]]$args
      weight.out[[i]]$estimand <- args[["estimand"]]
      weight.out[[i]]$method <- args[["method"]]
    }
  }
  
  return(weight.out)
}

eval_weights <- function(weights, args) {
  eval.method <- match.arg(args[["eval.method"]], c("cross.validation", "bootstrap") )
  estimand <- args[["estimand"]]
  
  if (estimand == "ATE") {
    
    weights.sep <- separate_weights_ATE(weights)
    
    args0 <- c(list(weights = weights.sep$w0), args)
    args1 <- c(list(weights = weights.sep$w1), args)
    
    args0$cost <- args$cost[[1]]
    args1$cost <- args$cost[[2]]
    
    args1$sample_weight$a  <- args1$sample_weight$b
    
    args0$sample_weight$b      <- args1$sample_weight$b <- args1$sample_weight$total
    args0$sample_weight$total  <- renormalize(c(args0$sample_weight$a, args0$sample_weight$b))
    args1$sample_weight$total  <- renormalize(c(args1$sample_weight$a, args1$sample_weight$b))
    
    args1$x0 <- args1$x1
    args0$x1 <- args1$x1 <- args1$x
    
    if (eval.method == "cross.validation") {
      wp0 <- f.call.list(fun = "wass_cv", list.args = args0)
      wp1 <- f.call.list(fun = "wass_cv", list.args = args1)
    } else if (eval.method == "bootstrap") {
      wp0 <- f.call.list(fun = "wass_boot", list.args = args0)
      wp1 <- f.call.list(fun = "wass_boot", list.args = args1)
    }
    
    sel0 <- which(wp0 == min(wp0, na.rm = TRUE))
    sel1 <- which(wp1 == min(wp1, na.rm = TRUE))
    
    sel.weights0 <- clean_up_weights(weights.sep$w0, sel0, args)
    sel.weights1 <- clean_up_weights(weights.sep$w1, sel1, args)
    
    weight <- combine_weight_ATE(sel.weights0, sel.weights1)
    
  } else  {
    eval.args <- args
    if (args[["estimand"]] == "ATC") {
      eval.args$estimand <- "ATT"
      
      weight.eval <- weight_list_ATT_to_ATC(weights)
      
      eval.args$cost <- cost_list_ATT_to_ATC(eval.args$cost)
      
      names(eval.args$sample_weight)[1:2] <- c("b","a")
      eval.args$z <- 1 - eval.args$z
      
      temp.argn <- names(eval.args)
      sel.args  <- grep("x1|x0", temp.argn)
      switch.names <- temp.argn[sel.args]
      names(eval.args)[sel.args] <- rev(switch.names)
    } else {
      weight.eval <- weights
    }
    if (eval.method == "cross.validation") {
      cv.args <- c(list(weights = weight.eval), eval.args)
      wp <- f.call.list(fun = "wass_cv", list.args = cv.args)
    } else if (eval.method == "bootstrap") {
      boot.args <- c(list(weights = weight.eval), eval.args)
      wp <- f.call.list(fun = "wass_boot", list.args = boot.args)
    }
    
    sel <- which(wp == min(wp, na.rm = TRUE))
    weight <- clean_up_weights(weights, sel, args)
  }
  return(weight)
}

weight_list_ATT_to_ATC <- function(weights) {
  weights.new <- weights
  for (i in seq_along(weights)) {
    
    weights.new[[i]]$estimand <- "ATT"
    weights.new[[i]]$w0 <- weights[[i]]$w1
    weights.new[[i]]$w1 <- weights[[i]]$w0
    weights.new[[i]]$gamma <- t(weights.new[[i]]$gamma)
  }
  
  return(weights.new)
}

cost_list_ATT_to_ATC <- function(cost) {
  
  if(is.list(cost)) {
    return(lapply(cost, t))
  } else {
    return(t(cost))
  }
}

separate_weights_ATE <- function(weights) {
  weights0 <- weights1 <- weights
  
  # separate weights
  for (i in seq_along(weights)) {
    # if (is.na(weights1[[i]]$w0) || is.na(weights1[[i]]$w1)) {
    #   next
    # }
    weights1[[i]]$w0 <- weights1[[i]]$w1
    weights0[[i]]$gamma <- weights0[[i]]$gamma[[1]]
    weights1[[i]]$gamma <- weights1[[i]]$gamma[[2]]
    
    weights0[[i]]$args$constraint <- weights0[[i]]$args$constraint[1]
    weights1[[i]]$args$constraint <- weights1[[i]]$args$constraint[2]
    
    weights0[[i]]$estimand <- "ATT"
    weights1[[i]]$estimand <- "ATT"
    
    weights0[[i]]$w1 <- colSums(weights0[[i]]$gamma)
    weights1[[i]]$w1 <- colSums(weights1[[i]]$gamma)
  }
  
  return(list(w0 = weights0, w1 = weights1))
}

combine_weight_ATE <- function(weight0, weight1) {
    weight0$w1 <- weight1$w0
    weight0$gamma <- list(weight0$gamma,
                          weight1$gamma)
    weight0$args$constraint <- list(weight0$args$constraint,
                                    weight1$args$constraint)
    weight0$estimand <- "ATE"
  return(weight0)
}

clean_up_weights <- function(weights, selection, args) {
  if (isTRUE(args[["cgd"]]) ) {
    check.weights <- weights[selection]
    grid <- lapply(check.weights, function(w) w$args$constraint)
    args[["cgd"]] <- FALSE
    args[["grid"]] <- grid
    # args[["estimand"]] <-"ATT"
    new.weights <- weight_est_fun(args)
    sel.weight <- eval_weights(new.weights, args)
    return(sel.weight)
  } else {
    return(weights[[selection[[1]]]])
  }
}

mean_bal_grid <- function(bootIdx, weight, data, tx_ind, ...) {
  if(all(is.na(weight$w0)) || all(is.na(weight$w1))) return(NA_real_)
  wvec <- c(weight$w0, weight$w1 )
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

marg.cwass.fun.grid <- function(x, z, grid.length, p, data, cost, estimand, metric, wass.iter, add.joint, 
                                remove.joint = FALSE, ...) {
  if (estimand == "ATE") {
    D <- length(cost[[1]])
    D_plus <- D
    if (add.joint) {
      D <- D - 1
    } else if (remove.joint) {
      D_plus <- D <- D - 1
    }
  } else {
    D <- length(cost)
    D_plus <- D
    if (add.joint) {
      D <- D - 1
    } else if (remove.joint) {
      D_plus <- D <- D - 1
    }
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
                                                    p = p,
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
    adj.factor <- 1.25
    
    #for control
    nnm.adjusted0 <- (wass_nnm_0[1:D]^p * 1/scales)^(1/p)
    full.adjusted0 <- (wass_full_0[1:D]^p * 1/scales)^(1/p)
    min.grid0 <- pmax((max(nnm.adjusted0[keep])^p * scales)^(1/p), 1e-4)
    max.grid0 <- (max(full.adjusted0[keep])^p * scales)^(1/p) * adj.factor
    
    # min.grid0 <- wass_nnm_0[1:D]
    # max.grid0 <- wass_full_0[1:D]
    
    #for treated
    nnm.adjusted1 <- (wass_nnm_1[1:D]^p * 1/scales)^(1/p)
    full.adjusted1 <- (wass_full_1[1:D]^p * 1/scales)^(1/p)
    min.grid1 <- pmax((max(nnm.adjusted1[keep])^p * scales[keep])^(1/p), 1e-4)
    max.grid1 <- (max(full.adjusted1[keep])^p * scales[keep])^(1/p) * adj.factor

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
      
      grid_0 <- do.call("cbind", lapply(jgrid0, function(joint) rbind(grid_0, joint)))
      grid_1 <- do.call("cbind", lapply(jgrid1, function(joint) rbind(grid_1, joint)))
      
        
      
      grid_0 <- lapply(1:ncol(grid_0), function(g) list(margins = grid_0[1:D,g], joint = grid_0[D_plus,g]))
      grid_1 <- lapply(1:ncol(grid_1), function(g) list(margins = grid_1[1:D,g], joint = grid_1[D_plus,g]))
      
    } else {
      grid_0 <- t(apply(cbind(min.grid0, max.grid0), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      grid_1 <- t(apply(cbind(min.grid1, max.grid1), 1, function(x) exp(seq(log(x[1]), log(x[2]), length.out = grid.length))))
      
      grid_0 <- lapply(1:ncol(grid_0), function(g) list(margins = grid_0[,g]))
      grid_1 <- lapply(1:ncol(grid_1), function(g) list(margins = grid_1[,g]))
    }
    
    grid <- lapply(1:length(grid_0), function(i) list(grid_0[[i]],
                                          grid_1[[i]]))
    
    
    
  } else {
    nnm <- lapply(cost, function(cc) calc_weight(data, estimand = estimand, method = "NNM", 
                                                 cost = cc, p = p,
                                                 ...))
    wass_nnm <- sapply(1:D_plus, function(d) wass_dist_helper(a = nnm[[d]], cost = cost[[d]], p = p, method = "networkflow", niter = wass.iter, ...) )
    wass_full <- 
      # c(
      sapply(1:D_plus, function(d) 
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
    adj.factor <- 1.25
    
    nnm.adjusted  <- (wass_nnm[1:D]^p / scales)^(1/p)
    full.adjusted <- (wass_full[1:D]^p / scales)^(1/p)
    min.grid      <- pmax((max(nnm.adjusted[keep])^p * scales)^(1/p), 1e-4)
    max.grid      <- (max(full.adjusted[keep])^p * scales)^(1/p) * adj.factor
    
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
      
      grid <- lapply(1:ncol(grid_temp), function(i) list(margins = grid_temp[1:D,i],
                                                         joint = grid_temp[D_plus,i]))
      
    } else {
      grid_temp <- t(apply(cbind(min.grid, max.grid), 1, function(x) seq(x[1], x[2], length.out = grid.length)))
      grid <- lapply(1:ncol(grid_temp), function(i) list(margins = grid_temp[,i] ))
        
    }
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
    grid <- lapply(1:grid.length, function(i) list(list(joint = grid[1,i]),
                                                   list(joint = grid[2,i])))
    
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
    grid <- lapply(grid, function(g) list(joint = g))
  }
  return(grid)
}

cwass.fun.grid <- function(x, z, 
                           grid.length,
                           p, data, cost, estimand, metric, wass.iter, add.margins, add.joint, 
                           joint.mapping, penalty, ...) {
  
  
  
  grid <- if (add.margins) {
    marg.cwass.fun.grid(x, z, grid.length, p, data, cost, estimand, metric, wass.iter, add.joint, ...)
  } else if (add.joint & !add.margins) {
    joint.cwass.fun.grid(data, cost, grid.length, p, estimand, wass.iter, ...)
  }
  
  if (joint.mapping) {
    
    pen.grid <- pen.fun.grid(x, z, grid.length, p, data, 
                             cost, estimand, metric, 
                             wass.iter, add.margins, 
                             joint.mapping, penalty, ...)
    
    if (estimand == "ATE") {
      for (i in length(grid)) {
        grid[[i]][[1]] <- c(grid[[i]][[1]], penalty = pen.grid[[i]][[1]])
        grid[[i]][[2]] <- c(grid[[i]][[2]], penalty = pen.grid[[i]][[2]])
        
      }
    } else {
      for (i in length(grid)) {
        grid[[i]] <- c(grid[[i]], penalty = pen.grid[[i]])
      }
    }
    
  }
  
  return(grid)
  
}

pen.fun.grid <- function(x, z, 
                          grid.length,
                          p, 
                          data, cost, estimand, metric, wass.iter, add.margins, 
                          joint.mapping, penalty, ...) {
  
  
  n <- nrow(x)
  
  # x0 <- x[z == 0, , drop = FALSE]
  # x1 <- x[z == 1, , drop = FALSE]
  # 
  # n0 <- nrow(x0)
  # n1 <- nrow(x1)
  if (penalty == "none") {
    if (estimand == "ATE") {
      return(list(list(list(penalty = 0.0), list(penalty = 0.0))))
    } else {
      return(list(list(penalty = 0.0)))
    }
  } 
  
  if (estimand == "ATE") {
    cost1 <- cost[[1]]
    cost0 <- cost[[2]]
    D <- length(cost1)
    
    if (add.margins) {
      cost1 <- cost1[[D]]
      cost0 <- cost0[[D]]
    }
    mc1 <- median(cost1)
    mc0 <- median(cost0)
    
    if (joint.mapping) {
      mc1 <- mc1/median(cost1)
      mc0 <- mc0/median(cost0)
    }
    
    grid0 <- lapply(c(0, exp(seq(log(1e-6 * mc0), log(1e3 * mc0), length.out = grid.length ))),
                   function(nn) nn)
    
    grid1 <- lapply(c(0, exp(seq(log(1e-6 * mc1) , log(1e3 * mc1), length.out = grid.length ))),
                   function(nn) nn)
    
    grid <- lapply(1:grid.length, function(gg) list(list(penalty = grid0[[gg]]), 
                                                    list(penalty = grid1[[gg]])))
  } else {
    
    D <- length(cost)
    
    if (add.margins) {
      cost <- cost[[D]]
    }
    mc <- median(cost)
    
    if (joint.mapping) {
      mc <- mc/median(cost)
    }
    grid <- lapply(c(0,exp(seq(log(1e-6 * mc), 
                               log(1e3 * mc) , 
                               length.out = grid.length ))),
                   function(nn) list(penalty = nn))
    
  }
  
  return(grid)
  
}

wass.fun.grid <- function(x, z, 
                          grid.length,
                          p, 
                          data, cost, estimand, metric, wass.iter, add.margins, 
                          joint.mapping, 
                          penalty, ...) {
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
                 wass.iter = wass.iter, add.joint = FALSE, 
                 remove.joint = TRUE, ...)
    args <- args[!duplicated(names(args))]
    argn <- lapply(names(args), as.name)
    names(argn)  <- names(args)
    f.call <- as.call(c(list(as.name("marg.cwass.fun.grid")), argn) )
    marg.grid <- eval(f.call, envir = args)
    grid <- pen.fun.grid(x, z, 
                        grid.length,
                        p, 
                        data, cost, estimand, 
                        metric, wass.iter, add.margins, 
                        joint.mapping, penalty,
                        ...)
    # keep <- round(seq.int(1L,length(marg.grid), length.out = grid.length))
    # marg.grid <- marg.grid[keep]
    # grid <- grid[keep]
    if (estimand == "ATE") {
      d_c <- length(cost[[1]]) - 1
      grid <- unlist(lapply(marg.grid, function(m) 
        lapply(grid, function(g) list(c(m[[1]], g[[1]]),
                                      c(m[[2]], g[[2]])))
               ), recursive = FALSE)
    } else {
      d_c <- length(cost) - 1
      grid <- unlist(lapply(marg.grid, function(m) 
        lapply(grid, function(g) c(list(margins = m$margins[1:d_c]), 
                                   g))),
        recursive = FALSE)
    }
    
  } else {
    # grid <- lapply(c(0, exp(seq(log(1e-3), log(1e6), length.out = grid.length - 1))),
    #                function(nn) nn)
    grid <- pen.fun.grid(x, z, 
                         grid.length,
                         p, 
                         data, cost, estimand, metric, 
                         wass.iter, add.margins, 
                         joint.mapping, penalty, ...)
    # if (estimand == "ATE") {
    #   grid <- lapply(grid, function(gg) list(gg, gg))
    # }
  }
  
  if (joint.mapping) {
    jm <- rep(exp(seq(log(1e-3), log(10), length.out = grid.length)), each = length(grid))
    
    grid <- rep(grid, grid.length)
    
    if (estimand == "ATE") {
      for (i in 1:length(grid)) {
        grid[[i]][[1]] <- c(grid[[i]][[1]], joint = list(jm[i]))
        grid[[i]][[2]] <- c(grid[[i]][[2]], joint = list(jm[i]))
        
      }
    } else {
      for (i in 1:length(grid)) {
        grid[[i]] <- c(grid[[i]], joint = list(jm[i]))
      }
    }
  }
  
  return(grid)
  
}

scm.fun.grid <- function(x, z, 
                          grid.length,
                          p, 
                          data, cost, estimand, metric, wass.iter, add.margins, 
                          joint.mapping, penalty, ...) {
  D <- ncol(x)
  
  n <- nrow(x)
  
  x0 <- x[z == 0, , drop = FALSE]
  x1 <- x[z == 1, , drop = FALSE]
  
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  cost <- cost_fun(x, z, power = 2, metric = "Lp", estimand = estimand)
  
  grid <- pen.fun.grid(x, z, 
                       grid.length,
                       p = 2, 
                       data, cost = cost, estimand, metric = "Lp", 
                       wass.iter, add.margins = FALSE, 
                       joint.mapping = FALSE, penalty, ...)
    
  return(grid)
  
}

wass_grid_default <- function(x, z, grid.length,
                              p, data, cost, estimand, method, metric, wass.iter, 
                              add.joint, add.margins, joint.mapping, 
                              penalty, ...) {
  
    
  get_defaults <- switch(method,
                         "Constrained Wasserstein" = as.name("cwass.fun.grid"),
                         "Wasserstein" = as.name("wass.fun.grid"),
                         "SCM" = as.name("scm.fun.grid"),
                         "Sliced Wasserstein" = NULL)
  args <- list(x = x, z = z, 
               grid.length = grid.length,
               p = p, data = data, 
               cost = cost, estimand = estimand, metric = metric, 
               wass.iter = wass.iter, 
               add.joint = add.joint,
               add.margins = add.margins,
               joint.mapping = joint.mapping,
               penalty = penalty,
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

get_boot <- function(n.boot, n0, n1, sample_weight, wass.method, ...) {
  
  # bootIdx.rows <- lapply(1:n.boot, function(ii) {sample.int(n0,n0, replace = TRUE)})
  # bootIdx.cols <- lapply(1:n.boot, function(ii) {sample.int(n1,n1, replace = TRUE)})
  if (wass.method == "networkflow" | wass.method == "exact" | wass.method == "hilbert") {
    rowCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n0), prob = sample_weight$a)), simplify = FALSE)
    colCount <- replicate(n.boot, c(rmultinom(1, adjust_m_of_n_btstrp(n1), prob = sample_weight$b)), simplify = FALSE)
  } else {
    rowCount <- replicate(n.boot, c(rmultinom(1, n0, prob = sample_weight$a)), simplify = FALSE)
    colCount <- replicate(n.boot, c(rmultinom(1, n1, prob = sample_weight$b)), simplify = FALSE)
  }
  
  return(list(rowCount = rowCount, colCount = colCount))
}

setup_boot_args <- function(boot.idx, weight.list, wass.dat, cost, p,
                      metric, wass.method, wass.iter,add.joint,
                      x0, x1, cost_a, cost_b, unbiased, verbose, ...){
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  rowCount <- boot.idx$rowCount
  colCount <- boot.idx$colCount
  
  if ( is.list(cost) ) {
    cc <- cost[[length(cost)]]#[rowCount > 0, colCount > 0]
  } else {
    cc <- cost#[rowCount > 0, colCount > 0]
  }
  
  entropy.meth.sel <- wass.method %in% c("sinkhorn","greenkhorn")
  
  if (is.null(cost_a) && isTRUE(unbiased) && entropy.meth.sel) {
    cost_a <- cost_fun(rbind(x0,x0), z = c(rep(1,n0),
                                           rep(0,n0)),
                       power = p,
                       metric = metric, estimand = "ATT")
  } else {
    if (is.list(cost_a)) cost_a <- cost_a[[1]]
  }
  
  if (is.null(cost_b) && isTRUE(unbiased) && entropy.meth.sel) {
    cost_b <- cost_fun(rbind(x1,x1), z = c(rep(1, n1),
                                           rep(0, n1)),
                       power = p,
                       metric = metric, estimand = "ATT")
  } else {
    if (is.list(cost_b)) cost_b <- cost_b[[1]]
  }
  
  boot.args <- list(FUN = "wass_grid", 
                    # bootIdx.row = bootIdx.rows, 
                    # bootIdx.col = bootIdx.cols,
                    rowCount = rowCount, 
                    colCount = colCount,
                    MoreArgs    = list(weight = weight.list[[1]], 
                                       # data = wass.dat, 
                                       # tx_idx = tx_idx, 
                                       cost = cc,
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
                                       unbiased = unbiased,
                                       ...)
  )
  boot.args$MoreArgs <- boot.args$MoreArgs[!duplicated(names(boot.args$MoreArgs))]
  
  boot.f.call <- f.call.list.no.eval(fun = "mapply",
                                     list.args = boot.args)
  
  return(boot.f.call)
}

wass_boot <- function(weights, n.boot, x0, x1, cost, p, metric = metric,
                      wass.method, wass.iter, add.joint,
                      sample_weight, cost_a = NULL, cost_b = NULL,
                      unbiased = FALSE,
                      verbose, ...) {
    
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    
    boot.idx <- get_boot(n.boot, n0, n1, sample_weight, wass.method, ...)
    
    boot.args <- setup_boot_args(boot.idx = boot.idx, 
                                 weight.list = weights, 
                                 wass.dat = wass.dat, 
                                 cost = cost, p,
                                 metric = metric, 
                                 wass.method = wass.method, 
                                 wass.iter = wass.iter, add.joint = add.joint,
                                 x0 = x0, x1 = x1, 
                                 cost_a = cost_a, cost_b = cost_b,
                                 unbiased = unbiased, ...)
    
    
    if (verbose ) {
      message("\nCalculating out of sample balance via bootstrap:")
      # pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
      pbapply::pboptions(type = "timer", style = 3, char = "=")
    } else {
      pbapply::pboptions(type = "none")
    }
    output <- pbapply::pbsapply(weights, function(ww) {
      boot.args$envir$MoreArgs$weight <- ww
      return(mean(eval(boot.args$expr, envir = boot.args$envir)^p))
    })
    
    pbapply::pboptions(type = "none")
    if (all(is.na(output))) stop("wass_grid_search: All grid values generated errors")
    
    return(output)
}

cv_sep <- function(fold, weight, cost, p, wass.method, wass.iter, cost_a, cost_b, 
                   unbiased, ...) {
  
  w0_raw    <- rowSums(weight$gamma[,-fold, drop = FALSE])
  w1_raw    <- colSums(weight$gamma[,fold, drop = FALSE])
  weight$w0 <- renormalize(w0_raw)
  weight$w1 <- renormalize(w1_raw)
  cost      <- cost[,fold, drop = FALSE]
  
  if (!is.null(cost_b)) {
    cost_b <- cost_b[fold,fold, drop = FALSE]
  }
  
  
  return(wass_dist_helper(a = weight, b = NULL,
                          cost = cost,
                          p = p, method = wass.method, niter = wass.iter, 
                          cost_a = cost_a, cost_b = cost_b, 
                          unbiased = unbiased, ...))
}

cv_eval <- function(weight, cv.list, cost, p, wass.method, wass.iter, sample_weight, 
                    cost_a, cost_b, unbiased, ...) {
  if (all(is.na(weight$w1) | all(is.na(weight$w0)))) return(NA_real_)
  n1 <- length(weight$w1)
  n0 <- length(weight$w0)
  
  p.temp <- weight$args$power
  if (is.null(p.temp)) {
    p <- list(...)[["p"]]
  } else {
    p <- p.temp
  }
  
  
  
  wp <- vapply(X = cv.list, FUN = cv_sep, FUN.VALUE = 1,
               weight = weight, cost = cost, p = p, wass.method = wass.method,
               wass.iter = wass.iter, sample_weight = sample_weight,
               cost_a = cost_a, cost_b = cost_b, unbiased = unbiased,
               ...)
  return(mean(wp^p))
}

cv_get <- function(n, K, R) {
  n <- as.integer(n)
  K <- as.integer(K)
  R <- as.integer(R)
  
  stopifnot(R >= 1)
  stopifnot(K > 1)
  stopifnot(n > K)
  
  sep.samples <- function(s, K, sets) {
    lapply(1:K, function(k) s[sets == k])
  }
  
  samples <- replicate(R, 
                       expr = sample.int(n = n, size = n, replace = FALSE, prob = NULL, useHash = FALSE),
                       simplify = FALSE)
  sets <- rep(1:K, length.out = n)
  
  folds <- unlist(lapply(samples, sep.samples, K = K, sets = sets), recursive = FALSE)
  
  return(folds)
}

wass_cv <- function(weights, K, R, cost, p,
                    wass.method, 
                    wass.iter, 
                    x0, x1,
                    cost_a = NULL, cost_b = NULL, 
                    unbiased = FALSE, verbose, 
                     ...) {
  
  
  
  if (verbose ) {
    message("\nCalculating out of sample balance via cross validation:")
    # pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
    pbapply::pboptions(type = "timer", style = 3, char = "=")
  } else {
    pbapply::pboptions(type = "none")
  }
  
  if ( is.list(cost) ) {
    cc <- cost[[length(cost)]]#[rowCount > 0, colCount > 0]
  } else {
    cc <- cost#[rowCount > 0, colCount > 0]
  }
  
  entropy.meth.sel <- wass.method %in% c("sinkhorn","greenkhorn")
  
  if (is.null(cost_a) && isTRUE(unbiased) && entropy.meth.sel) {
    n0 <- nrow(x0)
    cost_a <- cost_fun(rbind(x0,x0), z = c(rep(1,n0),
                                           rep(0,n0)),
                       power = p,
                       metric = metric, estimand = "ATT")
  } else {
    if (is.list(cost_a)) cost_a <- cost_a[[1]]
  }
  
  if (is.null(cost_b) && isTRUE(unbiased) && entropy.meth.sel) {
    n1 <- nrow(x1)
    cost_b <- cost_fun(rbind(x1,x1), z = c(rep(1, n1),
                                           rep(0, n1)),
                       power = p,
                       metric = metric, estimand = "ATT")
  } else {
    if (is.list(cost_b)) cost_b <- cost_b[[1]]
  }
  
  n_cv    <- ncol(cc)
  
  cv.list <- cv_get(n_cv, K, R)
  m_wp    <- pbapply::pbsapply(X = weights, FUN = cv_eval, 
                    cv.list = cv.list, cost = cc, p = p,
                    wass.method = wass.method, wass.iter = wass.iter, 
                    cost_a = cost_a, cost_b = cost_b, 
                    unbiased = unbiased, ...)
  
  return(m_wp)
  
}


wass_grid_eval <- function(data, grid = NULL, 
                             grid.length = 7,
                             estimand = c("ATT", "ATC","cATE","ATE"),
                             K = 10, R = 10,
                             n.boot = 1000,
                             eval.method = c("cross.validation", "bootstrap"),
                             method = c("Wasserstein","Constrained Wasserstein", "SCM"),
                             sample_weight = NULL,
                             wass.method = "networkflow", wass.iter = 0,
                             add.joint = TRUE,
                             add.margins = FALSE,
                             joint.mapping = FALSE,
                             verbose = FALSE,
                             neg.weights = FALSE,
                             cgd = FALSE,
                             ...) 
{

  get_calc_dist <- function(w, cost) {
    (sum(w$gamma * cost^(w$args$power)))^(1/(w$args$power))
  }
  
  get_outcome <- function(w, data, estimand) {
    estimate_effect(data = data, weights = w, hajek = TRUE,
                    doubly.robust = FALSE, estimand = estimand,
                    split.model = TRUE, matched = FALSE
                    )$estimate
  }
  
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  dots <- list(...)
  solver <- dots[["solver"]]
  p      <- dots[["p"]]
  metric <- dots[["metric"]]
  penalty <- dots[["penalty"]]
  eval.method <- match.arg(eval.method, c("cross.validation", "bootstrap"))
  if (is.null(solver)) solver <- "gurobi"
  if (is.null(p)) p <- 2
  if (is.null(metric)) metric <- "mahalanobis"
  if (is.null(penalty)) penalty <- "L2"
  add.margins <- isTRUE(add.margins)
  add.joint  <- isTRUE(add.joint)
  joint.mapping <- isTRUE(joint.mapping)
  neg.weights <- isTRUE(neg.weights)
  cgd <- isTRUE(cgd)
  
  if (method == "SCM") {
    add.margins <- add.joint <- joint.mapping <- FALSE
  }
  
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
  
  # wass.dat <- cbind(z = z, x)
  
  cost <- wg_cost_setup(cost = dots$cost, p = p,
                        estimand = estimand, x = x, z = z,
                        metric = metric, method = method,
                        rkhs.args = dots$rkhs.args,
                        add.margins = add.margins)
  
  
  if (all(is.null(grid)) | all(is.na(grid))) {
    grid <- f.call.list(fun = "wass_grid_default", 
                        list.args = list(x = x, z = z, 
                                         grid.length = grid.length,
                                         p = p, data = data, cost = cost, estimand = estimand,
                                         method = method, metric = metric, wass.iter = wass.iter, 
                                         add.joint = add.joint, add.margins = add.margins,
                                         joint.mapping = joint.mapping, penalty = penalty,
                                         ...))
  }
  
  args <- list(data = data, 
               x0 = x0,
               x1 = x1,
               x = x,
               z = z,
               grid = grid,  
               n.boot = n.boot,
               K = K, 
               R = R,
               eval.method = "cross.validation",
               wass.method = wass.method,
               wass.iter = wass.iter,
               sample_weight = sample_weight,
               estimand = estimand, 
               method = method, 
               solver = solver, 
               metric = metric,
               p = p, 
               cost = cost, 
               add.joint = add.joint,
               add.margins = add.margins, 
               joint.mapping = joint.mapping,
               neg.weights = neg.weights,
               cgd = cgd, 
               verbose = verbose,
               ...)
  
  weight.list <- weight_est_fun(args)
  
  taus <- vapply(X = weight.list,
                     FUN = get_outcome, FUN.VALUE = 1,
                     data = data, estimand = estimand)
  
  output.weight <- eval_weights(weight.list, args)
  
  cv.chosen <- output.weight$args$constraint
  cv.chosen$margins <- cv.chosen$margins[1]
  cv.chosen$joint <- cv.chosen$joint
  cv.chosen$penalty <- cv.chosen$penalty
  
  if (is.null(cv.chosen$margins)) cv.chosen$margins <- NA_real_
  if (is.null(cv.chosen$joint)) cv.chosen$joint <- NA_real_
  if (is.null(cv.chosen$penalty)) cv.chosen$penalty <- NA_real_
  
  args$eval.method <- "bootstrap"
  output.weight.b <- eval_weights(weight.list, args)
  
  boot.chosen <- output.weight.b$args$constraint
  boot.chosen$margins <- boot.chosen$margins[1]
  boot.chosen$joint <- boot.chosen$joint
  boot.chosen$penalty <- boot.chosen$penalty
  
  if (is.null(boot.chosen$margins)) boot.chosen$margins <- NA_real_
  if (is.null(boot.chosen$joint)) boot.chosen$joint <- NA_real_
  if (is.null(boot.chosen$penalty)) boot.chosen$penalty <- NA_real_
  
  grid.margins <- unlist(sapply(grid, function(x) x$margins[1]))
  grid.joint   <- unlist(sapply(grid, function(x) x$joint[1]))
  grid.penalty <- unlist(sapply(grid, function(x) x$penalty[1]))
  
  if (is.null(grid.margins)) grid.margins <- NA_real_
  if (is.null(grid.joint)) grid.joint <- NA_real_
  if (is.null(grid.penalty)) grid.penalty <- NA_real_
  
  overall.cost <- if (is.list(cost)) {
    cost[[length(cost)]]
  } else {
    cost
  }
  
  if (verbose) message("Evaluating Wasserstein distance on full sample")
  evaluated.distance <- vapply(X = weight.list,  FUN = function(X) wass_dist_helper(a = X,
                                                                                    p = p,
                                                                                    cost = overall.cost, method = wass.method,
                                                                                    niter = wass.iter, ...)
                               , FUN.VALUE = 1)
  
  wp  <- vapply(X = weight.list,  FUN = get_calc_dist,
                FUN.VALUE = 1,
                cost = overall.cost)
  
  weight.var <- vapply(X = weight.list, FUN = function(w) var(c(w$gamma)), FUN.VALUE = 1)
  weight.var.a <- vapply(X = weight.list, FUN = function(w) var(w$w0), FUN.VALUE = 1)
  weight.var.b <- vapply(X = weight.list, FUN = function(w) var(w$w1), FUN.VALUE = 1)
  
  weight.ent <- vapply(X = weight.list, FUN = function(w) entropy(c(w$gamma)), FUN.VALUE = 1)
  weight.ent.a <- vapply(X = weight.list, FUN = function(w) entropy(w$w0), FUN.VALUE = 1)
  weight.ent.b <- vapply(X = weight.list, FUN = function(w) entropy(w$w1), FUN.VALUE = 1)
  
  
  output <- data.frame(n = n, d = ncol(x),
                       estimates = taus,
                       cv.chosen.margins = cv.chosen$margins,
                       cv.chosen.joint = cv.chosen$joint, 
                       cv.chosen.penalty = cv.chosen$penalty,
                       boot.chosen.margins = boot.chosen$margins,
                       boot.chosen.joint = boot.chosen$joint, 
                       boot.chosen.penalty = boot.chosen$penalty,
                       grid.margins = grid.margins,
                       grid.joint  = grid.joint,
                       grid.penalty  = grid.penalty,
                       reg.wp = evaluated.distance, wp = wp,
                       var.gamma = weight.var,
                       var.a = weight.var.a,
                       var.b = weight.var.b,
                       ent.gamma = weight.ent,
                       ent.a = weight.ent.a,
                       ent.b = weight.ent.b,
                       p = p, eval.method = eval.method,
                       wass.method = wass.method,
                       wass.iter = wass.iter,
                       method = method,
                       metric = metric,
                       penalty = penalty,
                       estimand = estimand)
  
  
  return(output)
  
}
