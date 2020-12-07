sbw_grid_search <- function(data, grid = NULL, 
                            estimand = c("ATT", "ATC","cATE","ATE","feasible"),
                            n.boot = 100,
                            ...) 
{
  # if(is.null(grid) & !is.null(list(...)$constraint)) grid <- constraint
  if(all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 1/sqrt(get_p(data)), length.out = 100)
  estimand <- match.arg(estimand)
  
  args <- list(data = data, constraint = grid[1],  estimand = estimand, 
       method = "SBW",
       ...)
  args <- args[!duplicated(names(args))]
  argn <- lapply(names(args), as.name)
  names(argn) <- names(args)
  
  f.call <- as.call(c(as.name("calc_weight_bal"), argn))
  
  pd <- prep_data(data, ...)
  x1 <- as.matrix(pd$df[pd$z == 1,-1])
  x0 <- as.matrix(pd$df[pd$z == 0,-1])
  x <- rbind(x0,x1)
  mean.bal.dat <- cbind(z = pd$z, x)
  
  n <- nrow(pd$df)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
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
  weight.list[[min.idx]]$addl.args$standardized.mean.difference <- grid[min.idx]
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
  return(weight.list[[min.idx]])
}

mean_bal_grid <- function(bootIdx, weight, data, tx_ind,...) {
  wvec <- c(weight$w0, weight$w1)
  dataResamp <- data[bootIdx,]
  weightResamp <- wvec[bootIdx]
  z <- dataResamp[,tx_ind]
  wl <- list(w0 = weightResamp[z == 0], w1 = weightResamp[z == 1])
  bals <- mean_bal(dataResamp, weights = wl, treatment.indicator = tx_ind, ...)
  return(mean(bals))
}
