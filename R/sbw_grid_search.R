sbw_grid_search <- function(data, grid = NULL, 
                            estimate = c("ATT", "ATC","ATE","feasible"),
                            n.boot = 100,
                            ...) 
{
  if(all(is.null(grid)) | all(is.na(grid))) grid <- seq(0, 0.5, length.out = 100)
  weight.list <- lapply(grid, function(delta) {
    out <- do.call("calc_weight_bal", list(data = data, constraint = delta,  estimate = estimate, 
                                             method = "SBW",
                                             ...))
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
  
  min.idx <- which(output == min(output, na.rm=TRUE))
  weight.list[[min.idx]]$standardized.mean.difference <- grid[min.idx]
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
