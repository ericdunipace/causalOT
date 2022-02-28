confint.causalWeights <- function(object, parm, level = 0.95, data, method = c("bootstrap", "asymptotic"),...) {
  
 if(method == "bootstrap") {
   return(ci_boot_cw(data,object, parm, level, ...))
 } else {
   stop("asymptotic confidence intervals not currently implemented")
 }
  
  
}


ci_boot_cw <- function(data, object, parm, level, n.boot = 1000, balance.covariates = NULL,
                    treatment.indicator = NULL, outcome = NULL, ...) {
  boot.fun <- function(idx, data, object,  ...) {
    data <- get_subset(data, idx, ...)
    
    wt.args <- c(list(data = data, constraint = object$args$constraint,
                      estimand = object$estimand, method = object$method,
                      transport.matrix = !is.null(object$gamma),
                      grid.search = isTRUE(object$args$grid.search), 
                      formula = object$args$formula), 
                 object$args,
                 ...)
    
    wt.args <- wt.args[!duplicated(names(wt.args))]
    wtargn <- names(wt.args)
    wf.call <- as.call(c(list(as.name("calc_weight")), wtargn))
    
    weight <- eval(wf.call, envir = wt.args) 
    
    est.args <- c(list(data = data, weights = weight),
                  ...)
    
    est.args <- est.args[!duplicated(names(est.args))]
    estargn <- names(est.args)
    estf.call <- as.call(c(list(as.name("estimate_effect")), estargn))
    
    est <- eval(estf.call, envir = est.args)$estimate
    clear_subset(data)
    return(est)
  }
  ns <- get_n(data, ...)
  N  <- sum(ns)
  
  boot.idx <- replicate(n = n.boot, sample.int(N,N, replace = TRUE))
  
  boots <- lapply(boot.idx, boot.fun, data = data, object = object, ...)
  
  ci <- quantile(boots, probs = c((1-level)/2, level + (1-level)/2))
  names(ci) <- c("lower", "upper")
  
  return(ci)
  
}



