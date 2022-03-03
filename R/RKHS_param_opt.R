## TODO: add linear kernel...

setClass("RKHS_param", slots = c(theta = "numeric", 
                                 gamma = "numeric",
                                 p = "numeric",
                                 sigma_2 = "numeric",
                                 kernel = "character",
                                 metric = "character",
                                 is.dose = "logical"
                                 
                                 ),
                       prototype = list(theta = numeric(0),
                                        gamma = numeric(0),
                                        p = numeric(0),
                                        sigma_2 = numeric(0),
                                        kernel = character(0),
                                        metric = character(0),
                                        is.dose = logical(0)))

RKHS_param_opt <- function(x, y, z, power = 2:3, metric = c("mahalanobis", "Lp"), is.dose = FALSE, 
                           opt.method = c("stan", "optim"), 
                           kernel = c("RBF","polynomial","linear"),
                           estimand = c("ATC","ATT","ATE"), ...) {
  
  opt.method <- match.arg(opt.method)
  metric <- match.arg(metric)
  estimand <- match.arg(estimand)
  kernel <- match.arg(kernel)
  
  is.dose <- isTRUE(is.dose)
  power <- as.integer(power)
  
  y_std <- rep(NA, length(y))
  y_std[z==0] <- scale(y[z==0])
  y_std[z==1] <- scale(y[z==1])
  
  similarity_mats <- calc_similarity(x, z, metric = metric, kernel = kernel, is.dose, estimand)
  
  kern_dose <- function(theta_0, theta_1, gamma_0, gamma_1, sigma2) {
    K <- gamma_0*(1.0 + theta_0 * similarity_mats[["Z"]])^p *  gamma_1*(1.0 + theta_1 * similarity_mats[["X"]])^p
    n <- nrow(K)
    return(K + diag(sigma2, n, n))
  }
  kern_ <- function(theta_0, theta_1, gamma_0, gamma_1, sigma2) {
    s2 <- as.double(rep(sigma2, length(z)))
    return(kernel_update_(sim_ = similarity_mats, z_ = z, theta_ = c(theta_0, theta_1),
                          gamma_ =c(gamma_0, gamma_1), p = p, sigma_2_ = s2,
                          kernel_ = kernel))
    
  }
  cmb_kern <- switch(as.character(is.dose),
                     "TRUE" = kern_dose,
                     "FALSE" = kern_)
  
  # if(opt.method == "bayesian.optimization") {
  #   if(kernel != "polynomial") {
  #     kernel <- "polynomial"
  #     warning("RBF and linear kernels aren't supported for method optim. Switching to polynomial kernel")
  #   }
  #   scoring_fun <- function(theta_0, theta_1, gamma_0, gamma_1, sigma2) {
  #     K <- cmb_kern(theta_0, theta_1, gamma_0, gamma_1, sigma2)
  #     score <- marginal_lik_gp_(y_std, K)
  #     if(is.nan(score)) score <- -.Machine$double.xmax
  #     if(is.infinite(score)) score <- -.Machine$double.xmax
  #     return(list(Score = score))
  #   }
  #   args <- list(FUN = scoring_fun, ...)
  #   args <- args[!duplicated(names(args))]
  #   
  #   if(is.null(args$bounds)) args$bounds <- list(theta_0 = c(.Machine$double.xmin,100),
  #                                                theta_1 = c(.Machine$double.xmin,100),
  #                                                gamma_0 = c(.Machine$double.xmin,1000),
  #                                                gamma_1 = c(.Machine$double.xmin,1000),
  #                                                sigma2 = c(.Machine$double.xmin,1))
  #   if(is.null(args$initPoints)) args$initPoints <- 10
  #   if(is.null(args$iters.n)) args$iters.n <- 5
  #   if(is.null(args$iters.k)) args$iters.k <- 1
  #   res <- vals <- vector("list", length(power))
  #   param <- list()
  #   for(i in seq_along(power) ) {
  #     p <- power[i]
  #     res[[i]] <- do.call(ParBayesianOptimization::bayesOpt, args)
  #     param <- ParBayesianOptimization::getBestPars(res[[i]])
  #     vals[[i]] <- do.call(scoring_fun, param)
  #   }
  #   idx <- which.max(unlist(vals))
  #   param <- ParBayesianOptimization::getBestPars(res[[idx]])
  #   param$p <- power[idx]
  # } 
  # else 
  if(opt.method == "optim") {
    if(kernel != "polynomial" ) {
      kernel <- "polynomial"
      warning("RBF and linear kernels aren't supported for method optim. Switching to polynomial kernel")
    }
    scoring_fun <- function(param) {
      exp.param <- exp(param)
      K <- cmb_kern(exp.param["theta_0"], 
                    exp.param["theta_1"], 
                    exp.param["gamma_0"], 
                    exp.param["gamma_1"], 
                    exp.param["sigma2"])
      score <- marginal_lik_gp_(y_std, K)
      if(is.nan(score)) score <- -.Machine$double.xmax
      return(-score)
    } 
    # grad_fun <- function(theta_0, theta_1, gamma_0, gamma_1, p, sigma2) {
    #   K <- cmb_kern(theta_0, theta_1, gamma_0, gamma_1, p, sigma2)
    #   val <- marg_lik_gp_grad_(y, K)
    #   return(-val)
    # }
    
    args <- list(par = c(theta_0 = 0, theta_1 = 0, 
                         gamma_0=0, gamma_1=0, sigma2=0),
                 fn = scoring_fun, gr = NULL, ...)
    args <- args[!duplicated(names(args))]
    args <- args[names(args) %in% names(formals(optim))]
    res <- vector("list", length(power))
    
    for(i in seq_along(power) ) {
      p <-  power[i]
      res[[i]] <- do.call("optim", args)
      if(res[[i]]$convergence !=0) {
        if(res[[i]]$convergence ==1) {
          warning("max iter hit")
        }
      }
      # samp <- rstan::vb(stan_mod, data = data_stan)
    }
    # res <- do.call("optim", list(par = c(theta_0 = 0, theta_1 = 0, 
    #                                      gamma_0=0, gamma_1=0, sigma2=0),
    #                              fn = scoring_fun, gr = NULL, ...))
    idx <- which.min(sapply(res, function(r) r$value[1]))
    param <- lapply(exp(res[[idx]]$par), function(r) r)
    param$p <- power[idx]
  } 
  else if (opt.method == "stan") {
    if(is.list(similarity_mats)) {
      dis_x <- similarity_mats[["X"]]
      dis_z <- similarity_mats[["Z"]]
    } else {
      dis_x <- similarity_mats
      dis_z <- matrix(0, nrow(dis_x), ncol(dis_x))
    }
    
    arguments <- list(
      data = list(N = as.integer(length(y)),
                  y = c(y_std),
                  discrep = dis_x,
                  discrep_z = dis_z,
                  z = z,
                  p = power[1],
                  is_dose = as.integer(is.dose),
                  kernel = switch(kernel,
                                  "RBF" = 2L,
                                  "polynomial" = 1L,
                                  "linear" = 3L,
                                  2L)),
      ...,
      as_vector = FALSE
    )
    if(is.null(arguments$data$kernel)) arguments$data$kernel <- 2L
    formals.stan <- c("iter", "save_iterations", 
                      "refresh", "init_alpha", "tol_obj", "tol_grad", "tol_param", 
                      "tol_rel_obj", "tol_rel_grad", "history_size",
                      "object", "data", "seed", "init", "check_data",
                      "sample_file", "algorithm", "verbose", "hessian", 
                      "as_vector", "draws", "constrained", "importance_resampling")
    arguments <- arguments[!duplicated(names(arguments))]
    arguments <- arguments[names(arguments) %in% formals.stan]
    if(is.null(arguments$object)) {
      if(!is.dose) {
        arguments$object <- stanbuilder("gp_hyper")
      } else{
        arguments$object <- stanbuilder("gp_hyper_dose")
      }
    }
    if(is.null(arguments$algorithm) & kernel == "RBF") {
      arguments$algorithm <- "Newton"
    }
    
    argn <- lapply(names(arguments), as.name)
    names(argn) <- names(arguments)
    f.call <- as.call(c(list(call("::", as.name("rstan"), 
                                  as.name("optimizing"))), argn))
    # res <- eval(f.call, envir = arguments)
    if (kernel == "polynomial") {
      res <- vector("list", length(power))
      for(i in seq_along(power) ) {
        arguments$data$p <-  as.double(power[i])
        res[[i]] <- if(is.dose) {
          eval(f.call, envir = arguments)
        } else {
          rkhs_stan_binary_helper(f.call, arguments)
        }
        # samp <- rstan::vb(stan_mod, data = data_stan)
      }
      if (is.dose) {
        idx <- which.max(sapply(res, function(r) r$par$marg_lik))
        param <- res[[idx]]$par
      } else {
        idx <- which.max(sapply(res, function(r) r$marg_lik))
        param <- res[[idx]]
      }
      # param$sigma2 <- c(param$sigma_0, param$sigma_1)
      param$p <- power[idx]
    } else if (kernel == "RBF" | kernel == "linear") {
      if (!is.dose) {
        param <- rkhs_stan_binary_helper(f.call, arguments)
      } else {
        res <- eval(f.call, envir = arguments)
        param <- res$par
      }
    }
  }
  out <- list(theta = c(param$theta_0, param$theta_1),
              gamma = c(param$gamma_0, param$gamma_1),
              p = param$p,
              sigma_2 = param$sigma,
              kernel = kernel,
              metric = metric,
              is.dose = is.dose,
              is.standardized = TRUE)
  if (is.null(out$p)) out$p <- NA_real_
  if (is.null(out$theta)) out$theta <- c(NA_real_, NA_real_)
  if (is.null(out$gamma)) out$gamma <- c(NA_real_, NA_real_)
  class(out) <- "RKHS_param"
  return(out)
}

rkhs_stan_binary_helper <- function(f.call, arguments) {
  param <- list()
  tune.fun <- function(x, l, u){
    abs(0.01 - pgamma(l, exp(x[1]), exp(x[2]))) + abs(0.01 - pgamma(u, exp(x[1]), exp(x[2]), lower.tail = FALSE))
  }
  
  y_std <- arguments$data$y
  z <- arguments$data$z
  dis_x <- arguments$data$discrep
    
  #control gp
  arguments$data$N = sum(z==0)
  arguments$data$y = c(y_std[z==0])
  arguments$data$discrep = dis_x[z==0,z==0]
  
  arguments$data$a <- 0.
  arguments$data$b <- 0.
  
  if(arguments$data$kernel == 2) {
    l = sqrt(min(arguments$data$discrep[lower.tri(arguments$data$discrep)]))
    u = sqrt(max(arguments$data$discrep))
    tuned <- optim(par = runif(2), fn = tune.fun, l = l, u = u)
    if (tuned$convergence == 0) {
      arguments$data$a <- exp(tuned$par[1])
      arguments$data$b <- exp(tuned$par[2])
    }
    if (arguments$algorithm != "Newton") {
      res <- lapply(1:10,function(i) eval(f.call, envir = arguments))
      res <- res[[which.max(sapply(res, function(r) r$par$marg_lik))]]
    } else {
      res <- eval(f.call, envir = arguments)
    }
  } else {
    res <- eval(f.call, envir = arguments)
  }
  if (!is.null(res$par$theta_half)) {
    param$theta_0 <- if (arguments$data$kernel == 2) {
      1/res$par$theta_half^2
    } else {
      res$par$theta_half^2
    }
  }
  if (!is.null(res$par$gamma_half)) param$gamma_0 <- res$par$gamma_half^2
  param$sigma   <- res$par$sigma_half^2
  param$marg_lik <- res$par$marg_lik
  
  #treated gp
  arguments$data$N = sum(z==1)
  arguments$data$y = c(y_std[z==1])
  arguments$data$discrep = dis_x[z==1,z==1]
  
  arguments$data$a <- 0.
  arguments$data$b <- 0.
  
  if(arguments$data$kernel == 2) {
    l = sqrt(min(arguments$data$discrep[lower.tri(arguments$data$discrep)]))
    u = sqrt(max(arguments$data$discrep))
    tuned <- optim(par = runif(2), fn = tune.fun, l = l, u = u)
    if(tuned$convergence == 0) {
      arguments$data$a <- exp(tuned$par[1])
      arguments$data$b <- exp(tuned$par[2])
    }
    if(arguments$algorithm != "Newton"){
      res <- lapply(1:10,function(i) eval(f.call, envir = arguments))
      res <- res[[which.max(sapply(res, function(r) r$par$marg_lik))]]
    } else {
      res <- eval(f.call, envir = arguments)
    }
  } else {
    res <- eval(f.call, envir = arguments)
  }
    
  
  if(!is.null(res$par$theta_half)) {
    param$theta_1 <- if(arguments$data$kernel == 2) {
      1/res$par$theta_half^2
    } else {
      res$par$theta_half^2
    }
  }
  if(!is.null(res$par$gamma_half)) param$gamma_1 <- res$par$gamma^2
  param$sigma   <- c(param$sigma, res$par$sigma^2)
  param$marg_lik <- res$par$marg_lik + param$marg_lik
  return(param)
}
  