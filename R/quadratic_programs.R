quadprog.default <- function(x, z, y = NULL, constraint,  estimand = c("ATT", "ATC", "ATE","cATE","feasible"), 
                             method = supported.methods(),
                             sample_weight = NULL,
                             soc = FALSE,
                             ...) {
  meth <- match.arg(method, supported.methods())
  est <- match.arg(estimand)
  soc <- isTRUE(soc)
  
  dots <- list(...)
  
  form <- dots[["formula"]]
  
  sw <- get_sample_weight(sample_weight, z)
  
  if (est == "ATE") {
    sw.comb <- ate_sample_weight(sw)
    sw0 <- sw.comb[[1]]
    sw1 <- sw.comb[[2]]
  } else if (est == "cATE") {
    sw0 <- switch_sample_weight(sw,z)
    sw1 <- sw
  } else if (est == "ATC") {
    sw <- switch_sample_weight(sw,z)
  }
  
  if ( isTRUE(!is.null(form)) & isTRUE(!is.na(form)) ) {
    form <- form_all_squares(form, colnames(x))
    
    form.temp <- as.character(form[length(form)])
    form <- as.formula(paste0("~ 0 +", form.temp))
    mm <- model.matrix(form, data = data.frame(x))
    if ( all(mm[,1] == 1)) mm <- mm[,-1]
    if (method == "Wasserstein" | method == "Constrained Wasserstein" | method == "SCM") {
      bf <- list(mm = mm,
                 K = dots[["balance.constraints"]])
      if (is.null(bf$K)) {
        bf <- NULL
        warning("Formula provided but not constraints. Not using SBW constraints")
      }
    }
  } else {
    if (method == "SBW") {
      form <- formula(~ . + 0)
      mm <- model.matrix(form, data = data.frame(x))
    } else {
      bf <- NULL
    }
  }
  # margmass <- get_sample_weight(sample_weight, z = data$get_z())
  
  qp <- if (meth == "SBW") {
    # if (length(constraint) != ncol(mm))  {
    #   K <- rep(constraint, ncol(mm))[1:ncol(mm)]
    # } else {
    K <- constraint
    # }
    if (est == "cATE") {
      list(qp_sbw(x = mm, z = z, K = K, estimand = "ATC"),
           qp_sbw(x = mm, z = z, K = K, estimand = "ATT"))
    } else if (est == "ATE") {
      qp_sbw(x = mm, z = z, K = K, estimand = est)
    } else {
      list(qp_sbw(x = mm, z = z, K = K, estimand = est))
    } 
  } else if (meth == "SCM") {
    if (is.null(dots[["p"]])) dots[["p"]] <- 2
    if (is.null(dots[["metric"]])) dots[["metric"]] <- "mahalanobis"
    if (is.null(dots[["add.margins"]])) dots[["add.margins"]] <- FALSE
    if (is.null(dots[["penalty"]])) dots[["penalty"]] <- "L2"
    if (is.null(dots[["joint.mapping"]])) dots[["joint.mapping"]] <- FALSE
    if (is.null(dots[["neg.weights"]])) dots[["neg.weights"]] <- FALSE
    
    
    if (est == "cATE") {
      if(is.list(constraint)) {
        if(length(constraint) == 2 && !any(names(constraint) %in% c("joint","penalty","margins")) ) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
        
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_scm(x = x, z = 1 - z, K = constraint0,
                  p=dots[["p"]], dist = dots[["metric"]], 
                  cost = dots[["cost"]],
                  penalty = dots[["penalty"]],
                  joint.mapping = dots[["joint.mapping"]],
                  rkhs.args = dots[["rkhs.args"]], 
                  add.margins = dots[["add.margins"]], 
                  neg.weights = dots[["neg.weights"]],
                  bf = bf,
                  sample_weight = sw0,
                  soc = soc),
           qp_scm(x=x, z = z, K = constraint1,
                  p=dots[["p"]], dist = dots[["metric"]], cost = dots[["cost"]],
                  penalty = dots[["penalty"]],
                  joint.mapping = dots[["joint.mapping"]],
                  rkhs.args = dots[["rkhs.args"]], 
                  add.margins = dots[["add.margins"]], 
                  neg.weights = dots[["neg.weights"]],
                  bf = bf,
                  sample_weight = sw1,
                  soc = soc))
    } else if (est == "ATE") {
      bf0 <- bf1 <- bf
      if (!is.null(bf)) {
        bf0$mm <- rbind(bf$mm[z == 0, ,drop = FALSE], bf$mm)
        bf1$mm <- rbind(bf$mm[z == 1, ,drop = FALSE], bf$mm)
      }
      if(is.list(constraint)) {
        if(length(constraint) == 2) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_scm(x = rbind(x[z == 0, ,drop = FALSE],x), 
                  z = c(rep(0, sum(1 - z)), rep(1, nrow(x))),
                  K = constraint0,
                  p = dots[["p"]], dist = dots[["metric"]], 
                  cost = dots[["cost"]][[1]],
                  penalty = dots[["penalty"]],                   
                  joint.mapping = dots[["joint.mapping"]],
                  rkhs.args = dots[["rkhs.args"]], 
                  add.margins = dots[["add.margins"]], 
                  neg.weights = dots[["neg.weights"]],
                  bf = bf0,
                  sample_weight = sw0,
                  soc = soc),
           qp_scm(x = rbind(x[z == 1, ,drop = FALSE],x), 
                  z = c(rep(0, sum(z)), rep(1, nrow(x))),
                  K = constraint1,
                  p = dots[["p"]], dist = dots[["metric"]], 
                  cost = dots[["cost"]][[2]],
                  penalty = dots[["penalty"]],
                  joint.mapping = dots[["joint.mapping"]],
                  rkhs.args = dots[["rkhs.args"]], 
                  add.margins = dots[["add.margins"]], 
                  neg.weights = dots[["neg.weights"]],
                  bf = bf1,
                  sample_weight = sw1,
                  soc = soc))
    } else if (est == "ATC") {
      list(qp_scm(x=x, z=1-z, K=constraint, 
                  p=dots[["p"]], dist = dots[["metric"]], 
                  cost = dots[["cost"]], 
                  add.margins = dots[["add.margins"]],
                  penalty = dots[["penalty"]],
                  joint.mapping = dots[["joint.mapping"]],
                  add.joint = dots[["add.joint"]],
                  rkhs.args = dots[["rkhs.args"]],  
                  neg.weights = dots[["neg.weights"]],
                  bf = bf,
                  sample_weight = sw,
                  soc = soc))
    } else {
      list(qp_scm(x = x, z = z, K = constraint,
                  p = dots[["p"]], dist = dots[["metric"]], 
                  cost = dots[["cost"]],
                  rkhs.args = dots[["rkhs.args"]], 
                  penalty = dots[["penalty"]],
                  joint.mapping = dots[["joint.mapping"]],
                  add.margins = dots[["add.margins"]], 
                  neg.weights = dots[["neg.weights"]],
                  bf = bf,
                  sample_weight = sw,
                  soc = soc))
    }
  # } else if (meth == "EBW") {
  #     if (is.null(dots[["p"]])) dots[["p"]] <- 2
  #     if (is.null(dots[["metric"]])) dots[["metric"]] <- "mahalanobis"
  #     
  #     if (est == "cATE") {
  #       list(qp_ebw(x = x, z = 1 - z,
  #                   p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw0),
  #            qp_ebw(x=x, z = z, p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw1))
  #     } else if (est == "ATE") {
  #       
  #       list(qp_ebw(x = rbind(x[z == 0, ,drop = FALSE],x), 
  #                   z = c(rep(0, sum(1 - z)), rep(1, nrow(x))),
  #                   p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw0),
  #            qp_ebw(x = rbind(x[z == 1, ,drop = FALSE],x), 
  #                   z = c(rep(0, sum(z)), rep(1, nrow(x))),
  #                   p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw1))
  #     } else if (est == "ATC") {
  #       list(qp_ebw(x=x, z= 1 - z, 
  #                   p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw))
  #     } else {
  #       list(qp_ebw(x = x, z = z, p = dots[["p"]],
  #                   dist = dots[["metric"]], cost = dots[["cost"]],
  #                   sample_weight = sw))
  #     }
  } else if (meth == "Wasserstein") {
    
    # if (is.null(dots[["p"]])) dots[["p"]] <- 2
    if (is.null(dots[["metric"]])) dots[["metric"]] <- "mahalanobis"
    if (is.null(dots[["add.margins"]])) dots[["add.margins"]] <- FALSE
    if (is.null(dots[["penalty"]])) dots[["penalty"]] <- "L2"
    if (is.null(dots[["joint.mapping"]])) dots[["joint.mapping"]] <- FALSE
    if (is.null(dots[["neg.weights"]])) dots[["neg.weights"]] <- FALSE
    if (is.null(dots[["add.divergence"]])) dots[["add.divergence"]] <- FALSE
    
    if (dots[["metric"]] == "RKHS" & is.null(dots[["rkhs.args"]])) {
      if (is.null(dots[["opt.method"]])) dots[["opt.method"]] <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
      )
      if (is.null(dots[["kernel"]])) dots[["kernel"]] <- "RBF"
      args <- list(x = x, 
                   z = z, 
                   y = y,
                   metric = switch(dots[["metric"]],
                                   "RKHS" = "mahalanobis",
                                   "mahalanobis" = "mahalanobis",
                                   "Lp" = "Lp"),
                   kernel = dots[["kernel"]],
                   is.dose = dots[["is.dose"]],
                   opt.method = dots[["opt.method"]],
                   estimand = temp.est,
                   ...)
      args <- args[!duplicated(names(args))]
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
      dots[["rkhs.args"]] <- eval(f.call, args)
    }
    if (est == "cATE") {
      if(is.list(constraint)) {
        if(length(constraint) == 2 && !any(names(constraint) %in% c("joint","penalty","margins"))) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_wass(x = x, z = 1-z, K = constraint0,
                   p=dots[["p"]], dist = dots[["metric"]], 
                   cost = dots[["cost"]],
                   penalty = dots[["penalty"]],
                   joint.mapping = dots[["joint.mapping"]],
                   rkhs.args = dots[["rkhs.args"]], 
                   add.margins = dots[["add.margins"]], 
                   neg.weights = dots[["neg.weights"]],
                   bf = bf,
                   sample_weight = sw0,
                   soc = soc, divergence = dots[["add.divergence"]]),
           qp_wass(x=x, z = z, K = constraint1,
                   p=dots[["p"]], dist = dots[["metric"]], cost = dots[["cost"]],
                   penalty = dots[["penalty"]],
                   joint.mapping = dots[["joint.mapping"]],
                   rkhs.args = dots[["rkhs.args"]], 
                   add.margins = dots[["add.margins"]], 
                   neg.weights = dots[["neg.weights"]],
                   bf = bf,
                   sample_weight = sw1,
                   soc = soc, divergence = dots[["add.divergence"]]))
    } else if (est == "ATE") {
      bf0 <- bf1 <- bf
      if (!is.null(bf)) {
        bf0$mm <- rbind(bf$mm[z == 0, ,drop = FALSE], bf$mm)
        bf1$mm <- rbind(bf$mm[z == 1, ,drop = FALSE], bf$mm)
      }
      if(is.list(constraint)) {
        if(length(constraint) == 2) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_wass(x = rbind(x[z == 0, ,drop = FALSE],x), 
                   z = c(rep(0, sum(1 - z)), rep(1, nrow(x))),
                   K = constraint0,
                   p = dots[["p"]], dist = dots[["metric"]], 
                   cost = dots[["cost"]][[1]],
                   penalty = dots[["penalty"]],                   
                   joint.mapping = dots[["joint.mapping"]],
                   rkhs.args = dots[["rkhs.args"]], 
                   add.margins = dots[["add.margins"]], 
                   neg.weights = dots[["neg.weights"]],
                   bf = bf0,
                   sample_weight = sw0,
                   soc = soc, divergence = dots[["add.divergence"]]),
           qp_wass(x = rbind(x[z == 1, ,drop = FALSE],x), 
                   z = c(rep(0, sum(z)), rep(1, nrow(x))),
                   K = constraint1,
                   p = dots[["p"]], dist = dots[["metric"]], 
                   cost = dots[["cost"]][[2]],
                   penalty = dots[["penalty"]],
                   joint.mapping = dots[["joint.mapping"]],
                   rkhs.args = dots[["rkhs.args"]], 
                   add.margins = dots[["add.margins"]], 
                   neg.weights = dots[["neg.weights"]],
                   bf = bf1,
                   sample_weight = sw1,
                   soc = soc, divergence = dots[["add.divergence"]]))
    } else if (est == "ATC") {
      list(qp_wass(x=x, z=1-z, K=constraint, 
                   p=dots[["p"]], dist = dots[["metric"]], 
                   cost = dots[["cost"]], 
                   add.margins = dots[["add.margins"]],
                   penalty = dots[["penalty"]],
                   joint.mapping = dots[["joint.mapping"]],
                   add.joint = dots[["add.joint"]],
                   rkhs.args = dots[["rkhs.args"]],  
                   neg.weights = dots[["neg.weights"]],
                   bf = bf,
                   sample_weight = sw,
                   soc = soc, divergence = dots[["add.divergence"]]))
    } else {
      list(qp_wass(x = x, z = z, K = constraint,
                   p=dots[["p"]], dist = dots[["metric"]], 
                   cost = dots[["cost"]],
                   rkhs.args = dots[["rkhs.args"]], 
                   penalty = dots[["penalty"]],
                   joint.mapping = dots[["joint.mapping"]],
                   add.margins = dots[["add.margins"]],
                   neg.weights = dots[["neg.weights"]],
                   bf = bf,
                   sample_weight = sw,
                   soc = soc, divergence = dots[["add.divergence"]]))
    }
    
  } else if (meth == "Constrained Wasserstein") {
    dots <- list(...)
    # if(is.null(dots[["p"]])) dots[["p"]] <- 2
    if(is.null(dots[["metric"]])) dots[["metric"]] <- "mahalanobis"
    if (is.null(dots[["add.joint"]])) dots[["add.joint"]] <- TRUE
    if (is.null(dots[["add.margins"]])) dots[["add.margins"]] <- FALSE
    if(!dots[["add.joint"]] & !dots[["add.margins"]]) stop("must run margins or joint")
    if (is.null(dots[["penalty"]])) dots[["penalty"]] <- "L2"
    if (is.null(dots[["joint.mapping"]])) dots[["joint.mapping"]] <- FALSE
    if (is.null(dots[["neg.weights"]])) dots[["neg.weights"]] <- FALSE
    
    if(dots[["metric"]] == "RKHS" & is.null(dots[["rkhs.args"]])) {
      if(is.null(dots[["opt.method"]])) dots[["opt.method"]] <- "stan"
      temp.est <- switch(estimand,
                         "ATT" = "ATT",
                         "ATC" = "ATC",
                         "ATE"
      )
      if(is.null(dots[["kernel"]])) dots[["kernel"]] <- "RBF"
      args <- list(x=x, 
                   z=z, 
                   y = y,
                   metric = switch(dots[["metric"]],
                                   "RKHS" = "mahalanobis",
                                   "mahalanobis" = "mahalanobis",
                                   "Lp" = "Lp"),
                   kernel = dots[["kernel"]],
                   is.dose = dots[["is.dose"]],
                   opt.method = dots[["opt.method"]],
                   estimand = temp.est,
                   ...)
      args <- args[!duplicated(names(args))]
      argn <- lapply(names(args), as.name)
      names(argn) <- names(args)
      f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
      dots[["rkhs.args"]] <- eval(f.call, args)
    }
    if (est == "cATE") {
      if(is.list(constraint)) {
        if(length(constraint) == 2 && !any(names(constraint) %in% c("joint","penalty","margins"))) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_wass_const(x=x, z=1-z, K=constraint0, 
                         p=dots[["p"]], dist = dots[["metric"]], 
                         cost = dots[["cost"]], 
                         add.margins = dots[["add.margins"]],
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.joint = dots[["add.joint"]],
                         rkhs.args = dots[["rkhs.args"]],  
                         neg.weights = dots[["neg.weights"]],
                         bf = bf,
                         sample_weight = sw0,
                         soc = soc),
           qp_wass_const(x=x, z=z, K=constraint1, 
                         p=dots[["p"]], 
                         dist = dots[["metric"]], 
                         cost = dots[["cost"]],
                         add.margins = dots[["add.margins"]],
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.joint = dots[["add.joint"]],
                         rkhs.args = dots[["rkhs.args"]],  
                         neg.weights = dots[["neg.weights"]],
                         bf = bf,
                         sample_weight = sw1,
                         soc = soc))
    } else if (est == "ATE") {
      bf0 <- bf1 <- bf
      if (!is.null(bf)) {
        bf0$mm <- rbind(bf$mm[z == 0, ,drop = FALSE], bf$mm)
        bf1$mm <- rbind(bf$mm[z == 1, ,drop = FALSE], bf$mm)
      }
      if(is.list(constraint)) {
        if(length(constraint) == 2) {
          constraint0 <- constraint[[1]]
          constraint1 <- constraint[[2]]
        } else {
          constraint0 <- constraint
          constraint1 <- constraint
        }
      } else {
        constraint0 <- constraint1 <- constraint
      }
      list(qp_wass_const(x = rbind(x[z == 0, , drop = FALSE],x), 
                         z = c(rep(0, sum(1 - z )), rep(1, nrow(x))),
                         K = constraint0,
                         p = dots[["p"]], dist = dots[["metric"]], 
                         cost = dots[["cost"]][[1]],
                         rkhs.args = dots[["rkhs.args"]], 
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.margins = dots[["add.margins"]], 
                         neg.weights = dots[["neg.weights"]],
                         bf = bf0,
                         sample_weight = sw0,
                         soc = soc),
           qp_wass_const(x = rbind(x[z == 1, ,drop = FALSE],x), 
                         z = c(rep(0, sum(z)), rep(1, nrow(x))),
                         K = constraint1,
                         p = dots[["p"]], dist = dots[["metric"]], 
                         cost = dots[["cost"]][[2]],
                         rkhs.args = dots[["rkhs.args"]], 
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.margins = dots[["add.margins"]], 
                         neg.weights = dots[["neg.weights"]],
                         bf = bf1,
                         sample_weight = sw1,
                         soc = soc))
    } else if (est == "ATC") {
      list(qp_wass_const(x=x, z=1-z, K=constraint, 
                         p=dots[["p"]], dist = dots[["metric"]], 
                         cost = dots[["cost"]], 
                         add.margins = dots[["add.margins"]],
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.joint = dots[["add.joint"]],
                         rkhs.args = dots[["rkhs.args"]],  
                         neg.weights = dots[["neg.weights"]],
                         bf = bf,
                         sample_weight = sw,
                         soc = soc))
    } else {
      list(qp_wass_const(x=x, z=z, K=constraint, 
                         p=dots[["p"]], 
                         dist = dots[["metric"]], cost = dots[["cost"]],
                         penalty = dots[["penalty"]],
                         joint.mapping = dots[["joint.mapping"]],
                         add.margins = dots[["add.margins"]],
                         add.joint = dots[["add.joint"]],
                         rkhs.args = dots[["rkhs.args"]], 
                         neg.weights = dots[["neg.weights"]],
                         bf = bf,
                         sample_weight = sw,
                         soc = soc))
    }
  } else if (meth == "RKHS" | meth == "RKHS.dose") {
    dots <- list(...)
    if(is.null(dots[["p"]])) dots[["p"]] <- 2
    if(is.null(dots$metric)) dots$metric <- "mahalanobis"
    if(is.null(dots$theta)) dots$theta <- c(1,1)
    if(is.null(dots$gamma)) dots$gamma <- c(1,1)
    if(is.null(dots$lambda)) dots$lambda <- 0
    if(is.null(dots$sigma_2)) dots$sigma_2 <- 0
    
    if(est == "cATE") {
      list(qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$metric, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = "ATT"),
           qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$metric, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = "ATT"))
    } else {
      list(qp_rkhs(x=x, z=z,
                   p=dots$p, theta = dots$theta, gamma = dots$gamma,
                   lambda = dots$lambda, sigma_2 = dots$sigma_2,
                   dist = dots$metric, cost = dots$cost,
                   is.dose = isTRUE(meth == "RKHS.dose"),
                   estimand = estimand))
    }
  }
  return(qp)
}

quadprog.DataSim <- function(data, constraint,  estimand = c("ATT", "ATC", "ATE","cATE","feasible"), 
                             method = supported.methods(),
                             sample_weight = NULL,
                             soc = FALSE,
                             ...) {
  # meth <- match.arg(method, supported.methods())
  # est <- match.arg(estimand)
  # 
  x <- data$get_x()
  z <- data$get_z()
  y <- data$get_y()
  
  return(quadprog.default(x = x, z = z, y = y, 
                          constraint = constraint,  
                          estimand = estimand, 
                          method = method,
                          sample_weight = NULL,
                          soc = soc,
                          ...))
  
  # dots <- list(...)
  # 
  # form <- dots$formula
  # if ( isTRUE(!is.null(form)) & isTRUE(!is.na(form)) ) {
  #   form <- form_all_squares(form, colnames(data$get_x()))
  #   
  #   # if (is.character(form)) {
  #   #   form.temp <- strsplit(form, "~")[[1]][2]
  #   # 
  #   # } else if (inherits(form,"formula")) {
  #   #   form.temp <- as.character(form[3])
  #   # }
  #   # form.terms <- strsplit(form.temp, "\\+")[[1]]
  #   # is.square  <- grepl("I\\(\\s*\\.\\^2\\s*\\)", form.terms)
  #   # form.terms <- form.terms[!is.square]
  #   # form.nsq   <- paste0(form.terms, collapse = "+")
  #   # square.terms <- NULL
  #   # if ( any(is.square) ) {
  #   #   square.terms <- paste0("I(",colnames(data$get_x()), "^2)", collapse = " + ")
  #   # }
  #   # form <- as.formula(paste0("~ 0 + ",
  #   #                           paste0(c(form.nsq, square.terms), 
  #   #                                  collapse = " + "))
  #   #                    )
  #   form.temp <- as.character(form[length(form)])
  #   form <- as.formula(paste0("~ 0 +", form.temp))
  #   mm <- model.matrix(form, data = data.frame(data$get_x()))
  #   if ( all(mm[,1] == 1)) mm <- mm[,-1]
  #   if (method == "Wasserstein" | method == "Constrained Wasserstein") {
  #     bf <- list(mm = mm,
  #                K = dots$balance.constraints)
  #   }
  # } else {
  #   if (method == "SBW") {
  #     form <- formula(~ . + 0)
  #     mm <- model.matrix(form, data = data.frame(data$get_x()))
  #   } else {
  #     bf <- NULL
  #   }
  # }
  # # margmass <- get_sample_weight(sample_weight, z = data$get_z())
  # 
  # qp <- if (meth == "SBW") {
  #   # if (length(constraint) != ncol(mm))  {
  #   #   K <- rep(constraint, ncol(mm))[1:ncol(mm)]
  #   # } else {
  #     K <- constraint
  #   # }
  #   if (est == "cATE") {
  #     list(qp_sbw(x = mm, z = data$get_z(), K = K, estimand = "ATC"),
  #          qp_sbw(x = mm, z = data$get_z(), K = K, estimand = "ATT"))
  #   } else if (est == "ATE") {
  #     qp_sbw(x = mm, z = data$get_z(), K = K, estimand = est)
  #   } else {
  #     list(qp_sbw(x = mm, z = data$get_z(), K = K, estimand = est))
  #   } 
  # } else if (meth == "Wasserstein") {
  #   
  #   if (is.null(dots$p)) dots$p <- 2
  #   if (is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if (is.null(dots$add.margins)) dots$add.margins <- FALSE
  #   
  #   if (dots$metric == "RKHS" & is.null(dots$rkhs.args)) {
  #     if (is.null(dots$opt.method)) dots$opt.method <- "stan"
  #     temp.est <- switch(estimand,
  #                        "ATT" = "ATT",
  #                        "ATC" = "ATC",
  #                        "ATE"
  #                        )
  #     if (is.null(dots$kernel)) dots$kernel <- "RBF"
  #     args <- list(x = data$get_x(), 
  #                  z = data$get_z(), 
  #                  y = data$get_y(),
  #                  metric = switch(dots$metric,
  #                                  "RKHS" = "mahalanobis",
  #                                  "mahalanobis" = "mahalanobis",
  #                                  "Lp" = "Lp"),
  #                  kernel = dots$kernel,
  #                  is.dose = dots$is.dose,
  #                  opt.method = dots$opt.method,
  #                  estimand = temp.est,
  #                  ...)
  #     args <- args[!duplicated(names(args))]
  #     argn <- lapply(names(args), as.name)
  #     names(argn) <- names(args)
  #     f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
  #     dots$rkhs.args <- eval(f.call, args)
  #   }
  #   if (est == "cATE") {
  #     list(qp_wass(x = data$get_x(), z=data$get_z(), K = constraint,
  #                  p=dots$p, estimand = "ATC", dist = dots$metric, cost = dots$cost,
  #                  rkhs.args = dots$rkhs.args, 
  #                  add.margins = dots$add.margins, bf = bf,
  #                  sample_weight = sample_weight),
  #          qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
  #                  p=dots$p, estimand = "ATT", dist = dots$metric, cost = dots$cost,
  #                  rkhs.args = dots$rkhs.args, 
  #                  add.margins = dots$add.margins, bf = bf,
  #                  sample_weight = sample_weight))
  #   } else if (est == "ATE") {
  #     qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
  #             p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #             rkhs.args = dots$rkhs.args, 
  #             add.margins = dots$add.margins, bf = bf,
  #             sample_weight = sample_weight)
  #   } else {
  #     list(qp_wass(x=data$get_x(), z=data$get_z(), K = constraint,
  #                  p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #                  rkhs.args = dots$rkhs.args, 
  #                  add.margins = dots$add.margins, bf = bf,
  #                  sample_weight = sample_weight))
  #   }
  #   
  # } else if (meth == "Constrained Wasserstein") {
  #   dots <- list(...)
  #   if(is.null(dots$p)) dots$p <- 2
  #   if(is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if (is.null(dots$add.joint)) dots$add.joint <- TRUE
  #   if (is.null(dots$add.margins)) dots$add.margins <- FALSE
  #   if(!dots$add.joint & !dots$add.margins) stop("must run margins or joint")
  #   if(dots$metric == "RKHS" & is.null(dots$rkhs.args)) {
  #     if(is.null(dots$opt.method)) dots$opt.method <- "stan"
  #     temp.est <- switch(estimand,
  #                        "ATT" = "ATT",
  #                        "ATC" = "ATC",
  #                        "ATE"
  #     )
  #     if(is.null(dots$kernel)) dots$kernel <- "RBF"
  #     args <- list(x=data$get_x(), 
  #                  z=data$get_z(), 
  #                  y = data$get_y(),
  #                  metric = switch(dots$metric,
  #                                  "RKHS" = "mahalanobis",
  #                                  "mahalanobis" = "mahalanobis",
  #                                  "Lp" = "Lp"),
  #                  kernel = dots$kernel,
  #                  is.dose = dots$is.dose,
  #                  opt.method = dots$opt.method,
  #                  estimand = temp.est,
  #                  ...)
  #     args <- args[!duplicated(names(args))]
  #     argn <- lapply(names(args), as.name)
  #     names(argn) <- names(args)
  #     f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
  #     dots$rkhs.args <- eval(f.call, args)
  #   }
  #   if (est == "cATE") {
  #     list(qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
  #                  p=dots$p, estimand = "ATC", dist = dots$metric, 
  #                  cost = dots$cost, 
  #                  add.margins = dots$add.margins,
  #                  add.joint = dots$add.joint,
  #                  rkhs.args = dots$rkhs.args,  bf = bf,
  #                  sample_weight = sample_weight),
  #          qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
  #                  p=dots$p, estimand = "ATT", 
  #                  dist = dots$metric, cost = dots$cost,
  #                  add.margins = dots$add.margins,
  #                  add.joint = dots$add.joint,
  #                  rkhs.args = dots$rkhs.args,  bf = bf,
  #                  sample_weight = sample_weight))
  #   } else if (est == "ATE") {
  #     qp_wass_const(x=data$get_x(), z=data$get_z(),K=constraint, 
  #             p=dots$p, estimand = est, 
  #             dist = dots$metric, cost = dots$cost,
  #             add.margins = dots$add.margins,
  #             add.joint = dots$add.joint,
  #             rkhs.args = dots$rkhs.args,  bf = bf,
  #             sample_weight = sample_weight)
  #   } else {
  #     list(qp_wass_const(x=data$get_x(), z=data$get_z(), K=constraint, 
  #                   p=dots$p, estimand = est, 
  #                   dist = dots$metric, cost = dots$cost,
  #                   add.margins = dots$add.margins,
  #                   add.joint = dots$add.joint,
  #                   rkhs.args = dots$rkhs.args, bf = bf,
  #                   sample_weight = sample_weight))
  #   }
  # } else if (meth == "RKHS" | meth == "RKHS.dose") {
  #   dots <- list(...)
  #   if(is.null(dots$p)) dots$p <- 2
  #   if(is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if(is.null(dots$theta)) dots$theta <- c(1,1)
  #   if(is.null(dots$gamma)) dots$gamma <- c(1,1)
  #   if(is.null(dots$lambda)) dots$lambda <- 0
  #   if(is.null(dots$sigma_2)) dots$sigma_2 <- 0
  #   
  #   if(est == "cATE") {
  #     list(qp_rkhs(x=data$get_x(), z=data$get_z(),
  #           p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #           lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #           dist = dots$metric, cost = dots$cost,
  #           is.dose = isTRUE(meth == "RKHS.dose"),
  #           estimand = "ATT"),
  #          qp_rkhs(x=data$get_x(), z=data$get_z(),
  #                  p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #                  lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #                  dist = dots$metric, cost = dots$cost,
  #                  is.dose = isTRUE(meth == "RKHS.dose"),
  #                  estimand = "ATT"))
  #   } else {
  #     list(qp_rkhs(x=data$get_x(), z=data$get_z(),
  #             p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #             lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #             dist = dots$metric, cost = dots$cost,
  #             is.dose = isTRUE(meth == "RKHS.dose"),
  #             estimand = estimand))
  #   }
  # }
  # return(qp)
}

quadprog.data.frame <- function(data, constraint,  
                                estimand = c("ATT", "ATC", "ATE", "cATE", "feasible"), 
                                method = supported.methods(),
                                sample_weight = NULL,
                                soc = FALSE,
                                ...) {
  meth <- match.arg(method)
  est <- match.arg(estimand)
  soc <- isTRUE(soc)
  
  df <- prep_data(data, ...)
  z <- as.numeric(df$z)
  x.df <- df$df[,!(colnames(df$df) == "y")]
  x <- as.matrix(x.df)
  cn <- colnames(x.df)
  x.df <- as.data.frame(x.df)
  colnames(x) <- colnames(x.df) <- cn
  y <- df$df[,(colnames(df$df) == "y")]
  
  return(quadprog.default(x = x, z = z, y = y, 
                          constraint = constraint,  
                          estimand = estimand, 
                          method = method,
                          sample_weight = NULL,
                          soc = soc,
                          ...))
  # col_x <- ncol(x)
  
  # dots <- list(...)
  # form <- dots$formula
  # 
  # if (isTRUE(!is.null(form)) & isTRUE(!is.na(form))) {
  #   form <- form_all_squares(form, colnames(x.df))
  #   # if (is.character(form)) {
  #   #   form.temp <- strsplit(form, "~")[[1]][2]
  #   #   
  #   # } else if (inherits(form,"formula")) {
  #   #   form.temp <- as.character(form[3])
  #   # }
  #   # form.terms <- strsplit(form.temp, "\\+")[[1]]
  #   # is.square  <- grepl("I\\(\\s*\\.\\^2\\s*\\)", form.terms)
  #   # form.terms <- form.terms[!is.square]
  #   # form.nsq   <- paste0(form.terms, collapse = "+")
  #   # square.terms <- NULL
  #   # if ( any(is.square) ) {
  #   #   square.terms <- paste0("I(",colnames(x.df), "^2)", collapse = " + ")
  #   # }
  #   # form <- as.formula(paste0("~ 0 + ",
  #   #                           paste0(c(form.nsq, square.terms), 
  #   #                                  collapse = " + "))
  #   # )
  #   
  #   form.temp <- as.character(form[length(form)])
  #   form <- as.formula(paste0("~ 0 +", form.temp))
  #   mm <- model.matrix(form, data = x.df)
  #   if ( all(mm[,1] == 1)) mm <- mm[,-1]
  #   if (method == "Wasserstein" | method == "Constrained Wasserstein") {
  #     bf <- list(mm = mm,
  #                K = dots$balance.constraints)
  #   }
  # } else {
  #   if (method == "SBW") {
  #     form <- formula(~ . + 0)
  #     mm <- model.matrix(form, data = x.df)
  #   } else {
  #     mm <- NULL
  #     bf <- NULL
  #   }
  # }
  # margmass <- get_sample_weight(sample_weight, z)
  # 
  # qp <- if (meth == "SBW") {
  #   if (length(constraint) != ncol(mm))  {
  #     K <- rep(constraint, ncol(mm))[1:ncol(mm)]
  #   } else {
  #     K <- constraint
  #   }
  #   if (est == "cATE") {
  #     list(qp_sbw(x = mm, z = z, K = K, estimand = "ATC"),
  #          qp_sbw(x = mm, z = z, K = K, estimand = "ATT"))
  #   } else if (est == "ATE") {
  #     qp_sbw(x = mm, z = z, K = K, estimand = est)
  #   } else {
  #     list(qp_sbw(x = mm, z = z, K = K, estimand = est))
  #   } 
  # } else if (meth == "Wasserstein") {
  #   # dots <- list(...)
  #   if (is.null(dots$p)) dots$p <- 2
  #   if (is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if (is.null(dots$add.margins)) dots$add.margins <- FALSE
  #   if (dots$metric == "RKHS" & is.null(dots$rkhs.args)) {
  #     if (is.null(dots$opt.method)) dots$opt.method <- "stan"
  #     temp.est <- switch(estimand,
  #                        "ATT" = "ATT",
  #                        "ATC" = "ATC",
  #                        "ATE"
  #     )
  #     if (is.null(dots$kernel)) dots$kernel <- "RBF"
  #     args <- list(x = x, 
  #                  z = z, 
  #                  y = y,
  #                  metric = switch(dots$metric,
  #                                  "RKHS" = "mahalanobis",
  #                                  "mahalanobis" = "mahalanobis",
  #                                  "Lp" = "Lp"),
  #                  kernel = dots$kernel,
  #                  is.dose = dots$is.dose,
  #                  opt.method = dots$opt.method,
  #                  estimand = temp.est,
  #                  ...)
  #     args <- args[!duplicated(names(args))]
  #     argn <- lapply(names(args), as.name)
  #     names(argn) <- names(args)
  #     f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
  #     dots$rkhs.args <- eval(f.call, args)
  #   }
  #   if(est == "cATE") {
  #     list(qp_wass(x=x, z=z, K = constraint,
  #                  p=dots$p, estimand = "ATC", dist = dots$metric, cost = dots$cost,
  #                  add.margins = dots$add.margins,
  #                  bf = bf,
  #                  sample_weight = sample_weight),
  #          qp_wass(x=x, z=z, K = constraint,
  #                  p=dots$p, estimand = "ATT", dist = dots$metric, cost = dots$cost,
  #                  add.margins = dots$add.margins,
  #                  bf = bf,
  #                  sample_weight = sample_weight))
  #   } else if (est == "ATE") {
  #     qp_wass(x=x, z=z, K = constraint,
  #             p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #             add.margins = dots$add.margins,
  #             bf = bf,
  #             sample_weight = sample_weight)
  #   } else {
  #     list(qp_wass(x=x, z=z, K = constraint,
  #                  p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #                  add.margins = dots$add.margins,
  #                  bf = bf,
  #                  sample_weight = sample_weight))
  #   }
  #   
  # } else if (meth == "Constrained Wasserstein") {
  #   dots <- list(...)
  #   if(is.null(dots$p)) dots$p <- 2
  #   if(is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if (is.null(dots$add.joint)) dots$add.joint <- TRUE
  #   if (is.null(dots$add.margins)) dots$add.margins <- FALSE
  #   if(!dots$add.joint & !dots$add.margins) stop("must run margins or joint")
  #   if(dots$metric == "RKHS" & is.null(dots$rkhs.args)) {
  #     if(is.null(dots$opt.method)) dots$opt.method <- "stan"
  #     temp.est <- switch(estimand,
  #                        "ATT" = "ATT",
  #                        "ATC" = "ATC",
  #                        "ATE"
  #     )
  #     if(is.null(dots$kernel)) dots$kernel <- "RBF"
  #     args <- list(x = x, 
  #                  z = z, 
  #                  y = y,
  #                  metric = switch(dots$metric,
  #                                  "RKHS" = "mahalanobis",
  #                                  "mahalanobis" = "mahalanobis",
  #                                  "Lp" = "Lp"),
  #                  kernel = dots$kernel,
  #                  is.dose = dots$is.dose,
  #                  opt.method = dots$opt.method,
  #                  estimand = temp.est,
  #                  ...)
  #     args <- args[!duplicated(names(args))]
  #     argn <- lapply(names(args), as.name)
  #     names(argn) <- names(args)
  #     f.call <- as.call(c(list(as.name("RKHS_param_opt")), argn))
  #     dots$rkhs.args <- eval(f.call, args)
  #   }
  #   if (est == "cATE") {
  #     list(qp_wass_const(x=x, z=z,K=constraint, 
  #                  p=dots$p, estimand = "ATC", dist = dots$metric, cost = dots$cost,
  #                  bf = bf,
  #                  add.margins = dots$add.margins,
  #                  add.joint = dots$add.joint,
  #                  sample_weight = sample_weight),
  #          qp_wass_const(x=x, z=z,K=constraint, 
  #                  p=dots$p, estimand = "ATT", dist = dots$metric, cost = dots$cost,
  #                  bf = bf,
  #                  add.margins = dots$add.margins,
  #                  add.joint = dots$add.joint,
  #                  sample_weight = sample_weight))
  #   } else if (est == "ATE") {
  #     qp_wass_const(x=x, z=z,K=constraint, 
  #             p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #             bf = bf,
  #             add.margins = dots$add.margins,
  #             add.joint = dots$add.joint,
  #             sample_weight = sample_weight)
  #   } else {
  #     list(qp_wass_const(x=x, z=z, K=constraint, 
  #                        p=dots$p, estimand = est, dist = dots$metric, cost = dots$cost,
  #                        bf = bf,
  #                        add.margins = dots$add.margins,
  #                        add.joint = dots$add.joint,
  #                        sample_weight = sample_weight))
  #   }
  # } else if (meth == "RKHS" | meth == "RKHS.dose") {
  #   dots <- list(...)
  #   if(is.null(dots$p)) dots$p <- 2
  #   if(is.null(dots$metric)) dots$metric <- "mahalanobis"
  #   if(is.null(dots$theta)) dots$theta <- c(1,1)
  #   if(is.null(dots$gamma)) dots$gamma <- c(1,1)
  #   if(is.null(dots$lambda)) dots$lambda <- 0
  #   if(is.null(dots$sigma_2)) dots$sigma_2 <- 0
  #   
  #   if(est == "cATE") {
  #     list(qp_rkhs(x=x, z=z,
  #                  p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #                  lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #                  dist = dots$metric, cost = dots$cost,
  #                  is.dose = isTRUE(meth == "RKHS.dose"),
  #                  estimand = "ATC"),
  #          qp_rkhs(x=x, z=z,
  #                  p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #                  lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #                  dist = dots$metric, cost = dots$cost,
  #                  is.dose = isTRUE(meth == "RKHS.dose"),
  #                  estimand = "ATT"))
  #   } else {
  #     list(qp_rkhs(x=x, z=z,
  #                  p=dots$p, theta = dots$theta, gamma = dots$gamma,
  #                  lambda = dots$lambda, sigma_2 = dots$sigma_2,
  #                  dist = dots$metric, cost = dots$cost,
  #                  is.dose = isTRUE(meth == "RKHS.dose"),
  #                  estimand = est))
  #   }
  # }
  # return(qp)
}


qp_sbw <- function(x, z, K, estimand = c("ATT", "ATC", 
                                         "ATE", "feasible")) {
  est <- match.arg(estimand)
  if (est == "ATC") {
    z <- 1 - z
  }
  
  x1 <- x[z == 1,,drop = FALSE]
  x0 <- x[z == 0,,drop = FALSE]
  x1_var <- colVar(x1)
  x0_var <- colVar(x0)
  
  if (est != "ATE") {
    pool_sd <- sqrt((x1_var + x0_var) / 2)
    
    if (length(K) != 1 & length(K) != ncol(x)) stop("Constraint must be length 1 or length of number of columns for estimand ", estimand, " in SBW")
    
    K_sd <- K * pool_sd
    K_upr <- K_sd
    K_lwr <- (-K_sd)
  } else {
    x_var <- colVar(x)
    xm <- colMeans(x)
    pool_sd_0 <- sqrt((x_var + x0_var) / 2)
    pool_sd_1 <- sqrt((x1_var + x_var) / 2)
    
    if (length(K) != 2) {
      if (is.numeric(K) & length(K) == 1) {
        K <- c(K, K)
      } else if (is.numeric(K) & length(K) == ncol(x)) {
        K <- list(K, K)
      } else {
        stop("Constraint must be a length 1 or 2 numeric vector or a length 2 list of numeric vectors or appropriate length ",est, " in SBW")
      }
    }
    
    K_sd_0 <- K[[1]] * pool_sd_0
    K_sd_1 <- K[[2]] * pool_sd_1
    
    K_upr_0  <- K_sd_0 + xm
    K_lwr_0  <- (-K_sd_0 + xm)
    K_upr_1  <- K_sd_1 + xm
    K_lwr_1  <- (-K_sd_1 + xm)
  }
  
  if (est == "ATE") {
    n <- nrow(x)
    n1 <- nrow(x1)
    n0 <- nrow(x0)
    
    d <- ncol(x)
    Q0_0 <- Matrix::Diagonal(n0, 1) 
    Q0_1 <- Matrix::Diagonal(n1, 1)
    
    A1_0 <- Matrix::sparseMatrix(i = rep(1, n0),
                                 j = 1:n0,
                                 x = 1,
                                 dims = c(1, n0))
    A1_1 <- Matrix::sparseMatrix(i = rep(1, n1),
                                 j = 1:n1,
                                 x = 1,
                                 dims = c(1, n1))
    
    x_constraint_0 <- t(x0)
    x_constraint_1 <- t(x1)
    
    
    L0_0 <- c(w = rep(0, n0))
    L0_1 <- c(w = rep(0, n1))
    
    
    A2_0 <- x_constraint_0
    A2_1 <- x_constraint_1
    
    
    A_1 <- rbind(A1_1, A2_1)
    lc_1 <- c(rep(1, nrow(A1_1)), K_lwr_1)
    uc_1 <- c(rep(1, nrow(A1_1)), K_upr_1)
    
    
    A_0 <- rbind(A1_0, A2_0)
    lc_0 <- c(rep(1, nrow(A1_0)), K_lwr_0)
    uc_0 <- c(rep(1, nrow(A1_0)), K_upr_0)
    
    # A_1 <- rbind(A1_1, A2_1, A2_1)
    # A_0 <- rbind(A1_0, A2_0, A2_0)
    # 
    # vals_0 <- c(rep(1, nrow(A1_0)), K_lwr_0, K_upr_0)
    # vals_1 <- c(rep(1, nrow(A1_1)), K_lwr_1, K_upr_1)
    # 
    # dir_0 <- c(rep("E", nrow(A1_0)),rep("G",length(K_lwr_0)), rep("L",length(K_upr_0)))
    # dir_1 <- c(rep("E", nrow(A1_1)),rep("G",length(K_lwr_1)), rep("L",length(K_upr_1)))
    # 
    program_0 <- list(obj = list(Q = Q0_0, L = L0_0),
                      LC = list(A = A_0, 
                                lc = lc_0,
                                uc = uc_0
                                # dir = dir_0,
                                # vals = vals_0
                                ),
                      nvar = n0)
    
    program_1 <- list(obj = list(Q = Q0_1, L = L0_1),
                      LC = list(A = A_1, 
                                lc = lc_1,
                                uc = uc_1
                                # dir = dir_1,
                                # vals = vals_1
                                ),
                      nvar = n1)
    
    program_1 <- add_bounds(program_1, FALSE)
    program_0 <- add_bounds(program_0, FALSE)
    
    return(list(program_0, program_1))
    
  } else if (est != "feasible") {
    n <- nrow(x0)
    d <- ncol(x0)
    x_constraint <- t(x0)
    x1m <- colMeans(x1)
    K_upr <- K_upr + x1m
    K_lwr <- K_lwr + x1m
    Q0 <- Matrix::Diagonal(n, 1) 
    A1 <- matrix(1,nrow = 1, ncol = n)
  } else if (est == "feasible") {
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
                               dims = c(n,n), repr = "T"
    )
    A1 <- Matrix::sparseMatrix(i = c(rep.int(1,n0), rep.int(2,n1)),
                               j = c(1:n0, n0 + 1:n1),
                               x = 1,
                               dims = c(2,n))
  }
  
  
  L0 <- c(w = rep(0,n))
  
  A2 <- x_constraint
  
  A <- rbind(A1, A2)
  
  lc <- c(rep(1, nrow(A1)), K_lwr)
  uc <- c(rep(1, nrow(A1)), K_upr)
  
  program <- list(obj = list(Q = Q0, L = L0),
                  LC = list(A = A, lc = lc, uc = uc))
  program$nvar <- n
  # add bounds
  program <- add_bounds(program, FALSE)
  return(program)
}


set_K <- function(K, x1, x0, d_c, add.joint, add.margins, penalty, joint.mapping,
                  method) {
  lc0 <- function(x) length(x) == 0
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  n <- n0 + n1
  d <- ncol(x0)
  # list(penalty = NULL,
  #      joint   = NULL,
  #      margins = NULL,
  #      balance = NULL)
  K_joint <- K_margins <- K_penalty <- NULL
  if (is.list(K)) {
    if(add.joint || add.margins || penalty != "none") {
      if (!any(names(K) %in% c("joint", "margins","penalty") ) ) {
        stop("constraint should be a list with optional named slots `joint`, `margins`, and/or `penalty`. Specify the one needed for goal.",)
      }
    }
    K_joint <- as.numeric(K$joint)
    K_margins <- as.numeric(K$margins)
    K_penalty <- as.numeric(K$penalty)
    
  } else if (length(K) == 1 & is.numeric(K) & method == "Constrained Wasserstein") {
    K_joint <- K
  } else if (length(K) == 1 & is.numeric(K) & method == "Wasserstein") {
    K_penalty <- K
  } else if (length(K) == (d_c + 1) & is.numeric(K) & method == "Constrained Wasserstein") {
    K_joint <- K[d_c + 1]
    K_margins <- K[1:d_c]
  } else if (length(K) == (d_c + 1) & is.numeric(K) & method == "Wasserstein") {
    K_penalty <- K[d_c + 1]
    K_margins <- K[1:d_c]
  } else if (add.joint || add.margins || penalty != "none") {
    stop("constraint should be a list with optional named slots `joint`, `margins`, `penalty`. Specify the one needed for goal.",)
  }
  
  if (!add.joint) K_joint <- NULL
  if (!add.margins) K_margins <- NULL
  
  if (add.margins | add.joint) {
    if ( (lc0(K_joint) & add.joint) | (lc0(K_margins) & add.margins)) {
      sigma1 <- cov(x1)
      sigma0 <- cov(x0)
      v1 <- diag(sigma1)
      v0 <- diag(sigma0)
      shared_var <- n1/n * v1 + n0/n * v0
      
      if (add.margins & lc0(K_joint)) {
        K_margins <- sqrt(0.2^2 * shared_var * 2  + v1 + v0 - 2 * sqrt(sqrt(v1) * v0 * sqrt(v1)))
      } else if (add.joint & !add.margins) {
        K_joint <- sum(0.2^2 * shared_var * 2) + sum(v1) + sum(v0) -
          2 * sum(diag(sqrt_mat(sqrt_mat(sigma1) %*% sigma0 %*% sqrt_mat(sigma1))))
      } else if (add.margins & !lc0(K_joint)) {
        K_margins <- sqrt(K_joint^2/d)
      }
    }
    if (add.joint) {
      if (length(K_joint) == d_c + 1 & (add.margins &  lc0(K_margins))) {
        if (add.margins & lc0(K_margins)) {
          K_margins <- K_joint[1:d_c]
          K_joint <- K_joint[d_c + 1]
        } else {
          K_margins <- NULL
          K_joint <- K_joint[d_c + 1]
        }
      } else if (length(K_joint) != 1) {
        stop("constraint$joint must be length 1 ")
      }
    }
    if (add.margins) {
      if(length(K_margins) != d_c){
        K_margins <- rep(K_margins, d_c)[1:d_c]
      }
    } 
  }
  
  if (method == "Constrained Wasserstein" & length(K_penalty) == 0) {
    K_penalty <- 1
  } else if (method == "Wasserstein") {
    if (penalty != "none" && (length(K_penalty) != 1 || !is.numeric(K_penalty)) ) {
      stop("slot `penalty` in constraint list has an error.")
    }
    if (joint.mapping) {
      if (length(K_joint) != 1 || !is.numeric(K_joint))  {
        stop("Must specify slot `joint` in constraint list if using joint mapping")
      }
    }
  }
  
  K <- list(joint = K_joint,
            margins = K_margins,
            penalty = K_penalty)
  
  return(K)
}


add_bc <- function(op, bf, z, K) {
  # handle mm
  if (!is.null(bf)) {
    if (!is.list(bf)) stop("Must specify balance contraints for balance functions.
                            calc_weight(data, constraints = ...,
                            balance.constraints = ..., formula = ...)")
    mm  <- bf$mm
    
    mm0 <- mm[z == 0,, drop = FALSE]
    mm1 <- mm[z == 1,, drop = FALSE ]
    
    n0 <- nrow(mm0)
    n1 <- nrow(mm1)
    n <- n0 + n1
    
    mmse  <- sqrt(matrixStats::colVars(mm0) * n0/n +
                    matrixStats::colVars(mm1) * n1/n )
    
    mmbal <-  Matrix::crossprod(mm0, vec_to_row_constraints(n0,n1))
    mmtarg  <- colMeans(mm1)
    
    Kmm <- bf$K * mmse
    Kmm_low <- -Kmm + mmtarg
    Kmm_high <- Kmm + mmtarg
    
    op$LC$A <- rbind(op$LC$A, mmbal)
    
    op$LC$lc <- c(op$LC$lc, sbw = Kmm_low)
    op$LC$uc <- c(op$LC$uc, sbw = Kmm_high)
    
    # op$LC$vals <- c(op$LC$vals, "bc_up" = Kmm_high, "bc_low" = -Kmm_low)
    # op$LC$dir <- c(op$LC$dir, rep("L", 2 * length(Kmm_high)))
    # op$LC$A <- rbind(op$LC$A, mmbal, -mmbal)
  } 
  
  return(op)
  
}

pen_var <- function(qp, n0, n1, lambda) {
  
  nvar <- length(qp$obj$L)
  
  # if (divergence) {
  #   extra_col <- nvar - (n0 * n1) + 3
  #   Fnew <- rbind(
  #     c(rep(0, nvar) , 1),
  #     c(rep(0, nvar + 1)),
  #     cbind(vec_to_row_constraints(n0,n1),
  #           zero_mat_sp(n0, extra_col)
  #     ))
  #   gnew <- c(0, 1, rep(0,n0))
  #   cones_new <- matrix(list("RQUAD",  n0 + 2, NULL),
  #                       nrow = 3, ncol = 1)
  #   
  #   if (is.null(qp$cones)) {
  #     qp$cones <- list()
  #     
  #     
  #     qp$cones$F <- Fnew
  #     qp$cones$g <- gnew
  #     qp$cones$cones <- cones_new
  #     
  #   } else {
  #     
  #     qp$cones$F <- rbind(cbind(qp$cones$F, rep(0, nrow(qp$cones$F))),
  #                         Fnew)
  #     qp$cones$g <- c(qp$cones$g, gnew)
  #     qp$cones$cones <- cbind(qp$cones$cones, 
  #                       cones_new)
  # 
  # 
  #   }
  #   qp$obj$L <- c(qp$obj$L, pen = lambda)
  #   qp$LC$A <- cbind(qp$LC$A, rep(0, nrow(qp$LC$A)))
  #   qp$bounds$lb <- c(qp$bounds$lb, 0)
  #   qp$bounds$ub <- c(qp$bounds$ub, Inf)
  # } else {
    extra_col <- nvar - (n0 * n1) + 1
    Fnew <- rbind(
      c(rep(0, nvar) , 1),
      c(rep(0, nvar + 1)),
      cbind(vec_to_row_constraints(n0,n1),
            zero_mat_sp(n0, extra_col)
      ))
    gnew <- c(0, 1, rep(0,n0))
    cones_new <- matrix(list("RQUAD",  n0 + 2, NULL),
                        nrow = 3, ncol = 1)
    
    if (is.null(qp$cones)) {
      qp$cones <- list()
      
      
      qp$cones$F <- Fnew
      qp$cones$g <- gnew
      qp$cones$cones <- cones_new
      
    } else {
      
      qp$cones$F <- rbind(cbind(qp$cones$F, rep(0, nrow(qp$cones$F))),
                          Fnew)
      qp$cones$g <- c(qp$cones$g, gnew)
      qp$cones$cones <- cbind(qp$cones$cones, 
                        cones_new)


    }
    qp$obj$L <- c(qp$obj$L, pen = lambda)
    qp$LC$A <- cbind(qp$LC$A, rep(0, nrow(qp$LC$A)))
    qp$bounds$lb <- c(qp$bounds$lb, 0)
    qp$bounds$ub <- c(qp$bounds$ub, Inf)
  # }
  return(qp)
}

qp_pen <- function(qp, n0, n1, a, b, penalty, lambda, soc, divergence) {
  nvar <- n0 * n1
  clength <- length(qp$obj$L)
  if (penalty == "L2") {
    if (soc) {
      if (is.null(qp$cones)) {
        qp$cones <- list()
        
        qp$cones$F <- rbind(
          c(rep(0, clength) , 1),
          c(rep(0, clength + 1)),
          cbind(Matrix::Diagonal(n = nvar, x = 1),
                Matrix::sparseMatrix(i = integer(0),
                                     j = integer(0), x = 0,
                                     dims = c(nvar, 
                                              1 + clength - nvar))))
        qp$cones$g <- c(0, 1, rep(0,nvar))
        qp$cones$cones <- matrix(list("RQUAD",  nvar + 2, NULL),
                                 nrow = 3, ncol = 1)
        
      } else {
        qp$cones$F <- rbind(cbind(qp$cones$F, rep(0, nrow(qp$cones$F))),
                            c(rep(0, ncol(qp$cones$F)) , 1),
                            c(rep(0, ncol(qp$cones$F) + 1)),
                            cbind(Matrix::Diagonal(n = nvar, x = 1),
                                  Matrix::sparseMatrix(i = integer(0),
                                                       j = integer(0), x = 0,
                                                       dims = c(nvar, 
                                                                ncol(qp$cones$F) - nvar + 1 ))))
        qp$cones$g <- c(qp$cones$g, c(0,1, rep(0,nvar)))
        qp$cones$cones <- cbind(qp$cones$cones, 
                                matrix(list("RQUAD",  nvar + 2, NULL),
                                       nrow = 3, ncol = 1))
        
        
      }
      qp$obj$L <- c(qp$obj$L, pen = lambda)
      qp$LC$A <- cbind(qp$LC$A, rep(0, nrow(qp$LC$A)))
      qp$bounds$lb <- c(qp$bounds$lb, 0)
      qp$bounds$ub <- c(qp$bounds$ub, Inf)
    } else {
      if (is.null(lambda)) stop("must specify a term in constraints list = penalty, eg contraint = list(penalty = 1)")
      qp$obj$Q <- if (is.null(qp$obj$Q)) {
        Matrix::Diagonal(n = nvar, x = lambda * 0.5 )
      } else {
        qp$obj$Q + Matrix::Diagonal(n = nvar, x = lambda * 0.5 )
      }
    }
    

  } else if (penalty == "entropy") {
    if (is.null(lambda)) stop("must specify a term in constraints list = penalty, eg contraint = list(penalty = 1)")
    
    qp$obj$L <- c(qp$obj$L, pen =  rep(-lambda,nvar))
    
    if (is.null(qp$cones) ) {
      qp$cones <- list()
      qp$cones$F <- Matrix::sparseMatrix(i = c(seq(2, 3 * nvar,by = 3),
                                               seq(3, 3 * nvar,by = 3)),
                                         j = c(1:nvar, (nvar + 1):(2 * nvar)),
                                         x = 1,
                                         dims = c(nvar * 3, 
                                                  clength + nvar))
      qp$cones$g <- rep(c(1, 0, 0), nvar)
      qp$cones$cones <- matrix(list("PEXP", 3, NULL), nrow = 3, ncol = nvar)
      rownames(qp$cones$cones) <- c("type","dim","conepar")
      
      
    } else {
      qp$cones$F <- rbind(cbind(qp$cones$F, 
                                Matrix::sparseMatrix(i = integer(0),
                                                     j = integer(0), x = 0,
                                                     dims = c(nrow(qp$cones$F), 
                                                              nvar))),
      Matrix::sparseMatrix(i = c(seq(2, 3 * nvar,by=3),
                                               seq(3, 3 * nvar,by=3)),
                                         j = c(1:nvar, (nvar + 2) : (2 * nvar + 1)),
                                         x = 1)
      )
      qp$cones$g <- c(qp$cones$g, rep(c(1, 0, 0), nvar))
      qp$cones$cones <- cbind(qp$cones$cones,
                              matrix(list("PEXP", 3, NULL), 
                                     nrow = 3, ncol = nvar))
      
    }
    
    qp$LC$A    <- cbind(qp$LC$A, 
                        Matrix::sparseMatrix(i = integer(0), 
                                             j = integer(0), x = 0, 
                                             dims = c(nrow(qp$LC$A), nvar)))
    
    qp$bounds$lb <- c(qp$bounds$lb, rep(0, nvar))
    qp$bounds$ub <- c(qp$bounds$ub, rep(Inf, nvar))
    
    if (!is.null(qp$obj$Q)) {
      L <- robust_sqrt_mat(as.matrix(qp$obj$Q[1:n0,1:n0]))
      Lmat <- Matrix::kronecker(X = Matrix::Diagonal(n = n1,
                                                     x = 1),
                                Y = L )
      sparse0 <- Matrix::sparseMatrix(i = integer(0),
                                      j = integer(0), x = 0,
                                      dims = c(nvar, nvar + 1))
      sparseL <- rbind(c(rep(0, nvar * 2), 1),
                       c(rep(0, nvar * 2), 0),
                       cbind(Lmat, sparse0)
                       )
      
      qp$cones$F <- rbind(cbind(qp$cones$F, 
                                Matrix::sparseMatrix(i = integer(0),
                                                                 j = integer(0), x = 0,
                                                                 dims = c(nrow(qp$cones$F), 1)))
                          , sparseL)
      qp$cones$g <- c(qp$cones$g, c(0, 0.5, rep(0, nvar)))
      
      qp$cones$cones <- cbind(qp$cones$cones, matrix(list("RQUAD", nvar + 2, NULL),
                                                     nrow = 3, ncol = 1))
      
      # sparse0 <- Matrix::sparseMatrix(i = integer(0),
      #                                 j = integer(0), x = 0,
      #                                 dims = c(nvar, nvar))
      # sparsebottom <- Matrix::sparseMatrix(i = integer(0),
      #                                      j = integer(0), x = 0,
      #                                      dims = c(nvar, 2 * nvar))
      # sparseL <- rbind(cbind(Lmat, sparse0), sparsebottom)
      # q_ineq <- Rmosek::mosek_qptoprob(F = sparseL, f = qp$obj$L,
      #                                  A = qp$LC$A,
      #                                  b = qp$LC$vals,
      #                                  lb = rep(-Inf, nvar * 2),
      #                                  ub = rep(Inf, nvar * 2))
      # qp$cones$cones <- cbind(qp$cones$cones, matrix(c(q_ineq$cones, list(NULL)),
      #                                                nrow = 3, ncol=1))
      # qp$cones$F <- cbind(qp$cones$F, Matrix::sparseMatrix(i = integer(0),
      #                                                      j = integer(0), x = 0,
      #                                                      dims = c(nrow(qp$cones$F), nvar * 2 + 2)))
      # qp$obj$L <- q_ineq$c
      qp$obj$Q <- NULL
      qp$obj$L <- c(qp$obj$L, pen = 1)
      qp$bounds$lb <- c(qp$bounds$lb, 0)
      qp$bounds$ub <- c(qp$bounds$ub, Inf)
      qp$LC$A <- cbind(qp$LC$A, Matrix::sparseMatrix(i = integer(0),
                                                     j = integer(0), x = 0,
                                                     dims = c(nrow(qp$obj$A),1)))
      # qp$LC$vals <- c(qp$LC$vals, rep(0, nvar * 2))
      # qp$LC$dir  <- c(qp$LC$dir, rep("E", nvar * 2))
      # stop("mosek can't handle non-linear problems with conic constraints.")
      
      
      # sparse0 <- Matrix::sparseMatrix(i = integer(0),
      #                                 j = integer(0), x = 0,
      #                                 dims = c(nvar, 2 * nvar))
      # sparsebottom <- Matrix::sparseMatrix(i = integer(0),
      #                                      j = integer(0), x = 0,
      #                                      dims = c(2*nvar, 3*nvar))
      # qp$obj$Q <- rbind(cbind(qp$obj$Q, sparse0), sparsebottom)
      
      # qp$cones$F <- rbind(cbind(qp$conesF, Matrix::sparseMatrix(i = integer(0),
      #                                                     j = integer(0), x = 0,
      #                                                     dims = c(nrow(F), nvar))),
      #                     Matrix::sparseMatrix(i = integer(0),
      #                                          j = integer(0), x = 0,
      #                                          dims = c(nrow(F), ncol(F)))
      # 
      # qp$LC$A    <- cbind(qp$LC$A,
      #                     Matrix::sparseMatrix(i = integer(0),
      #                                          j = integer(0), x = 0,
      #                                          dims = c(nrow(qp$LC$A), nvar)))
    }
    
    
  } else if (penalty == "variance") {
    if (is.null(lambda)) stop("must specify a term in constraints list = penalty, eg contraint = list(penalty = 1)")
    # qmat <- Matrix::kronecker(matrix(lambda, n1, n1), 
    #                           Matrix::Diagonal(n = n0, x = 1))
    # qmat <- Matrix::sparseMatrix(i = rep(1:nvar, each = n1),
    #                              j = rep(sapply(1:n0, function(i) c(seq(i, nvar, n0))),
    #                                      n1),
    #                              x = lambda)
    # qp$obj$Q <- if (is.null(qp$obj$Q)) {
    #   qmat
    # } else {
    #   qp$obj$Q + qmat
    # }
    # qp$obj$L <- if (is.null(qp$obj$L)) {
    #   -2 * rep(a, n1) * lambda
    # } else {
    #   qp$obj$L - 2 * rep(a, n1) * lambda
    # }
    qp <- pen_var(qp, n0, n1, lambda)
  } 
  qp$nvar <- nvar
  return(qp)
}

add_mapping <- function(op, x0, x1, p, b, lambda, penalty) {
  # currently only do p = 2
  # if (penalty == "entropy") {
  #   op$obj$Q <- if (is.null(op$obj$Q) ) {
  #     Matrix::kronecker(X = Matrix::Diagonal(n = nrow(x1),
  #                                            x = 1),
  #                       Y = tcrossprod(x0,x0)) / ncol(x0) / nrow(x0)
  #   } else {
  #     op$obj$Q + Matrix::kronecker(X = Matrix::Diagonal(n = nrow(x0),
  #                                                       x = 1),
  #                                  Y = tcrossprod(x0,x0)) / ncol(x0) / nrow(x0)
  #   }
  # 
  #   op$obj$L <- if (is.null(op$obj$L) ) {
  #     -2 *c(tcrossprod(x0,x1)) / ncol(x0) / nrow(x0)
  #   } else {
  #     op$obj$L * as.numeric(lambda) -2 *c(tcrossprod(x0,x1)) / ncol(x0) / nrow(x0)
  #   }
  # } else {
    d <- ncol(x1)
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    
    nvar <- n0 * n1
    
    # need to scale x0 and x1 so that different objectives are on the same scale
    # and so that the weights can approx x1
    scale_x1 <- sweep(x1, MARGIN = 1, STATS = sqrt(b)/sqrt(d), FUN = "*")
    # scale_x0 <- sweep(x0, MARGIN = 1, STATS = 1.0/(sqrt(b) * sqrt(d)), FUN = "*")
    op$cones <- list()
    
    nvar <- length(op$obj$L)
    op$cones$F <- rbind(c(rep(0, nvar), 1),
                        rep(0, nvar + 1),
                        cbind(Matrix::kronecker(X = Matrix::Diagonal(n = n1, x = 1.0/(sqrt(b) * sqrt(d))),
                                                Y = t(x0)), rep(0, d * n1))
    )
    op$cones$g <- c(0, 0.5, -c(t(scale_x1)))
    op$cones$cones <- matrix(list("RQUAD",  nrow(op$cones$F), NULL),
                             nrow = 3, ncol = 1)
    
    # adjust the cost objective if present and add in extra variable to all components
    op$obj$L <- c(op$obj$L * as.numeric(lambda), mapping = 1)
    op$LC$A <- cbind(op$LC$A, Matrix::sparseMatrix(i = integer(0),
                                                   j = integer(0),
                                                   x = 0,
                                                   dims = c(nrow(op$LC$A),
                                                            1)))
    if (!is.null(op$obj$Q)) {
      op$obj$Q <- Matrix::bdiag(op$obj$Q, 
                                # add in zero matrix to lower right diagonal
                                Matrix::sparseMatrix(i = integer(0),
                                                               j = integer(0),
                                                               x = 0,
                                                               dims = c(1,1)))
    }
    
    op$bounds$lb <- c(op$bounds$lb, 0)
    op$bounds$ub <- c(op$bounds$ub, Inf) 
    
    
    
  return(op)
}

add_bounds <- function(op, neg.weights, add.joint) {
  if (neg.weights) {
    op$bounds <- if (add.joint) {
      list(lb = rep(-Inf, length(op$obj$L)),
                      ub = rep(Inf, length(op$obj$L)))
    } else {
      list(lb = rep(-1, length(op$obj$L)),
           ub = rep(1, length(op$obj$L)))
    }
  } else {
    op$bounds <- list(lb = rep(0, length(op$obj$L)),
                      ub = rep(Inf, length(op$obj$L)))
  }
  return(op)
}

get_marg_const <- function(b, n0,n1, divergence = TRUE) {
  A_marg_const_b_col_joint <- vec_to_col_constraints(n0, n1)
  if (!divergence) {
    return(list(A = A_marg_const_b_col_joint, lc = b, uc = b))
  }
  A_marg_const_b_col <- vec_to_col_constraints(n1, n1)
  A_marg_const_b_row <- vec_to_row_constraints(n1, n1)
  A_marg_const_a_col <- vec_to_col_constraints(n0, n0)
  A_marg_const_a_row <- vec_to_row_constraints(n0, n0)
  A_marg_const_a_row_joint <- vec_to_row_constraints(n0, n1)
  marg_a_combine <- rbind(cbind(A_marg_const_a_row_joint, -A_marg_const_a_row
                                , zero_mat_sp(n0,n1*n1) #comment out
                                ),
                          cbind(A_marg_const_a_row_joint, -A_marg_const_a_col
                                  , zero_mat_sp(n0,n1*n1) # comment out
                                  ))
                          
  marg_b_combine <- rbind(cbind(A_marg_const_b_col_joint, zero_mat_sp(n1, n0*n0)
                                , zero_mat_sp(n1, n1*n1) # comment out
                                )
                          , cbind(zero_mat_sp(n1, n0 * n1), zero_mat_sp(n1, n0*n0), A_marg_const_b_col), #comment out
                          cbind(zero_mat_sp(n1, n0 * n1), zero_mat_sp(n1, n0*n0), A_marg_const_b_row) #comment out
                          )
  return(list(A = rbind(marg_b_combine, marg_a_combine),
       lc = c(b, 
              b, b, #comment out
              rep(0, n0 * 2)),
       uc = c(b, 
              b, b, #comment out
              rep(0, n0 * 2))))         
}

bc_to_dir_const <- function(qp) {
  lc <- qp$LC$lc
  uc <- qp$LC$uc
  A  <- qp$LC$A
  
  
  equal <- which(lc == uc)
  lt    <- which(lc != uc & !is.infinite(uc))
  gt    <- which(lc != uc & !is.infinite(lc))
  
  dir <- c(rep("E", length(equal)), rep("L", length(lt)),
           rep("G", length(gt)))
  vals <- c(lc[equal], uc[lt], lc[gt])
  
  A   <- rbind(A[equal,, drop = FALSE],
               A[lt,, drop = FALSE],
               A[gt,, drop = FALSE])
  
  qp$LC$A <- A
  qp$LC$dir <- dir
  qp$LC$vals <- vals
  return(qp)
}

bc_to_gt_const <- function(qp) { #for quadprog only
  lc <- qp$LC$lc
  uc <- qp$LC$uc
  A  <- qp$LC$A
  
  
  equal_l <- which(lc == uc)
  # equal_u <- which(uc == lc)
  lt    <- which(lc != uc & !is.infinite(uc))
  gt    <- which(lc != uc & !is.infinite(lc))
  
  vals <- c(lc[equal_l], -uc[equal_l], -uc[lt], lc[gt])
  
  A   <- rbind(A[equal_l,, drop = FALSE],
               -A[equal_l,, drop = FALSE],
               -A[lt,, drop = FALSE],
               A[gt,, drop = FALSE])
  
  qp$LC$A <- A
  qp$LC$vals <- vals
  return(qp)
}

bc_to_gt_const_quadprog <- function(qp) { #for quadprog only
  lc <- qp$LC$lc
  uc <- qp$LC$uc
  A  <- qp$LC$A
  
  
  equal_l <- which(lc == uc)
  # equal_u <- which(uc == lc)
  lt    <- which(lc != uc & !is.infinite(uc))
  gt    <- which(lc != uc & !is.infinite(lc))
  
  vals <- c(lc[equal_l], -uc[lt], lc[gt])
  
  A   <- rbind(A[equal_l,, drop = FALSE],
               -A[lt,, drop = FALSE],
               A[gt,, drop = FALSE])
  
  qp$LC$A <- A
  qp$LC$vals <- vals
  qp$LC$dir_eq_num <- length(equal_l)
  return(qp)
}

dir_to_bc <- function(qp) {
  qp$dir
  
  
}

qp_wass <- function(x, z, K = list(penalty = NULL,
                                   joint   = NULL,
                                   margins = NULL,
                                   balance = NULL), 
                    p = NULL,
                    penalty = c("L2","entropy",
                                "variance","none"),
                    dist = dist.metrics(), cost = NULL,
                    rkhs.args = NULL, add.joint = TRUE,
                    add.margins = FALSE,
                    joint.mapping = FALSE,
                    neg.weights = FALSE,
                    bf = NULL,
                    sample_weight = NULL,
                    soc = FALSE, divergence = FALSE) {
  
  # estimand <- match.arg(estimand)
  divergence <- isTRUE(divergence)
  
  # if (divergence) {
  #   return(qp_wass_div(x = x, z = z, K = K, 
  #                      p = p,
  #                      penalty = penalty,
  #                      dist = dist, cost = cost,
  #                      rkhs.args = rkhs.args, add.joint = TRUE,
  #                      add.margins = add.margins,
  #                      joint.mapping = joint.mapping,
  #                      neg.weights = neg.weights,
  #                      bf = bf,
  #                      sample_weight = sample_weight,
  #                      soc = soc))
  # }
  
  margmass = get_sample_weight(sample_weight, z = z)
  joint.mapping <- isTRUE(joint.mapping)
  penalty <- match.arg(penalty)
  neg.weights <- isTRUE(neg.weights)
  soc <- isTRUE(soc)
  
  if (neg.weights && penalty == "entropy") {
    stop("Can't have negative weights with the entropy penalty!")
  }
  
  x1 <- x[z == 1,,drop = FALSE]
  x0 <- x[z == 0,,drop = FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  weight.dim <- n0 * n1
  
  dist <- match.arg(dist)
  
  if(is.null(p)) {
    if(penalty != "entropy") {
      p <- floor(d/2 + 1)
    } else {
      p <- 2
    }
  }
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  if(p < d/2 && penalty != "entropy") warning("power of the wasserstein distance is less than rate optimal d/2 for non-entropy penalties.")
  
  if(is.null(cost) ) {
    if(add.margins) {
      cost <- lapply(1:d, function(i) cost_fun(x[,i,drop = FALSE], z, 
                                               ground_p = p, metric = dist,
                                               rkhs.args = rkhs.args, estimand = "ATT"))
      cost[d + 1] <- list(cost_fun(x, z, 
                                   ground_p = p, metric = dist,
                                   rkhs.args = rkhs.args, estimand = "ATT"))
      
    } else {
      cost <- cost_fun(x, z, 
                       ground_p = p, metric = dist,
                       rkhs.args = rkhs.args, estimand = "ATT")
    }
  }
  
  
  if (add.margins) {
    if (!is.list(cost) && is.matrix(cost)) {
      cost <- c(lapply(1:d, function(i) cost_fun(x[,i,drop = FALSE], z, 
                                                 ground_p = p, metric = dist,
                                                 rkhs.args = rkhs.args, estimand = "ATT")),
                cost)
    }
    d_cost <- length(cost) - 1
    cost.marg <- cost[1:d_cost]
    cost <- cost[[d_cost + 1]]
    
    sapply(cost.marg, function(cc) stopifnot(dim(cc) %in% c(n0,n1)))
    
    marg_cost_vec <- Matrix::sparseMatrix(i = rep(1:d, each = n0 * n1),
                                          j = rep(1:(n0 * n1), d),
                                          x = c(sapply(cost.marg, c)^p),
                                          dims = c(d, n0 * n1))
  } else {
    marg_cost_vec <- cost.marg <- NULL
  }
  
  stopifnot(dim(cost) %in% c(n0,n1))
  
  joint_cost_vec <- Matrix::sparseMatrix(i = rep(1, n0 * n1),
                                         j = 1:(n0*n1),
                                         x = c(cost)^p,
                                         dims = c(1, n0*n1))
  
  K <- set_K(K, x0, x1, d_c = length(cost.marg), joint.mapping, add.margins,
             penalty = penalty, joint.mapping = joint.mapping,
             method = "Wasserstein")
  
  #objective function
  obj <- list(L = c(cost = as.numeric(joint_cost_vec)))
  
  # linear constraints
  LC <- list()
  
  #linear constraint matrix A
  sum_const_A <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                      j = 1:(n0*n1),
                                      x = rep(1, n1 * n0),
                                      dims = c(1, n0*n1))
  marg_const_mat_A <- vec_to_col_constraints(n0,n1)
  
  LC$A <- rbind(sum_const_A, 
                marg_const_mat_A, 
                marg_cost_vec)
  
  # linear constraint values
  sum_const <- 1
  marg_const <- margmass$b 
  K_const  <- c(margins = K$margins)
  
  # set constraint bounds
  LC$uc <- c(sum_1 = 1, sum_b = marg_const, K_const^p)
  LC$lc <- c(sum_1 = 1, sum_b = marg_const, rep(0, length(K_const)))
  
  # if (neg.weights) {
  #   K_const <- c(K_const, joint = 0, margins = rep(0, length(K_const)))
  #   LC$A <- rbind(LC$A, -joint_cost_vec, -1 * marg_cost_vec)
  # }
  # 
  # LC$vals <- c(sum_1 = 1, sum_b = marg_const, K_const^p)
  # 
  # # set direction
  # LC$dir <- c(rep("E", 1 + length(marg_const)),
  #             rep("L", length(K_const)))
  
  # create lp/qp
  op <- list(obj = obj, LC = LC)
  
  # add bal constraints if necessary
  op <- add_bc(op, bf, z, K)
  
  # add bounds
  op <- add_bounds(op, neg.weights, joint.mapping)
  
  # add in mapping
  if (joint.mapping) op <- add_mapping(op, x0, x1, p, margmass$b, K$joint/median(cost), penalty) #currently just p = 2 method
  
  # add in penalty functions
  op <- qp_pen(op, n0, n1, margmass$a, margmass$b, penalty = penalty, K$penalty, soc, FALSE)
  # (qp, n0, n1, a, b, penalty, lambda, soc, divergence)
  
  # K_vals <- rep(K, d + 1)
  # L0 <- rep(0, n0 * n1) #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  # obj <- list(Q = Matrix::Diagonal(weight.dim, 0.5 * K),#Q0,
  #             L = LO)
  # 
  # LC <- list()
  # LC$A <- rbind(sum_const, 
  #               marg_const_mat, 
  #               cost_vec)
  # LC$vals <- c(1, marg_const, K_vals)
  # LC$dir <- c(rep("E", 1 + marg_const_n),
  #             rep("L", length(K_vals)))
  
  # op <- list(obj = obj, LC = LC)
  return(op)
}


qp_wass_div <- function(x, z, K = list(penalty = NULL,
                                   joint   = NULL,
                                   margins = NULL,
                                   balance = NULL), 
                    p = 2,
                    penalty = c("L2","entropy",
                                "variance"),
                    dist = dist.metrics(), cost = list(joint = NULL,
                                                       a = NULL,
                                                       b = NULL,
                                                       margins = NULL),
                    rkhs.args = NULL, add.joint = TRUE,
                    add.margins = FALSE,
                    joint.mapping = FALSE,
                    neg.weights = FALSE,
                    bf = NULL,
                    sample_weight = NULL,
                    soc = FALSE) {
  
  # estimand <- match.arg(estimand)
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  margmass = get_sample_weight(sample_weight, z = z)
  joint.mapping <- isTRUE(joint.mapping)
  penalty <- match.arg(penalty)
  neg.weights <- isTRUE(neg.weights)
  soc <- isTRUE(soc)
  
  if (neg.weights && penalty == "entropy") {
    stop("Can't have negative weights with the entropy penalty!")
  }
  
  x1 <- x[z == 1,,drop = FALSE]
  x0 <- x[z == 0,,drop = FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  weight.dim <- n0 * n1
  
  dist <- match.arg(dist)
  
  if(is.null(cost) ) {
    if(add.margins) {
      cost <- lapply(1:d, function(i) cost_fun(x[,i,drop = FALSE], z, 
                                               ground_p = p, metric = dist,
                                               rkhs.args = rkhs.args, estimand = "ATT"))
      cost[d + 1] <- list(cost_fun(x, z, 
                                   ground_p = p, metric = dist,
                                   rkhs.args = rkhs.args, estimand = "ATT"))
      
    } else {
      cost <- cost_fun(x, z, 
                       ground_p = p, metric = dist,
                       rkhs.args = rkhs.args, estimand = "ATT")
      cost_a <- cost_fun(rbind(x0,x0), c(rep(0,n0), rep(1,n0)), 
                       ground_p = p, metric = dist,
                       rkhs.args = rkhs.args, estimand = "ATT")
      cost_b <- cost_fun(rbind(x1,x1), c(rep(0,n1), rep(1,n1)), 
                         ground_p = p, metric = dist,
                         rkhs.args = rkhs.args, estimand = "ATT")
    }
  }
  
  
  if (add.margins) {
    if (!is.list(cost) && is.matrix(cost)) {
      cost <- c(lapply(1:d, function(i) cost_fun(x[,i,drop = FALSE], z, 
                                                 ground_p = p, metric = dist,
                                                 rkhs.args = rkhs.args, estimand = "ATT")),
                cost)
    }
    d_cost <- length(cost) - 1
    cost.marg <- cost[1:d_cost]
    cost <- cost[[d_cost + 1]]
    
    sapply(cost.marg, function(cc) stopifnot(dim(cc) %in% c(n0,n1)))
    
    marg_cost_vec <- Matrix::sparseMatrix(i = rep(1:d, each = n0 * n1),
                                          j = rep(1:(n0 * n1), d),
                                          x = c(sapply(cost.marg, c)^p),
                                          dims = c(d, n0 * n1))
  } else {
    marg_cost_vec <- cost.marg <- NULL
  }
  
  stopifnot(dim(cost) %in% c(n0,n1))
  
  median_cost_p <- median(cost^p)
  
  joint_cost_vec <- c( cost^p
                      , -0.5 *c(cost_a^p)
                      , -0.5 * c(cost_b^p)
                      ) #/ median_cost_p
  
  K_joint <- set_K(K, x0, x1, d_c = length(cost.marg), joint.mapping, add.margins,
             penalty = penalty, joint.mapping = joint.mapping,
             method = "Wasserstein")
  
  #objective function
  obj <- list(L = c(cost = joint_cost_vec))
  
  # linear constraints
  LC <- list()
  
  #linear constraint matrix A and bounds
  marg_const <- get_marg_const(margmass$b, n0, n1, divergence = TRUE)
  
  cost_const <- rbind(marg_cost_vec
                      # , joint_cost_vec
                      )
  
  LC$A <- rbind(marg_const$A, cost_const)
  # LC$A <- rbind(marg_const$A, 
  #               cbind(cost_const,
  #               zero_mat_sp(nrow(cost_const),
  #                           ncol(marg_const$A) - length(joint_cost_vec))))
  
  # linear constraint values
  K_const  <- c(margins = K$margins)
  
  # set constraint bounds
  LC$uc <- c(marg_const$lc, marginal_cost_const = K_const^p
             # , joint = Inf
             )
  LC$lc <- c(marg_const$uc, marginal_cost_const =rep(0, length(K_const))
             # ,joint = 0
             )
  
  # if (neg.weights) {
  #   K_const <- c(K_const, joint = 0, margins = rep(0, length(K_const)))
  #   LC$A <- rbind(LC$A, -joint_cost_vec, -1 * marg_cost_vec)
  # }
  # 
  # LC$vals <- c(sum_1 = 1, sum_b = marg_const, K_const^p)
  # 
  # # set direction
  # LC$dir <- c(rep("E", 1 + length(marg_const)),
  #             rep("L", length(K_const)))
  
  # create lp/qp
  op <- list(obj = obj, LC = LC)
  
  # add bal constraints if necessary
  op <- add_bc(op, bf, z, K)
  
  # add bounds
  op <- add_bounds(op, neg.weights, joint.mapping)
  
  # add in mapping
  if (joint.mapping) op <- add_mapping(op, x0, x1, p, margmass$b, K$joint/median(cost), penalty) #currently just p = 2 method
  
  # add in penalty functions
  op <- qp_pen(op, n0, n1, margmass$a, margmass$b, penalty = penalty, K$penalty, soc, FALSE)
  op <- qp_pen(op, n0, n0, margmass$a, margmass$a, penalty = penalty, K$penalty, soc, TRUE)
  op <- qp_pen(op, n1, n1, margmass$b, margmass$b, penalty = penalty, K$penalty, soc, TRUE)
  
  op$nvar <- n0*n1
  # K_vals <- rep(K, d + 1)
  # L0 <- rep(0, n0 * n1) #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  # obj <- list(Q = Matrix::Diagonal(weight.dim, 0.5 * K),#Q0,
  #             L = LO)
  # 
  # LC <- list()
  # LC$A <- rbind(sum_const, 
  #               marg_const_mat, 
  #               cost_vec)
  # LC$vals <- c(1, marg_const, K_vals)
  # LC$dir <- c(rep("E", 1 + marg_const_n),
  #             rep("L", length(K_vals)))
  
  # op <- list(obj = obj, LC = LC)
  return(op)
}


qp_wass_const <- function(x, z, K = list(penalty = NULL,
                                         joint   = NULL,
                                         margins = NULL,
                                         balance = NULL), 
                          p = NULL,
                          penalty = c("L2","entropy",
                                      "variance", "none"),
                          dist = dist.metrics(), cost = NULL,
                          rkhs.args = NULL, add.joint = TRUE,
                          add.margins = FALSE,
                          joint.mapping = FALSE,
                          neg.weights = FALSE,
                          bf = NULL,
                          sample_weight = NULL,
                          soc = FALSE) {
  
  # estimand <- match.arg(estimand)
  margmass = get_sample_weight(sample_weight, z = z)
  joint.mapping <- isTRUE(joint.mapping)
  penalty <- match.arg(penalty)
  neg.weights <- isTRUE(neg.weights)
  soc <- isTRUE(soc)
  
  if (neg.weights && penalty == "entropy") {
    stop("Can't have negative weights with the entropy penalty!")
  }
  
  x1 <- x[z == 1,,drop = FALSE]
  x0 <- x[z == 0,,drop = FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  weight.dim <- n0 * n1
  
  dist <- match.arg(dist)
  
  if(is.null(p)) {
    if(penalty != "entropy") {
      p <- floor(d/2 + 1)
    } else {
      p <- 2
    }
  }
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  if(p < d/2 && penalty != "entropy") warning("power of the wasserstein distance is less than rate optimal d/2 for non-entropy penalties.")
  
  if(is.null(cost) ) {
    if(add.margins) {
      cost.marg <- lapply(1:d, function(i) cost_fun(x[,i,drop = FALSE], z, 
                                                    ground_p = p, metric = dist,
                                                    rkhs.args = rkhs.args, estimand = "ATT"))
      if(add.joint) {
        cost.joint <- list(cost_fun(x, z, 
                                    ground_p = p, metric = dist,
                                    rkhs.args = rkhs.args, estimand = "ATT"))
      }
    } else if (add.joint) {
      cost.joint <- cost_fun(x, z, 
                             ground_p = p, metric = dist,
                             rkhs.args = rkhs.args, estimand = "ATT")
    }
  } else {
    if (add.joint & add.margins) {
      cost.marg <- cost[1:d]
      cost.joint <- cost[[d+1]]
    } else if (add.margins) {
      cost.marg <- cost[1:d]
    } else if (add.joint) {
      cost.joint <- cost
    }
  }
  
  
  if (add.margins) {
    sapply(cost.marg, function(cc) stopifnot(dim(cc) %in% c(n0,n1)))
    
    marg_cost_vec <- Matrix::sparseMatrix(i = rep(1:d, each = n0 * n1),
                                          j = rep(1:(n0 * n1), d),
                                          x = c(sapply(cost.marg, c)^p),
                                          dims = c(d, n0 * n1))
  } else {
    marg_cost_vec <- cost.marg <- NULL
  }
  
  if (add.joint) {
    stopifnot(dim(cost.joint) %in% c(n0,n1))
    
    joint_cost_vec <- Matrix::sparseMatrix(i = rep(1, n0 * n1),
                                           j = 1:(n0*n1),
                                           x = c(cost.joint)^p,
                                           dims = c(1, n0*n1))
  } else {
    joint_cost_vec <- NULL
  }
  
  K <- set_K(K, x1, x0, d_c = length(cost.marg),
             add.joint, add.margins, penalty = penalty, 
             joint.mapping = joint.mapping,
             method = "Constrained Wasserstein")
  
  #objective function
  obj <- list(L = c(cost = rep(0, n1 * n0)))
  
  # linear constraints
  LC <- list()
  #linear constraint matrix A
  sum_const_A <- Matrix::sparseMatrix(i = rep(1, n0*n1),
                                      j = 1:(n0*n1),
                                      x = rep(1, n1 * n0),
                                      dims = c(1, n0*n1))
  marg_const_mat_A <- vec_to_col_constraints(n0,n1)
  
  LC$A <- rbind(sum_const_A, 
                marg_const_mat_A, 
                joint_cost_vec,
                marg_cost_vec)
  
  
  # linear constraint values
  sum_const <- 1
  marg_const <- margmass$b 
  K_const  <- c(joint = K$joint, margins = K$margins)
  
  if (neg.weights) {
    K_const <- c(K_const, joint = 0, margins = rep(0, length(K_const) - 1))
    LC$A <- rbind(LC$A, -1 * joint_cost_vec, -1 * marg_cost_vec)
  }
  
  # set constraint bounds
  LC$uc <- c(sum_1 = 1, sum_b = marg_const, K_const^p)
  LC$lc <- c(sum_1 = 1, sum_b = marg_const, rep(0, length(K_const)))
  
  # create lp/qp
  op <- list(obj = obj, LC = LC)
  
  # add bal constraints if necessary
  op <- add_bc(op, bf, z, K)
  
  # add bounds
  op <- add_bounds(op, neg.weights, joint.mapping)
  
  # add in mapping
  if (joint.mapping) op <- add_mapping(op, x0, x1, p, margmass$b, K$joint, penalty) #currently just p = 2 method
  
  # add in penalty functions
  op <- qp_pen(op, n0, n1, margmass$a, margmass$b, penalty = penalty, K$penalty, soc, FALSE)
  
  # K_vals <- rep(K, d + 1)
  # L0 <- rep(0, n0 * n1) #simple_triplet_matrix(i=integer(0), j = integer(0), v = numeric(0),
  #nrow = 1, ncol = n0*n1) #c(rep(0, n0*n1)) #c( w = rep( -marg_mass, n0 * n1 ) )
  # obj <- list(Q = Matrix::Diagonal(weight.dim, 0.5 * K),#Q0,
  #             L = LO)
  # 
  # LC <- list()
  # LC$A <- rbind(sum_const, 
  #               marg_const_mat, 
  #               cost_vec)
  # LC$vals <- c(1, marg_const, K_vals)
  # LC$dir <- c(rep("E", 1 + marg_const_n),
  #             rep("L", length(K_vals)))
  
  # op <- list(obj = obj, LC = LC)
  return(op)
}

qp_scm <- function(x, z, K = list(penalty = NULL,
                                  joint   = NULL,
                                  margins = NULL,
                                  balance = NULL), 
                   p = 2,
                   penalty = c("L2","entropy",
                               "variance","none"),
                   dist = dist.metrics(), cost = NULL,
                   rkhs.args = NULL, add.joint = TRUE,
                   add.margins = FALSE,
                   joint.mapping = FALSE,
                   neg.weights = FALSE,
                   bf = NULL,
                   sample_weight = NULL,
                   soc = FALSE) {
  
  # estimand <- match.arg(estimand)
  stopifnot(is.numeric(p))
  stopifnot(length(p) == 1)
  margmass = get_sample_weight(sample_weight, z = z)
  joint.mapping <- isTRUE(joint.mapping)
  penalty <- match.arg(penalty)
  neg.weights <- isTRUE(neg.weights)
  soc <- isTRUE(soc)
  
  if (neg.weights && penalty == "entropy") {
    stop("Can't have negative weights with the entropy penalty!")
  }
  
  x1 <- x[z == 1,,drop = FALSE]
  x0 <- x[z == 0,,drop = FALSE]
  
  n <- nrow(x)
  d <- ncol(x)
  
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  weight.dim <- n0 * n1
  
  dist <- match.arg(dist)
  
  K <- set_K(K, x0, x1, d_c = 0, FALSE, add.margins,
             penalty = penalty, joint.mapping = joint.mapping,
             method = "SCM")
  
  #objective function
  obj <- list(L = rep(0, n0 * n1))
  
  # linear constraints
  LC <- list()
  
  #linear constraint matrix A
  # sum_const_A <- Matrix::sparseMatrix(i = rep(1, n0*n1),
  #                                     j = 1:(n0*n1),
  #                                     x = rep(1, n1 * n0),
  #                                     dims = c(1, n0*n1))
  marg_const_mat_A <- vec_to_col_constraints(n0,n1)
  
  LC$A <- rbind(marg_const_mat_A
                )
  
  
  # linear constraint values
  marg_const <- margmass$b 
  
  LC$lc <- c(sum_b = marg_const)
  LC$uc <- c(sum_b = marg_const)
  
  # LC$vals <- c(sum_b = marg_const)
  # # set direction
  # LC$dir <- c(rep("E", 1 
  #                 + length(marg_const)
  #                 ))
  
  # create lp/qp
  op <- list(obj = obj, LC = LC)
  
  # add bounds
  op <- add_bounds(op, neg.weights, TRUE)
  
  # add in mapping
  op <- add_mapping(op, x0, x1, p, margmass$b, 0.0, penalty) #currently just p = 2 method
  
  # add in penalty functions
  op <- qp_pen(op, n0, n1, margmass$a, margmass$b, penalty = penalty, K$penalty, soc, FALSE)
  #(qp, n0, n1, a, b, penalty, lambda, soc, divergence)
  
  return(op)
}

# qp_ebw <- function(x, z,
#                         p = 2,
#                         dist = dist.metrics(), cost = list(joint = NULL,
#                                                            a = NULL,
#                                                            b = NULL),
#                         sample_weight = NULL) {
#   
#   # estimand <- match.arg(estimand)
#   stopifnot(is.numeric(p))
#   stopifnot(length(p) == 1)
#   margmass = get_sample_weight(sample_weight, z = z)
#   
#   x1 <- x[z == 1,,drop = FALSE]
#   x0 <- x[z == 0,,drop = FALSE]
#   
#   n <- nrow(x)
#   d <- ncol(x)
#   
#   n1 <- nrow(x1)
#   n0 <- nrow(x0)
#   
#   weight.dim <- n0 * n1
#   
#   dist <- match.arg(dist)
#   
#   if(is.null(cost) ) {
#       cost <- cost_fun(x, z, 
#                        ground_p = p, metric = dist,
#                        rkhs.args = rkhs.args, estimand = "ATT")
#       cost_a <- cost_fun(rbind(x0,x0), c(rep(0,n0), rep(1,n0)), 
#                          ground_p = p, metric = dist,
#                          rkhs.args = rkhs.args, estimand = "ATT")
#       cost_b <- cost_fun(rbind(x1,x1), c(rep(0,n1), rep(1,n1)), 
#                          ground_p = p, metric = dist,
#                          rkhs.args = rkhs.args, estimand = "ATT")
#   }
#   
#   
#   stopifnot(dim(cost) %in% c(n0,n1))
#   
#   joint_cost_vec <- rowMeans(cost)
#   
#   #objective function
#   obj <- list(L = c(cost = 2 * joint_cost_vec^p/n0),
#               Q = Matrix::sparseMatrix(i = rep(1:n0, n0),
#                                        j = rep(1:n0, each = n0),
#                                        x = -cost_a^p/n0^2,
#                                        dims = c(n0,n0)))
#   
#   # linear constraints
#   LC <- list()
#   LC$A <- Matrix::sparseMatrix(i = rep(1, n0), j = 1:n0,
#                                x = 1, dims = c(1,n0))
#   
#   # set constraint bounds
#   LC$uc <- c(1)
#   LC$lc <- c(1)
#   
#   op <- list(obj = obj, LC = LC)
#   
#   return(op)
# }

qp_rkhs <- function(x, z, p = 1, estimand = c("ATC", "ATT", "ATE"),
                    theta = c(1,1),
                    gamma = c(1,1), lambda = 0, 
                    sigma_2 = c(1,1),
                    kernel = c("RBF", "polynomial"),
                    dist = c("mahalanobis","Lp"), cost = NULL, 
                    is.dose = FALSE, ...) {
  # est <- match.arg(estimand)
  # if (est == "ATC") {
  #   z <- 1 - z
  # }
  n <- nrow(x)
  kernel <- match.arg(kernel)
  dist <- match.arg(dist)
  # cost.fun <- switch(isTRUE(is.dose),
  #                    "TRUE" = "kernel_calculation",
  #                    "FALSE" = 
  
  if (is.null(cost)) { 
    kernels <- kernel_calculation(x, z, p = p, theta = theta, 
                                  gamma = gamma, sigma_2 = sigma_2, 
                                  metric = dist, kernel =  kernel,
                                  is.dose = is.dose,
                                  estimand = estimand)
    cost <- list(kernels$cov_kernel, kernels$mean_kernel)
  } else {
    stopifnot(all(dim(cost) %in% n))
  }
  
  cost[[1]] <- Matrix::Matrix(data = cost[[1]], sparse = TRUE)
  
  if (is.dose) {
    Q0 <- (1/(n^2)) * (cost[[1]] + Matrix::Diagonal(n, lambda/(n^2)))
    
    L0 <- c(w = -2/(n^2) * c(cost[[2]]))
    
    A <- t(rep(1.0/n, n))
    vals <- as.double(n)
    # dir <- "E"
    uc <- lc <- vals
    
  } else {
    Q0 <- cost[[1]]
    
    L0 <- c(w = -2 * c(cost[[2]]))
    
    A <- rbind(t(1/n * z),
               t(1/n * (1 - z)  ))
    vals <- as.double(c(n,n))
    # dir <- c("E", "E")
    uc <- lc <- vals
    
  }
  
  Q0 <- as(Q0, "dsTMatrix")
  
  
  quick_op <- list(obj = list(Q = Q0, L = L0),
                   LC = list(A = A, 
                             lc = lc,
                             uc = uc
                             # dir = dir,
                             # vals = vals
                             ))
  quick_op$nvar <- length(L0)
  quick_op$bounds <- list(lb = rep(0,   quick_op$nvar),
                          ub = rep(Inf, quick_op$nvar))
  return(quick_op)
}

qp_lin_comb <- function(x) {
  n <- length(x)
  
  #objective function
  obj <- list(L = c(x))
  
  # linear constraints
  LC <- list()
  LC$A <- Matrix::sparseMatrix(i = rep(1, n), j = 1:n,
                               x = 1, dims = c(1,n))
  
  # set constraint bounds
  LC$uc <- c(1)
  LC$lc <- c(1)
  
  op <- list(obj = obj, LC = LC)
  op$bounds <- list(lb = rep(0,n), ub = rep(Inf,n))
  op$nvar <- n
  return(op)
}

qp_dual <- function(f, g, b, dispersion = numeric(0)) {
  n <- length(f)
  
  #objective function
  obj <- list(L = c(f))
  
  # linear constraints
  LC <- list()
  LC$A <- rbind(Matrix::sparseMatrix(i = rep(1, n), j = 1:n,
                               x = 1, dims = c(1,n)),
                c(f))
  
  # set constraint bounds
  LC$uc <- c(1, Inf)
  LC$lc <- c(1, -sum(g * b))
  
  op <- list(obj = obj, LC = LC)
  op$bounds <- list(lb = rep(0,n), ub = rep(Inf,n))
  
  # dispersion
  if (length(dispersion) > 0) {
    op$cones <- list(F = rbind(c(rep(0, n ),1),
                               rep(0, n + 1),
                         cbind(Matrix::Diagonal(n, 1), zero_mat_sp(n,1))),
                     cones =  matrix(list("RQUAD",  n + 2, NULL),
                                     nrow = 3, ncol = 1),
                     g = c(0,1, rep(0, n)))
    op$obj$L <- c(op$obj$L, dispersion)
    op$LC$A <- cbind(op$LC$A, zero_mat_sp(nrow(op$LC$A),1))
    op$bounds$lb <- c(op$bounds$lb, 0)
    op$bounds$ub <- c(op$bounds$ub, Inf)
  }
  
  op$nvar <- n
  return(op)
}

qp_dual_max <- function(f, g, b, dispersion = numeric(0)) {
  n <- length(f)
  
  #objective function
  obj <- list(L = c(-f))
  
  # linear constraints
  LC <- list()
  LC$A <- rbind(Matrix::sparseMatrix(i = rep(1, n), j = 1:n,
                                     x = 1, dims = c(1,n)),
                c(-f))
  
  # set constraint bounds
  LC$uc <- c(1, sum(g * b))
  LC$lc <- c(1, -Inf)
  
  op <- list(obj = obj, LC = LC)
  op$bounds <- list(lb = rep(0,n), ub = rep(Inf,n))
  
  # dispersion
  if (length(dispersion) > 0) {
    op$cones <- list(F = rbind(c(rep(0, n ),1),
                               rep(0, n + 1),
                               cbind(Matrix::Diagonal(n, 1), zero_mat_sp(n,1))),
                     cones =  matrix(list("RQUAD",  n + 2, NULL),
                                     nrow = 3, ncol = 1),
                     g = c(0,1, rep(0, n)))
    op$obj$L <- c(op$obj$L, dispersion)
    op$LC$A <- cbind(op$LC$A, zero_mat_sp(nrow(op$LC$A),1))
    op$bounds$lb <- c(op$bounds$lb, 0)
    op$bounds$ub <- c(op$bounds$ub, Inf)
  }
  
  op$nvar <- n
  return(op)
}

qp_proj <- function(f, g, a, b, BC) {
  n <- length(a)
  
  #objective function
  obj <- list(L = c(-a),
              Q = Matrix::Diagonal(n,x = 1))
  
  # linear constraints
  LC <- list()
  LC$A <- rbind(1, Matrix::Matrix(data = t(BC$source), sparse = TRUE))
  
  # set constraint bounds
  sds <- matrixStats::colSds(BC$target)
  delta_star <- BC$K * sds
  E_target <- colMeans(BC$target)
  
  LC$uc <- c(1, delta_star + E_target)
  LC$lc <- c(1, -delta_star + E_target)
  
  op <- list(obj = obj, LC = LC)
  op$bounds <- list(lb = rep(0,n), ub = rep(Inf,n))
  
  
  op$nvar <- n
  return(op)
}

qp_proj_update <- function(f, g, a, b, op) {
  op$obj$L <- c(-a)
  
  return(op)
}

lp_min_constraint <- function(f, g, a, b, BC) {
  n <- length(a)
  
  #objective function
  obj <- list(L = c(f))
  
  # linear constraints
  LC <- list()
  LC$A <- rbind(1, Matrix::Matrix(data = t(BC$source), sparse = TRUE))
  
  # set constraint bounds
  sds <- matrixStats::colSds(BC$target)
  delta_star <- BC$K * sds
  E_target <- colMeans(BC$target)
  
  LC$uc <- c(1, delta_star + E_target)
  LC$lc <- c(1, -delta_star + E_target)
  
  op <- list(obj = obj, LC = LC)
  op$bounds <- list(lb = rep(0,n), ub = rep(Inf,n))
  
  
  op$nvar <- n
  return(op)
}

lp_min_constraint_update <- function(f, g, a, b, op) {
  op$obj$L <- c(f)
  
  return(op)
}


check_wass_const <- function(opt_problem) {
  cost <- c(opt_problem$LC$A[1,])
  const <- c(opt_problem$LC$vals[1])
  n <- length(cost)
  mass <- rep.int(1/n,n)
  wass <- sum(cost * mass)
  
  output <- if (wass < const) {
    list(res = mass, skip_cplex = TRUE)
  } else {
    list(res = NULL, skip_cplex = FALSE)
  }
  return(output)
}

# update.qp <- function(object, new, method, add.margins) {
#   if (method == "Wasserstein") {
#     object$obj$Q <- Matrix::Diagonal(nrow(object$obj$Q), new[length(new)] * 0.5)
#   } if (method == "Constrained Wasserstein") {
#     object$LC$
#   }
# }

# setOldClass("DataSim")
setGeneric("quadprog", function(data, ...) UseMethod("quadprog"))
setMethod("quadprog", "DataSim", quadprog.DataSim)
setMethod("quadprog", "data.frame", quadprog.data.frame)
setMethod("quadprog", "matrix", quadprog.data.frame)

