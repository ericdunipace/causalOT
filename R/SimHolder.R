{
  SimHolder <- R6::R6Class("SimHolder",
                           public = list(run = function() {
                             private$output <- vector("list", private$nsim)
                             
                             for(iter in 1:private$nsim) {
                               if(private$verbose) message("Iteration: ",iter," of ", private$nsim)
                               private$iter <- iter
                               private$output[[iter]] <- private$update()
                             }
                           },
                           get.output = function() {
                             full <- data.table::rbindlist(private$output)
                             exclude <- "options"
                             exclude.cn <- -which(colnames(private$output[[1]]) == exclude)
                             
                             output <- full[,exclude.cn]
                             excluded <- data.table::rbindlist(lapply(private$output, function(o) 
                               data.table::rbindlist(lapply(o[[exclude]], function(e) as.data.frame(lapply(unlist(e), function(i) i))), 
                                                     use.names = TRUE, 
                                                     fill = TRUE)),   
                               use.names = TRUE, 
                               fill = TRUE)
                             data.table::setDT(output)
                             for(j in colnames(excluded)) data.table::set(output, i = NULL, j = j,  value = excluded[,j])
                             
                             return(output)
                           },
                           get.outcome = function(output) {
                             exclude <- c("wasserstein", "ess.frac", "psis.k", "psis.ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             return(output[,exclude.cn])
                           },
                           get.ESS.frac = function(output) {
                             exclude <- c("wasserstein", "estimate", "ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             ESS.output <- output[,exclude.cn]
                             data.table::setDT(ESS.output)
                             
                             ESS <- data.table::rbindlist(lapply(output[,"ess.frac"], function(o) as.data.frame(t(o))))
                             colnames(ESS) <- c("ESS.frac.control", "ESS.frac.treated")
                             
                             for(j in colnames(ESS)) data.table::set(ESS.output, i = NULL, j = j, value = ESS[,j])
                             
                             return(ESS.output)
                           },
                           get.diagnostics = function(output) {
                             exclude <- c("wasserstein", "estimate", "ess.frac", "psis.k", "psis.ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             ESS.output <- output[,exclude.cn]
                             data.table::setDT(ESS.output)
                             
                             ESS <- data.table::rbindlist(lapply(output[,"ess.frac"], function(o) as.data.frame(t(o))))
                             colnames(ESS) <- c("ESS.frac.control", "ESS.frac.treated")
                             
                             psis.ess <- data.table::rbindlist(lapply(output[,"psis.ess.frac"], function(o) as.data.frame(t(o))))
                             colnames(psis.ess) <- c("psis.ESS.frac.control", "psis.ESS.frac.treated")
                             
                             psis.k <- data.table::rbindlist(lapply(output[,"psis.k"], function(o) as.data.frame(t(o))))
                             colnames(psis.k) <- c("psis.k.control", "psis.k.treated")
                             
                             for(j in colnames(ESS)) data.table::set(ESS.output, i = NULL, j = j, value = ESS[,j])
                             for(j in colnames(psis.ess)) data.table::set(ESS.output, i = NULL, j = j, value = psis.ess[,j])
                             for(j in colnames(psis.k)) data.table::set(ESS.output, i = NULL, j = j, value = psis.k[,j])
                             
                             return(ESS.output)
                           },
                           get.psis = function(output) {
                             exclude <- c("wasserstein", "estimate", "ess.frac", "psis.k", "psis.ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             psis.output <- output[,exclude.cn]
                             data.table::setDT(psis.output)
                             
                             
                             psis.ess <- data.table::rbindlist(lapply(output[,"psis.ess.frac"], function(o) as.data.frame(t(o))))
                             colnames(psis.ess) <- c("psis.ESS.frac.control", "psis.ESS.frac.treated")
                             
                             psis.k <- data.table::rbindlist(lapply(output[,"psis.k"], function(o) as.data.frame(t(o))))
                             colnames(psis.k) <- c("psis.k.control", "psis.k.treated")
                             
                             for(j in colnames(psis.ess)) data.table::set(psis.output, i = NULL, j = j, value = psis.ess[,j])
                             for(j in colnames(psis.k)) data.table::set(psis.output, i = NULL, j = j, value = psis.k[,j])
                             
                             return(psis.output)
                           },
                           get.wass = function(output) {
                             wass.dt <- data.table::rbindlist(lapply(output[,"wasserstein"], function(o) as.data.frame(o[1:3])), idcol = TRUE, fill = TRUE)
                             dist <- unlist(lapply(output[,"wasserstein"], function(o) as.data.frame(o[4])), recursive = FALSE)
                             wass.dt <- data.table::set(wass.dt, i = NULL, j = "dist", value = dist)
                             check.rep <- table(wass.dt[,".id"])
                             stopifnot(all(check.rep == check.rep[1]))
                             nrep <- check.rep[1]
                             
                             exclude <- c("wasserstein", "estimate", "ESS.frac","psis.k","psis.ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             wass.output <- output[rep(1:nrow(output), each = nrep), exclude.cn]
                             data.table::setDT(wass.output)
                             cn <- colnames(wass.dt)
                             # cn <- cn[cn!= ".id"]
                             
                             for(j in cn) data.table::set(wass.output, i = NULL, j = j, value = wass.dt[,j])
                             
                             return(wass.output)
                           },
                           initialize = function(nsim = 100,
                                                 dataSim,
                                                 grid.search = TRUE,
                                                 truncations = NULL,
                                                 standardized.difference.means = NULL,
                                                 estimands = c("ATT","ATC","cATE","ATE"),
                                                 RKHS = list(lambdas = NULL, theta = NULL, gamma = NULL, p = NULL,
                                                             metric = "malahanobis",
                                                             kernel = "RBF",
                                                             sigma_2 = NULL, opt = NULL, opt.method = NULL),
                                                 outcome.model = "lm",
                                                 outcome.formula = list(none = NULL,
                                                                        augmentation = NULL),
                                                 model.augmentation = "both",
                                                 match = "both",
                                                 calculate.feasible = FALSE,
                                                 solver = "gurobi",
                                                 wass_powers = 2,
                                                 ground_powers = 2,
                                                 metrics = "Lp", 
                                                 constrained.wasserstein.target = c("RKHS", "SBW"),
                                                 cluster = FALSE,
                                                 verbose = FALSE) {
                             private$nsim <- nsim
                             private$cluster <- cluster
                             private$verbose <- isTRUE(verbose)
                             
                             if(is.null(dataSim) | !inherits(dataSim, "DataSim")) {
                               stop("DataSim class must be given")
                             } else {
                               private$simulator <- dataSim
                             }
                             if(is.null(truncations)) {
                               private$truncations <- 0
                             } else {
                               private$truncations <- truncations
                             }
                             if(!is.null(standardized.difference.means)) {
                               private$standardized.difference.means <- standardized.difference.means
                             } else {
                               private$standardized.difference.means <- NULL
                             }
                             private$metric <- match.arg(metrics,
                                                         c("Lp", "mahalanobis", "RKHS"), several.ok = TRUE)
                             if(!is.null(RKHS) & !missing(RKHS)) {
                               if(!is.list(RKHS)) RKHS <- list(RKHS)
                               private$RKHS <- list()
                               
                               if(is.null(RKHS$lambdas)) RKHS$lambdas <- seq(0,100, length.out = 11)
                               private$RKHS$lambdas <-  RKHS$lambdas
                               
                               if(is.null(RKHS$theta)) RKHS$theta <- as.double(c(1,1))
                               private$RKHS$theta <- RKHS$theta
                               
                               if(is.null(RKHS$gamma)) RKHS$gamma <- as.double(c(1,1))
                               private$RKHS$gamma <- RKHS$gamma
                               
                               if(is.null(RKHS$opt)) RKHS$opt <- TRUE
                               private$RKHS$opt <- isTRUE(RKHS$opt)
                               
                               if(is.null(RKHS$p)) {
                                 if(private$RKHS$opt) {
                                   RKHS$p <- as.double(2:4)
                                 } else {
                                   RKHS$p <- as.double(2)
                                 }
                               }
                               private$RKHS$p <- RKHS$p
                               
                               
                               if(is.null(RKHS$sigma_2)) RKHS$sigma_2 <- as.double(1)
                               private$RKHS$sigma_2  <- RKHS$sigma_2 
                               
                               
                               
                               if(is.null(RKHS$opt.method)) RKHS$opt.method <- "stan"
                               private$RKHS$opt.method <- RKHS$opt.method
                               if(!is.null(RKHS$iter)) {
                                 private$RKHS$iter <- RKHS$iter
                               }
                               
                               if(is.null(RKHS$metric)) {
                                 RKHS$metric <- "mahalanobis"
                                 # if(all(private$metric == "RKHS")) {
                                 #   RKHS$metric <- "mahalanobis"
                                 # } else {
                                 #   RKHS$metric <- private$metric[private$metric != "RKHS"]
                                 # }
                               }
                               private$RKHS$metric <- RKHS$metric
                               if(is.null(RKHS$kernel)) RKHS$kernel <- "RBF"
                               private$RKHS$kernel <- RKHS$kernel
                             } else {
                               private$RKHS <- list()
                               private$RKHS$lambdas <-  as.double(seq(0,100, length.out = 11))
                               private$RKHS$lambdas[1] <- 1e-3
                               private$RKHS$theta <- as.double(c(1,1))
                               private$RKHS$gamma <- as.double(c(1,1))
                               private$RKHS$p <- as.double(2:4)
                               private$RKHS$metric <- "mahalanobis"
                               # if(all(private$metric == "RKHS")) {
                               #     private$RKHS$metric <- "mahalanobis"
                               # } else {
                               #   private$RKHS$metric <- private$metric[private$metric != "RKHS"]
                               # }
                               private$RKHS$opt <- TRUE
                               private$RKHS$sigma_2 <- as.double(1)
                               private$RKHS$opt.method <- "stan"
                               private$RKHS$kernel <- "RBF"
                             }
                             if(is.null(grid.search) & (is.null(standardized.difference.means) | is.null(RKHS))) {
                               private$grid.search <- TRUE
                             } else {
                               private$grid.search <- isTRUE(grid.search)
                             }
                             if(!is.null(outcome.model)) {
                               if(!is.character(outcome.model)) outcome.model <- as.character(outcome.model)
                               private$outcome.model <- outcome.model
                             } else {
                               private$outcome.model <- "lm"
                             }
                             if(!is.null(outcome.formula)) {
                               if(!is.list(outcome.formula)) outcome.formula <- list(outcome.formula)
                               private$outcome.formula <- outcome.formula
                             } else {
                               private$outcome.formula <- list(none = NULL, augmentation = NULL)
                             }
                             
                             private$solver <- match.arg(solver, c("gurobi","cplex","mosek"))
                             private$wass_powers <- as.numeric(wass_powers)
                             private$ground_powers <- as.numeric(ground_powers)
                             
                             private$calculate.feasible <- isTRUE(calculate.feasible)
                             options <- get_holder_options()
                             if(!is.null(estimands)) {
                               estimands <- match.arg(estimands, several.ok = TRUE)
                             } else {
                               estimands <- options$estimates
                             }
                             private$estimand <- options$estimates[options$estimates %in% estimands]
                             if(!private$calculate.feasible) private$estimand <- private$estimand[private$estimand != "feasible"]
                             #removing "wasserstein"
                             private$method <- c(options$weights[options$weights != "Wasserstein" &
                                                                   options$weights != "RKHS.dose"],
                                                 "gp")
                             private$weight_based_methods <- private$method[private$method != "gp"]
                             private$nmethod <- length(private$method)
                             private$model.augmentation <- match.arg(model.augmentation, c("both", "yes", "no"))
                             private$match <- match.arg(match, c("both", "yes", "no"))
                             # private$cost.setup()
                             
                             private$temp.output <- vector("list", length = private$nmethod)
                             names(private$temp.output) <- private$method
                             
                             n_g <- length(private$ground_powers)
                             n_m <- length(private$metric)
                             n_w <- length(private$wass_powers)
                             nrows <- n_g * n_m * n_w
                             
                             private$wass_df <- list(metric = rep(private$metric, each = n_g * n_w))
                             private$wass_df$ground_p <- rep_len(rep(private$ground_powers, each = n_w), length.out = nrows)
                             private$wass_df$wass_p <- rep_len(private$wass_powers, length.out = nrows)
                             private$wass_df$dist <- vector("list",nrows)
                             private$cwass.targ <- match.arg(constrained.wasserstein.target)
                             private$method.setup()
                             
                             private$output.dt <- data.table::data.table(
                               method = rep(private$method[1], private$max.conditions),
                               estimand = character(private$max.conditions),
                               model = character(private$max.conditions),
                               model.augmentation = logical(private$max.conditions),
                               match = logical(private$max.conditions),
                               solver = character(private$max.conditions),
                               delta = numeric(private$max.conditions),
                               options = lapply(1:private$max.conditions, function(i) list(NA)), #vector("list", private$max.conditions),
                               wasserstein = lapply(1:private$max.conditions, function(i) private$wass_df),
                               ess.frac = lapply(1:private$max.conditions, function(i) private$wass_df),
                               psis.ess.frac = lapply(1:private$max.conditions, function(i) private$wass_df),
                               psis.k = lapply(1:private$max.conditions, function(i) private$wass_df),
                               estimate = rep(NA_real_, private$max.conditions)
                             )
                             private$weights <- sapply(private$estimand, function(i) {NULL}, simplify = NULL)
                           }
                           ),
                           private = list(calculate.feasible = "logical",
                                          cluster = "logical",
                                          costs = "list",
                                          cwass.targ = "character",
                                          estimand = "vector",
                                          grid.search = "logical",
                                          ground_powers = "vector",
                                          iter = "integer",
                                          output = "list",
                                          output.dt = "data.table",
                                          outcome.formula = "list",
                                          outcome.model = "list",
                                          match = "character",
                                          max.conditions = "integer",
                                          method = "vector",
                                          method.lookup = "data.frame",
                                          metric = "vector",
                                          model.augmentation = "character",
                                          nmethod = "integer",
                                          nsim = "integer",
                                          RKHS = "list",
                                          RKHS.opt = "list",
                                          simulator = "DataSim",
                                          solver = "character",
                                          standardized.difference.means = "vector",
                                          temp.output = "list",
                                          truncations = "vector",
                                          verbose = "logical",
                                          wass_df = "list",
                                          wass_powers = "vector",
                                          weights = "list",
                                          weight_based_methods = "character",
                                          #   cost.fun = function(X, Y, ground_p, direction, metric) {
                                          #     return(switch(metric,
                                          #             "Lp" = causalOT::cost_calc_lp(X, Y, ground_p, direction),
                                          # "mahalanobis" = causalOT::cost_mahalanobis(X, Y, ground_p, direction)))
                                          #     },
                                          cost.setup = function() {
                                            n_g <- length(private$ground_powers)
                                            n_m <- length(private$metric)
                                            nrows <- n_g * n_m
                                            ns <- private$simulator$get_n()
                                            # private$costs <- vector("list", n_g)
                                            # names(private$costs) <- private$metric
                                            # list(metric = rep(private$metric, each = n_g))
                                            private$costs <- sapply(private$metric, 
                                                                    function(i) sapply(as.character(private$ground_powers), 
                                                                                       function(j) matrix(0, nrow = ns["n0"], ncol = ns["n1"] ), simplify = FALSE), simplify = FALSE)
                                            # private$costs$cost <- lapply(1:nrows, function(i) matrix(0, nrow = ns["n0"], ncol = ns["n1"] ))   
                                            x0 <- private$simulator$get_x0()
                                            x1 <- private$simulator$get_x1()
                                            for(i in as.character(private$metric)) {
                                              if(i == "RKHS") next
                                              for (j in as.character(private$ground_powers)) {
                                                # function(X, Y, ground_p, direction, metric)
                                                private$costs[[i]][[j]] <- cost_fun(x=x0, y=x1,
                                                                                    power = as.numeric(j), 
                                                                                    metric = i)
                                              }
                                            }
                                            private$RKHS.opt <- vector("list", length(private$estimand))
                                            names(private$RKHS.opt) <- as.character(private$estimand)
                                            for (i in as.character(private$estimand)) {
                                              # private$RKHS.opt[[i]] <- list()
                                              if(i == "feasible" | i == "cATE") {
                                                next
                                              } else {
                                                est <- i
                                              }
                                              # for(j in as.character(private$RKHS$metric)) {
                                              private$RKHS.opt[[i]] <- 
                                                if (private$RKHS$opt) {
                                                  RKHS_param_opt(x = private$simulator$get_x(),
                                                                 z = private$simulator$get_z(),
                                                                 y = private$simulator$get_y(),
                                                                 power = private$RKHS$p,
                                                                 metric = private$RKHS$metric,
                                                                 kernel = private$RKHS$kernel[1],
                                                                 is.dose = FALSE,
                                                                 opt.method = private$RKHS$opt.method,
                                                                 estimand = est)
                                                } else {
                                                  list(theta = private$RKHS$theta,
                                                       gamma = private$RKHS$gamma,
                                                       p = private$RKHS$p,
                                                       sigma_2 = private$RKHS$sigma_2)
                                                }
                                              private$RKHS.opt[[i]]$is.dose <- FALSE
                                              private$RKHS.opt[[i]]$estimand <- est
                                            }
                                            if("RKHS" %in% as.character(private$metric)) {
                                              private$costs[["RKHS"]] <- vector("list", length(private$estimand))
                                              names(private$costs[["RKHS"]] ) <- as.character(private$estimand)
                                              for (i in as.character(private$estimand)) {
                                                if(i == "feasible" | i == "cATE") {
                                                  next
                                                } else {
                                                  est <- i
                                                }
                                                private$costs[["RKHS"]][[i]] <- cost_fun(x=x0, y=x1,
                                                                                       power = 2, 
                                                                                       metric = "RKHS",
                                                                                       kernel = private$RKHS$kernel[1],
                                                                                       rkhs.args = private$RKHS.opt[[i]],
                                                                                       estimand = est
                                                    )
                                              }
                                              # }
                                              if("feasible" %in% private$estimand) {
                                                private$RKHS.opt[["feasible"]] <- private$RKHS.opt[["ATE"]]
                                                private$costs[["RKHS"]][["feasible"]] <- private$costs[["RKHS"]][["ATE"]]
                                              }
                                              private$RKHS.opt[["cATE"]] <- private$RKHS.opt[["ATE"]] # not correct parameters but functions don't need these parameters
                                              private$costs[["RKHS"]][["cATE"]] <- private$costs[["RKHS"]][["ATE"]]
                                            }
                                            
                                          },
                                          check.skip = function(weights) {
                                            skip <- FALSE
                                            if(isTRUE(all(is.na(weights$w0))) | isTRUE(length(weights$w0) == 0) | isTRUE(is.null(weights$w0))) {
                                              skip <- TRUE
                                            }
                                            if(isTRUE(all(is.na(weights$w1))) | isTRUE(length(weights$w1) == 0) | isTRUE(is.null(weights$w1))) {
                                              skip <- TRUE
                                            }
                                            return(skip)
                                          },
                                          estimate = function(method) {
                                            cur <- private$method.lookup[private$method.lookup$method == method,]
                                            ess.frac <- list()
                                            psis.output <- list()
                                            ns <- private$simulator$get_n()
                                            # if(method == "Constrained Wasserstein") cur$options <- private$get_delta(cur$options[[1]])
                                            iter <- 1L
                                            data.table::set(private$output.dt, i = NULL, j = "method" , value = as.character(cur$method))
                                            wass.df <- private$wass_df
                                            for (solver in cur$solver) {
                                              for (o in cur$options[[1]]) {
                                                for (est in cur$estimand[[1]]) {
                                                  delta <- private$get_delta(o, est, method)
                                                  if ( isTRUE(is.null(delta)) ) next
                                                  if ( isTRUE(method == "RKHS" | method == "RKHS.dose" | method == "Wasserstein") & isTRUE(est == "feasible") ) next
                                                  if ( isTRUE(method == "RKHS.dose") & isFALSE(est == "ATE") ) next
                                                  if ( isTRUE(method == "Constrained Wasserstein") & isTRUE(est == "feasible")) {
                                                    if(o$metric == "RKHS") next
                                                  }
                                                  private$weight.calc(cur = cur, 
                                                                      estimand = est, 
                                                                      solver = solver,
                                                                      delta = delta,
                                                                      cost = private$get_cost(o, method, est), 
                                                                      p = private$get_power(o),
                                                                      grid.search = isTRUE(o$grid.search),
                                                                      opt.hyperparam = isTRUE(o$opt),
                                                                      opt.method = o$opt.method,
                                                                      metric = o$metric
                                                  )
                                                  if(private$check.skip(private$weights[[est]])) next
                                                  
                                                  ess.frac <- list(ESS(private$weights[[est]])/ns)
                                                  psis.output <- PSIS_diag(private$weights[[est]])
                                                  psis.ess.frac <- list(sapply(psis.output, function(w) w$n_eff)/ns)
                                                  psis.k <- list(lapply(psis.output, function(w) w$pareto_k))
                                                  # if(method == "Constrained Wasserstein" & est == "ATE") if(o$metric=="RKHS") browser()
                                                  for (mods in cur$outcome.model[[1]]) {
                                                    for (aug in cur$model.aug[[1]]) {
                                                      opt.dist <- list(private$simulator$opt_weight_dist(weight = private$weights[[est]], 
                                                                                                         estimand = est, augment = aug, 
                                                                                                         solver = private$solver[1]))
                                                      for (match in cur$match[[1]]) {
                                                        data.table::set(private$output.dt, i = iter, j = "estimand" , value = est)
                                                        data.table::set(private$output.dt, i = iter, j = "model" , value = mods)
                                                        data.table::set(private$output.dt, i = iter, j = "model.augmentation" , value = aug)
                                                        data.table::set(private$output.dt, i = iter, j = "match" , value = match)
                                                        data.table::set(private$output.dt, i = iter, j = "solver" , value = solver)
                                                        data.table::set(private$output.dt, i = iter, j = "delta" , value = if(!is.null(o$delta)){o$delta} else {NA_real_})
                                                        data.table::set(private$output.dt, i = iter, j = "options" , value = list(list(o)))
                                                        private$wass.calc(iter, est)
                                                        data.table::set(private$output.dt, i = iter, j = "ess.frac", value = ess.frac)
                                                        #PSIS
                                                        data.table::set(private$output.dt, i = iter, j = "psis.ess.frac", value = psis.ess.frac)
                                                        data.table::set(private$output.dt, i = iter, j = "psis.k", value = psis.k)
                                                        #opt wt dist
                                                        data.table::set(private$output.dt, i = iter, j = "opt.dist", value = opt.dist)
                                                        data.table::set(private$output.dt, i = iter, j = "estimate", 
                                                                        value = estimate_effect(private$simulator, 
                                                                                                formula = cur$outcome.formula[[1]][[aug + 1]],
                                                                                                weights = private$weights[[est]],
                                                                                                hajek = TRUE,
                                                                                                doubly.robust = aug,
                                                                                                matched = match,
                                                                                                target = est,
                                                                                                model = match.fun(mods))
                                                        )
                                                        iter <- iter + 1L
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          },
                                          model_estimate = function(method) {
                                            cur <- private$method.lookup[private$method.lookup$method == method,]
                                            ess.frac <- list()
                                            psis.output <- list()
                                            ns <- private$simulator$get_n()
                                            # if(method == "Constrained Wasserstein") cur$options <- private$get_delta(cur$options[[1]])
                                            iter <- 1L
                                            data.table::set(private$output.dt, i = NULL, j = "method" , value = as.character(cur$method))
                                            # wass.df <- private$wass_df
                                            for (solver in cur$solver) {
                                              for (o in cur$options[[1]]) {
                                                for (est in cur$estimand[[1]]) {
                                                  if(est == "feasible") next
                                                  
                                                  private$weights[[est]] <- list(w0 = rep(1/ns[1], ns[1]),
                                                                                 w1 = rep(1/ns[2], ns[2]),
                                                                                 gamma = NULL,
                                                                                 estimand = est)
                                                  class(private$weights[[est]]) <- "causalWeights"
                                                  
                                                  # ess.frac <- list(ESS(private$weights[[est]])/ns)
                                                  # psis.output <- PSIS_diag(private$weights[[est]])
                                                  # psis.ess.frac <- list(sapply(psis.output, function(w) w$n_eff)/ns)
                                                  # psis.k <- list(lapply(psis.output, function(w) w$pareto_k))
                                                  # if(method == "Constrained Wasserstein" & est == "ATE") if(o$metric=="RKHS") browser()
                                                  # for (mods in cur$outcome.model[[1]]) {
                                                  # for (aug in cur$model.aug[[1]]) {
                                                  # opt.dist <- list(private$simulator$opt_weight_dist(weight = private$weights[[est]], 
                                                  #                                                    estimand = est, augment = aug, 
                                                  #                                                    solver = private$solver[1]))
                                                  # for (match in cur$match[[1]]) {
                                                  data.table::set(private$output.dt, i = iter, j = "estimand" , value = est)
                                                  data.table::set(private$output.dt, i = iter, j = "model" , value = "gp")
                                                  data.table::set(private$output.dt, i = iter, j = "model.augmentation" , value = NA)
                                                  data.table::set(private$output.dt, i = iter, j = "match" , value = NA)
                                                  data.table::set(private$output.dt, i = iter, j = "solver" , value = NA_character_)
                                                  data.table::set(private$output.dt, i = iter, j = "delta" , value = NA_real_)
                                                  data.table::set(private$output.dt, i = iter, j = "options" , value = list(list(o)))
                                                  # private$wass.calc(iter, est)
                                                  data.table::set(private$output.dt, i = iter, j = "ess.frac", value = list(c(Control = 1,
                                                                                                                            Treated = 1)))
                                                  #PSIS
                                                  data.table::set(private$output.dt, i = iter, j = "psis.ess.frac", value = list(c(w0 = 1,
                                                                                                                                      w1 = 1)))
                                                  data.table::set(private$output.dt, i = iter, j = "psis.k", value = list(list(w0 = NA_real_,
                                                                                                                               w1 = NA_real_)))
                                                  #opt wt dist
                                                  data.table::set(private$output.dt, i = iter, j = "opt.dist", value = NA_real_)
                                                  tau <- if(est == "cATE") {
                                                    match.fun(cur$outcome.model[[1]])(data = private$simulator, 
                                                                                      formula = cur$outcome.formula[[1]][[1]],
                                                                                      weights = private$weights[[est]],
                                                                                      param = private$RKHS.opt[["ATE"]],
                                                                                      estimand = est)
                                                    
                                                  } else {
                                                    match.fun(cur$outcome.model[[1]])(data = private$simulator, 
                                                                                      formula = cur$outcome.formula[[1]][[1]],
                                                                                      weights = private$weights[[est]],
                                                                                      param = private$RKHS.opt[[est]],
                                                                                      estimand = est)
                                                  }
                                                  data.table::set(private$output.dt, i = iter, j = "estimate", 
                                                                  value = c(tau))
                                                  iter <- iter + 1L
                                                  # }
                                                  # }
                                                  # }
                                                }
                                              }
                                            }
                                          },
                                          get_delta = function(opts, estimand, method) {
                                            if(method == "Constrained Wasserstein") {
                                              if(private$cwass.targ == "SBW"){
                                                tf.est <- private$temp.output[["SBW"]]$estimand == estimand
                                                tf.delta <- private$temp.output[["SBW"]]$delta == opts$delta
                                                if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                                select.rows <- which(tf.delta & tf.est)[1]
                                                temp.wass <- private$temp.output[["SBW"]][select.rows, "wasserstein"][[1]]
                                                idx <- which(temp.wass$metric == opts$metric & 
                                                               temp.wass$ground_p == opts$ground_p &
                                                               temp.wass$wass_p == opts$wass_p )
                                                # delta <- temp.wass$dist[[idx]]
                                              } else if (private$cwass.targ == "RKHS") {
                                                tf.est <- private$temp.output[["RKHS"]]$estimand == estimand
                                                # tf.delta <- private$temp.output[["RKHS"]]$delta == opts$delta
                                                # if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                                select.rows <- which(tf.est)[1]
                                                temp.wass <- private$temp.output[["RKHS"]][select.rows, "wasserstein"][[1]]
                                                idx <- which(temp.wass$metric == opts$metric & 
                                                               temp.wass$ground_p == opts$ground_p &
                                                               temp.wass$wass_p == opts$wass_p )
                                                # delta <- temp.wass$dist[[idx]]
                                              } 
                                              # else if (private$cwass.targ == "RKHS.dose") {
                                              #   tf.est <- private$temp.output[["RKHS.dose"]]$estimand == estimand
                                              #   tf.delta <- private$temp.output[["RKHS.dose"]]$delta == opts$delta
                                              #   if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                              #   select.rows <- which(tf.delta & tf.est)[1]
                                              #   temp.wass <- private$temp.output[["RKHS.dose"]][select.rows, "wasserstein"][[1]]
                                              #   idx <- which(temp.wass$metric == opts$metric & 
                                              #                  temp.wass$ground_p == opts$ground_p &
                                              #                  temp.wass$wass_p == opts$wass_p )
                                              # }
                                              delta <- temp.wass$dist[[idx]]
                                              # if(is.null(delta)) delta <- NA
                                              return(delta)
                                            } else if (method == "Wasserstein" | method == "RKHS") {
                                              return(NA)
                                            } else {
                                              return(opts$delta)
                                            }
                                          },
                                          get_cost = function(opt, method, estimand) {
                                            if(grepl("Wasserstein", method)){
                                              cost <- if(is.null(opt$metric) | is.null(opt$ground_p)){
                                                NULL
                                                # } else if (!is.null(opt$theta) & !is.null(opt$gamma)) {
                                                #   private$kernel
                                              } else {
                                                if(opt$metric != "RKHS") {
                                                  private$costs[[opt$metric]][[as.character(opt$ground_p)]] 
                                                } else {
                                                  private$costs[["RKHS"]][[estimand]]
                                                }
                                              }
                                            } else if (method == "RKHS"){
                                              # cost <- private$costs[["RKHS"]][[estimand]]
                                              cost <- NULL
                                            } else {
                                              cost <- NULL
                                            }
                                            return(cost)
                                          },
                                          get_power = function(opt) {
                                            if (is.null(opt$wass_p)) {
                                              return(opt$rkhs_p)
                                            } else {
                                              return(opt$wass_p)
                                            }
                                          },
                                          # kernel.setup = function() {
                                          #   private$kernel <- kernel_calculation(private$x, private$z, p = )
                                          # }
                                          method.setup = function() {
                                            private$method.lookup <- data.frame(method = private$method)
                                            nrows <- nrow(private$method.lookup)
                                            private$method.lookup$weight.fun <- sapply(private$method, function(mm) switch(mm,
                                                                                                                           Logistic = "calc_weight",
                                                                                                                           NNM = "calc_weight",
                                                                                                                           SBW = "calc_weight",
                                                                                                                           RKHS = "calc_weight",
                                                                                                                           'Constrained Wasserstein' = "calc_weight",    
                                                                                                                           Wasserstein = "calc_weight",
                                                                                                                           NA_character_
                                            ))
                                            private$method.lookup$estimand <- lapply(1:nrows, function(i) private$estimand)
                                            private$method.lookup$outcome.model <- if(length(private$outcome.model) == 1) {
                                              lapply(1:nrows, function(i) private$outcome.model)
                                            } else {
                                              lapply(private$outcome.model, function(i) i)
                                            }
                                            private$method.lookup$outcome.model[private$method == "gp"] <-  list("gp_pred")
                                            private$method.lookup$outcome.formula <- lapply(1:nrows, function(i) private$outcome.formula) 
                                            private$method.lookup$model.aug <- lapply(1:nrows, function(i) switch(private$model.augmentation, 
                                                                                                                  "both" = c(FALSE, TRUE),
                                                                                                                  "yes" = TRUE,
                                                                                                                  "no" = FALSE))
                                            private$method.lookup$match <- lapply(1:nrows, function(i) switch(private$match, 
                                                                                                              "both" = c(FALSE, TRUE),
                                                                                                              "yes" = TRUE,
                                                                                                              "no" = FALSE))
                                            private$method.lookup$solver <- sapply(private$method, function(mm) switch(mm,
                                                                                                                       Logistic = "glm",
                                                                                                                       NNM = NA_character_,
                                                                                                                       SBW = private$solver,
                                                                                                                       RKHS = private$solver,
                                                                                                                       RKHS.dose = private$solver,
                                                                                                                       'Constrained Wasserstein' = private$solver,    
                                                                                                                       Wasserstein = private$solver,
                                                                                                                       NA_character_
                                            ))
                                            sdm <- private$standardized.difference.means
                                            lambdas <- private$RKHS$lambdas
                                            theta = list(private$RKHS$theta)
                                            gamma = list(private$RKHS$gamma)
                                            rkhs_p = private$RKHS$p
                                            sigma_2 = private$RKHS$sigma_2
                                            kernel = private$RKHS$kernel
                                            grid.search = private$grid.search
                                            if(isTRUE(private$grid.search)) sdm <- NA
                                            if(isTRUE(private$grid.search)) lambdas <- NA
                                            if(isTRUE(private$RKHS$opt)) {
                                              theta = NA
                                              gamma = NA
                                              sigma_2 = NA
                                              rkhs_p = NA
                                            }
                                            RKHS_list <- RKHS.dose_list <- list(delta = lambdas,
                                                                                theta = theta,
                                                                                gamma = gamma,
                                                                                rkhs_p = rkhs_p,
                                                                                metric = private$RKHS$metric,
                                                                                kernel = kernel,
                                                                                sigma_2 = sigma_2,
                                                                                grid.search = private$grid.search,
                                                                                opt = private$RKHS$opt,
                                                                                opt.method = private$RKHS$opt.method)
                                            RKHS_list$delta <- NULL
                                            wass_list <- Cwass_list <- list(metric = private$metric,
                                                                            ground_p = private$ground_powers,
                                                                            wass_p = private$wass_powers,
                                                                            std_diff = switch(private$cwass.targ,
                                                                                              "SBW" = sdm, 
                                                                                              "RKHS" = NA,
                                                                                              "RKHS.dose" = lambdas),
                                                                            # RKHS.metric = private$RKHS$metric,
                                                                            delta = NA
                                            )
                                            # if(!any(private$metric == "RKHS")) wass_list$RKHS.metric <- Cwass_list$RKHS.metric <- NULL
                                            wass_list$std_diff <- NA
                                            private$method.lookup$options <- sapply(private$method, function(mm) switch(mm,
                                                                                                                        Logistic = list(delta = private$truncations),
                                                                                                                        NNM = list(delta = NA),
                                                                                                                        SBW = list(grid.search = private$grid.search,
                                                                                                                                   delta = sdm),   
                                                                                                                        RKHS = RKHS_list,
                                                                                                                        RKHS.dose = RKHS.dose_list,
                                                                                                                        'Constrained Wasserstein' = Cwass_list,    
                                                                                                                        Wasserstein = wass_list,
                                                                                                                        gp = list(NA)))
                                            if("Logistic" %in% private$method) private$method.lookup$estimand[private$method.lookup$method == "Logistic"][[1]] <- private$estimand[private$estimand != "feasible"]
                                            private$max.conditions <- private$max.cond.calc()
                                            for(i in private$method) private$method.lookup$options[[i]] <- private$set.opts(private$method.lookup[i == private$method.lookup$method,])
                                          },
                                          max.cond.calc = function() {
                                            dims <- sapply(private$method, function(mm){
                                              cur <- private$method.lookup[private$method.lookup$method == mm,]
                                              lens <- sapply(cur$options[[1]], length)
                                              noptions <- prod(lens[lens!=0])
                                              cur$options <- NULL
                                              dim <- prod( sapply(unlist(cur, recursive = FALSE), length) ) * noptions
                                              return(dim)
                                            })
                                            return(max(dims, na.rm= TRUE))
                                          },
                                          set.opts = function(cur) { #sets options for run in estimate function
                                            lens <- sapply(cur$options[[1]], function(o) if(!is.null(o)) {return(length(o))} else {return(1)})
                                            nrows <- prod(lens)
                                            
                                            iter.list <- if(cur$method == "Wasserstein") {
                                              # private$wass_df[1:3]
                                              lapply(cur$options[[1]], rep, length.out = nrows)
                                            } else if (cur$method == "Constrained Wasserstein") {
                                              cur$options[[1]]$delta <- cur$options[[1]]$std_diff
                                              # temp <- list()
                                              # if(!is.null(temp.delta)) {
                                              #   temp$delta <- rep_len(temp.delta, length(temp[[1]]))
                                              # } else {
                                              #   temp["delta"] <- NA
                                              # }
                                              # temp <- lapply(private$wass_df[1:3], function(o) if(!is.null(temp.delta)) {rep(o, each = lens["std_diff"])} else {o})
                                              # temp
                                              lapply(cur$options[[1]], rep, length.out = nrows)
                                            } else {
                                              lapply(cur$options[[1]], rep, length.out = nrows)
                                            }
                                            opts <- vector("list", nrows)
                                            for(i in 1:nrows) {
                                              opts[[i]] <- lapply(iter.list, function(o) o[[i]])
                                            }
                                            return(opts)
                                          },
                                          update = function() {
                                            private$simulator$gen_data()
                                            private$cost.setup()
                                            # private$kernel.setup()
                                            if(private$verbose) message("  Method: ",appendLF = FALSE)
                                            for(mm in private$method) {
                                              if(private$verbose) {
                                                if(mm == private$method[length(private$method)]) {
                                                  message(mm)
                                                } else {
                                                  message(mm,", ",appendLF = FALSE)
                                                }
                                              }
                                              if(mm %in% private$weight_based_methods) {
                                                private$estimate(mm) #updates private$estimate(mm)
                                              } else {
                                                private$model_estimate(mm) #updates private$estimate(mm)
                                              }
                                              
                                              private$temp.output[[mm]] <- private$output.dt[!is.na(private$output.dt$estimate),] 
                                              data.table::set(private$output.dt, i = NULL, j = "estimate", value = NA_real_)
                                            }
                                            
                                            return(data.table::rbindlist(private$temp.output))
                                          },
                                          wass.calc = function(iter, estimand) {
                                            wass.df <- private$wass_df
                                            wass_iter <- 1L
                                            for(metric in private$metric) {
                                              if(metric == "RKHS") next
                                              for(ground_p in private$ground_powers) {
                                                for(wass_p in private$wass_powers) {
                                                  wass.df[["dist"]][[wass_iter]] <- causalOT::wasserstein_p(private$weights[[estimand]], 
                                                                                                            p = wass_p, 
                                                                                                            cost = private$costs[[metric]][[as.character(ground_p)]])
                                                  # names(wass.df[["dist"]][[wass_iter]]) <-
                                                  wass_iter <- wass_iter + 1L
                                                }
                                              }
                                            }
                                            if("RKHS" %in% private$metric) {
                                              temp.cost <- temp.wass <- cost.list <- cost.0 <- cost.1 <- NULL
                                              # for(j in private$RKHS$metric){
                                              
                                              for(wass_p in private$wass_powers) {
                                                if(estimand == "cATE"|estimand == "ATE" | estimand == "feasible") {
                                                  cost.list <- private$costs[["RKHS"]][[estimand]]#[[j]]
                                                  # cost.0 <- cost.list[[1]][,private$simulator$get_z() == 1]
                                                  # cost.1 <- cost.list[[2]][,private$simulator$get_z() == 0]
                                                  # temp.cost <- cost.0/2 + t(cost.1)/2
                                                  w0 <- list(w0 = private$weights[[estimand]]$w0, w1 = renormalize(rep(1,ncol(cost.list[[1]]))),
                                                             gamma = NULL,
                                                             estimand = "ATE")
                                                  w1 <- list(w0 = private$weights[[estimand]]$w1, w1 = renormalize(rep(1,ncol(cost.list[[2]]))),
                                                             gamma = NULL,
                                                             estimand = "ATE")
                                                  class(w1) <- class(w0) <- "causalWeights"
                                                  temp.wass <- c(causalOT::wasserstein_p(w0, 
                                                                                         p = wass_p, 
                                                                                         cost = cost.list[[1]]),
                                                                 causalOT::wasserstein_p(w1, 
                                                                                         p = wass_p, 
                                                                                         cost = cost.list[[2]]))
                                                } else {
                                                  temp.wass <- causalOT::wasserstein_p(private$weights[[estimand]], 
                                                                                       p = wass_p, 
                                                                                       cost = private$costs[["RKHS"]][[estimand]]#[[j]]
                                                  )
                                                }
                                                
                                                
                                                for(ground_p in private$ground_powers) {
                                                  wass.df[["dist"]][[wass_iter]] <- temp.wass
                                                  wass_iter <- wass_iter + 1L
                                                }
                                                # }
                                              }
                                            }
                                            data.table::set(private$output.dt, i = iter, j = "wasserstein" , value = list(wass.df))
                                          },
                                          weight.calc = function(cur, estimand, 
                                                                 solver, delta,
                                                                 cost = NULL, 
                                                                 p = NULL,
                                                                 grid.search = FALSE,
                                                                 opt.hyperparam = TRUE,
                                                                 opt.method = c("stan", "optim", "bayesian.optimization"),
                                                                 metric = metric) {
                                            method <- as.character(cur$method[[1]])
                                            if(grid.search & method == "SBW") delta <- private$standardized.difference.means
                                            if(grid.search & method == "RKHS.dose") delta <- private$RKHS$lambdas
                                            if(method == "RKHS"){
                                              if(opt.hyperparam & !is.null(private$RKHS.opt[[estimand]]) ) {
                                                opt.hyperparam <- FALSE
                                                theta <- private$RKHS.opt[[estimand]]$theta
                                                gamma <- private$RKHS.opt[[estimand]]$gamma
                                                sigma_2 <- private$RKHS.opt[[estimand]]$sigma_2
                                                lambda <- 0
                                                power <- private$RKHS.opt[[estimand]]$p
                                                # kernel <- private$RKHS$kernel
                                              } else {
                                                theta <- private$RKHS$theta
                                                gamma <- private$RKHS$gamma
                                                sigma_2 <- private$RKHS$sigma_2
                                                lambda <- delta
                                                power <- p[[1]]
                                                # kernel <- private$RKHS$kernel
                                              }
                                            } else {
                                              theta <- private$RKHS$theta
                                              gamma <- private$RKHS$gamma
                                              sigma_2 <- private$RKHS$sigma_2
                                              lambda <- delta
                                              power <- p[[1]]
                                              # kernel <- private$RKHS$kernel
                                            }
                                            if (estimand != "cATE") {
                                              private$weights[[estimand]]<- 
                                                calc_weight(private$simulator,  
                                                            constraint = delta,
                                                            estimand = estimand, 
                                                            method = method,
                                                            cost = cost, p = power,
                                                            transport.matrix = FALSE,
                                                            solver = solver,
                                                            opt.hyperparam = opt.hyperparam,
                                                            opt.method = opt.method,
                                                            grid.search = grid.search,
                                                            grid = delta,
                                                            kernel = private$RKHS$kernel,
                                                            rkhs.args = private$RKHS.opt[[estimand]],
                                                            theta =  theta,
                                                            gamma = gamma,
                                                            sigma_2 = sigma_2,
                                                            lambda = lambda,
                                                            iter = if(is.null(private$RKHS$iter)) 2000 else private$RKHS$iter,
                                                            maxit = if(is.null(private$RKHS$iter)) 2000 else private$RKHS$iter,
                                                            metric = metric)
                                            } else {
                                              private$weights[[estimand]] <- 
                                                convert_ATE(private$weights[["ATT"]], 
                                                            private$weights[["ATC"]],
                                                            transport.matrix = FALSE,
                                                            cost = cost, p = p[[1]])
                                            }
                                          }
                           )
                           
  )
}

#options class
# {
  # setClass("simOptions", slots = c(method = "character",
  #                                  weight.fun = "character",
  #                                  estimand = "character",
  #                                  outcome.model = "character",
  #                                  outcome.formula = "character",
  #                                  model.aug = "logical",
  #                                  match = "logical",
  #                                  solver = "character",
  #                                  options = "character"
  #                                  ),
  #                    prototype = c(method = NA_character_,
  #                                  weight.fun = NA_character_,
  #                                  estimand = NA_character_,
  #                                  outcome.model = NA_character_,
  #                                  outcome.formula = NA_character_,
  #                                  model.aug = NA,
  #                                  match = NA,
  #                                  solver = NA_character_,
  #                                  options = NA_character_))
#   setClass("logisticOptions", slots =
#              c(delta = "numeric"),
#            )
#   
#   setClass("nnmOptions", slots =
#              c(delta = "numeric")
#   )
#   setClass("rkhsOptions", slots =
#              c(delta = "numeric")
#   )
#   setClass("rkhsDoseOptions", slots =
#              c(delta = "numeric")
#   )
#   setClass("constrainedWassersteinOptions", slots =
#              c(delta = "numeric")
#   )
#   setClass("wassersteinOptions", slots =
#              c(delta = "numeric")
#   )
#   setClass("gpOptions", slots =
#              c(delta = "numeric")
#   )
#   # switch(mm,
#   #        Logistic = list(delta = private$truncations),
#   #        NNM = list(NA),
#   #        SBW = list(grid.search = private$grid.search,
#   #                   delta = sdm),   
#   #        RKHS = RKHS_list,
#   #        RKHS.dose = RKHS.dose_list,
#   #        'Constrained Wasserstein' = Cwass_list,    
#   #        Wasserstein = wass_list,
#   #        gp = list(NA)))
#   # if("Logistic" %in% private$method) private$method.lookup$estimand[private$method.lookup$method == "Logistic"][[1]] <- private$estimand[private$estimand != "feasible"]
#   # private$max.conditions <- private$max.cond.calc()
#   # for(i in private$method) private$method.lookup$options[[i]] <- private$set.opts(private$method.lookup[i == private$method.lookup$method,])
#   
# }
