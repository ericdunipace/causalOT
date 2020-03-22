SimHolder <- R6::R6Class("SimHolder",
                       public = list(run = function() {
                                     private$output <- vector("list", private$nsim)
                                     
                                     for(iter in 1:private$nsim) {
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
                                       exclude <- c("wasserstein", "ess.frac")
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
                                     get.wass = function(output) {
                                       wass.dt <- data.table::rbindlist(lapply(output[,"wasserstein"], function(o) as.data.frame(o)), idcol = TRUE)
                                       check.rep <- table(wass.dt[,".id"])
                                       stopifnot(all(check.rep == check.rep[1]))
                                       nrep <- check.rep[1]
                                       
                                       exclude <- c("wasserstein", "estimate", "ESS.frac")
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
                                                           RKHS = list(lambdas = NULL, theta = NULL, gamma = NULL, p = NULL,
                                                                       metrics = "malahanobis",
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
                                                           metrics = "Lp") {
                                       private$nsim <- nsim

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
                                                                   c("Lp", "mahalanobis"), several.ok = TRUE)
                                       if(!is.null(RKHS) & !missing(RKHS)) {
                                         if(!is.list(RKHS)) RKHS <- list(RKHS)
                                         private$RKHS <- list()
                                         
                                         if(is.null(RKHS$lambdas)) RKHS$lambdas <- seq(0,100, length.out = 11)
                                         private$RKHS$lambdas <-  RKHS$lambdas
                                         
                                         if(is.null(RKHS$theta)) RKHS$theta <- as.double(c(1,1))
                                         private$RKHS$theta <- RKHS$theta
                                         
                                         if(is.null(RKHS$gamma)) RKHS$gamma <- as.double(c(1,1))
                                         private$RKHS$gamma <- RKHS$gamma
                                         
                                         if(is.null(RKHS$p)) RKHS$p <- as.double(2)
                                         private$RKHS$p <- RKHS$p
                                         
                                         if(is.null(RKHS$sigma_2)) RKHS$sigma_2 <- as.double(1)
                                         private$RKHS$sigma_2  <- RKHS$sigma_2 
                                         
                                         if(is.null(RKHS$opt)) RKHS$opt <- TRUE
                                         private$RKHS$opt <- RKHS$opt
                                         
                                         if(is.null(RKHS$opt.method)) RKHS$opt.method <- "stan"
                                         private$RKHS$opt.method <- RKHS$opt.method
                                         if(!is.null(RKHS$iter)) {
                                           private$RKHS$iter <- RKHS$iter
                                         }
                                         
                                         if(is.null(RKHS$metrics)) RKHS$metrics <- private$metric
                                         private$RKHS$metric <- RKHS$metrics
                                       } else {
                                         private$RKHS <- list()
                                         private$RKHS$lambdas <-  as.double(seq(0,100, length.out = 11))
                                         private$RKHS$theta <- as.double(c(1,1))
                                         private$RKHS$gamma <- as.double(c(1,1))
                                         private$RKHS$p <- as.double(2)
                                         private$RKHS$metric <- private$metric
                                         private$RKHS$sigma_2 <- as.double(1)
                                         private$RKHS$opt <- TRUE
                                         private$RKHS$opt.method <- "stan"
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
                                       private$estimand <- options$estimates
                                       if(!private$calculate.feasible) private$estimand <- private$estimand[private$estimand != "feasible"]
                                       private$method <- options$weights
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
                                       private$wass_df$dist <- numeric(nrows)
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
                                                            estimate = rep(NA_real_, private$max.conditions)
                                                            )
                                       private$weights <- sapply(private$estimand, function(i) {NULL}, simplify = NULL)
                                     }
                                     ),
                       private = list(calculate.feasible = "logical",
                                      costs = "list",
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
                                      simulator = "DataSim",
                                      solver = "character",
                                      standardized.difference.means = "vector",
                                      temp.output = "list",
                                      truncations = "vector",
                                      wass_df = "list",
                                      wass_powers = "vector",
                                      weights = "list",
                                      cost.fun = function(X, Y, ground_p, direction, metric) {
                                        return(switch(metric,
                                                "Lp" = causalOT::cost_calc_lp(X, Y, ground_p, direction),
                                    "mahalanobis" = causalOT::cost_mahalanobis(X, Y, ground_p, direction)))
                                        },
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
                                            for(j in as.character(private$ground_powers)) {
                                              private$costs[[i]][[j]] <- private$cost.fun(x0, x1,
                                                                                        ground_p = as.numeric(j), 
                                                                                        "rowwise", 
                                                                                        i)
                                            }
                                        }
                                      },
                                      check.skip = function(weights) {
                                        skip <- FALSE
                                        if(isTRUE(is.na(weights$w0)) | isTRUE(length(weights$w0) == 0) | isTRUE(is.null(weights$w0))) {
                                          skip <- TRUE
                                        }
                                        if(isTRUE(is.na(weights$w1)) | isTRUE(length(weights$w1) == 0) | isTRUE(is.null(weights$w1))) {
                                          skip <- TRUE
                                        }
                                        return(skip)
                                      },
                                      estimate = function(method) {
                                        cur <- private$method.lookup[private$method.lookup$method == method,]
                                        ess.frac <- list()
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
                                              if ( isTRUE(method == "RKHS" | method == "RKHS.dose") & isFALSE(est == "ATE" | est == "feasible") ) next
                                              private$weight.calc(cur = cur, 
                                                                  estimand = est, 
                                                                  solver = solver,
                                                                  delta = delta,
                                                                  cost = private$get_cost(o), 
                                                                  p = private$get_power(o),
                                                                  grid.search = isTRUE(o$grid.search),
                                                                  opt.hyperparam = isTRUE(o$opt),
                                                                  opt.method = o$opt.method
                                                                  )
                                              if(private$check.skip(private$weights[[est]])) next
                                              ess.frac <- list(ESS(private$weights[[est]])/ns)
                                              for (mods in cur$outcome.model[[1]]) {
                                                for (aug in cur$model.aug[[1]]) {
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
                                                      data.table::set(private$output.dt, i = iter, j = "estimate", 
                                                                      value = outcome_model(private$simulator, 
                                                                                            formula = cur$outcome.formula[[1]][[aug + 1]],
                                                                                            weights = private$weights[[est]],
                                                                                            hajek = TRUE,
                                                                                            doubly.robust = aug,
                                                                                            matched = match,
                                                                                            target = est,
                                                                                            model = match.fun(private$outcome.model[[1]])))
                                                      iter <- iter + 1L
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                        }
                                      },
                                      get_delta = function(opts, estimand, method) {
                                        if(method == "Constrained Wasserstein") {
                                          tf.est <- private$temp.output[["SBW"]]$estimand == estimand
                                          tf.delta <- private$temp.output[["SBW"]]$delta == opts$delta
                                          tf.delta <- if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) rep(TRUE, length(tf.delta))
                                          select.rows <- which(tf.delta & tf.est)[1]
                                          temp.wass <- private$temp.output[["SBW"]][select.rows, "wasserstein"][[1]]
                                          idx <- which(temp.wass$metric == opts$metric & 
                                                         temp.wass$ground_p == opts$ground_p &
                                                         temp.wass$wass_p == opts$wass_p )
                                          delta <- temp.wass$dist[idx]
                                          # if(is.null(delta)) delta <- NA
                                          return(delta)
                                        } else if (method == "Wasserstein") {
                                          return(NA)
                                        } else {
                                          return(opts$delta)
                                        }
                                      },
                                      get_cost = function(opt) {
                                        cost <- if(is.null(opt$metric) | is.null(opt$ground_p)){
                                          NULL
                                        # } else if (!is.null(opt$theta) & !is.null(opt$gamma)) {
                                        #   private$kernel
                                        } else {
                                          private$costs[[opt$metric]][[as.character(opt$ground_p)]] 
                                        }
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
                                                                                                                       SBW = "calc_weight",
                                                                                                                       RKHS = "calc_weight",
                                                                                                                       'Constrained Wasserstein' = "calc_weight",    
                                                                                                                       Wasserstein = "calc_weight"
                                        ))
                                        private$method.lookup$estimand <- lapply(1:nrows, function(i) private$estimand)
                                        private$method.lookup$outcome.model <- lapply(1:nrows, function(i) private$outcome.model) 
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
                                                                                                                   SBW = private$solver,
                                                                                                                   RKHS = private$solver,
                                                                                                                   RKHS.dose = private$solver,
                                                                                                                   'Constrained Wasserstein' = private$solver,    
                                                                                                                   Wasserstein = private$solver
                                        ))
                                        sdm <- private$standardized.difference.means
                                        lambdas <- private$RKHS$lambdas
                                        theta = list(private$RKHS$theta)
                                        gamma = list(private$RKHS$gamma)
                                        rkhs_p = private$RKHS$p
                                        sigma_2 = private$RKHS$sigma_2
                                        grid.search = private$grid.search
                                        if(isTRUE(private$grid.search)) sdm <- NA
                                        if(isTRUE(private$grid.search)) lambdas <- NA
                                        if(isTRUE(private$RKHS$opt)) {
                                          theta = NA
                                          gamma = NA
                                          sigma_2 = NA
                                        }
                                        RKHS_list <- list(delta = lambdas,
                                                          theta = theta,
                                                          gamma = gamma,
                                                          rkhs_p = rkhs_p,
                                                          metric = private$RKHS$metric,
                                                          sigma_2 = sigma_2,
                                                          grid.search = private$grid.search,
                                                          opt = private$RKHS$opt,
                                                          opt.method = private$RKHS$opt.method)
                                        wass_list <- Cwass_list <- list(metric = private$metric,
                                                          ground_p = private$ground_powers,
                                                          wass_p = private$wass_powers,
                                                          std_diff = sdm,
                                                          delta = NULL
                                        )
                                        wass_list$std_diff <- NULL
                                        private$method.lookup$options <- sapply(private$method, function(mm) switch(mm,
                                                                                                                    Logistic = list(delta = private$truncations),
                                                                                                                    SBW = list(grid.search = private$grid.search,
                                                                                                                               delta = sdm),   
                                                                                                                    RKHS = RKHS_list,
                                                                                                                    RKHS.dose = RKHS_list,
                                                                                                                    'Constrained Wasserstein' = Cwass_list,    
                                                                                                                    Wasserstein = wass_list))
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
                                      set.opts = function(cur) {
                                        lens <- sapply(cur$options[[1]], function(o) if(!is.null(o)) {return(length(o))} else {return(1)})
                                        nrows <- prod(lens)
                                        opts <- vector("list", nrows)
                                        iter.list <- if(cur$method == "Wasserstein") {
                                          private$wass_df[1:3]
                                        } else if (cur$method == "Constrained Wasserstein") {
                                          temp.delta <- cur$options[[1]]$std_diff
                                          temp <- lapply(private$wass_df[1:3], function(o) if(!is.null(temp.delta)) {rep(o, each = lens["std_diff"])} else {o})
                                          if(!is.null(temp.delta)) {
                                            temp$delta <- rep_len(temp.delta, length(temp[[1]]))
                                          } else {
                                            temp["delta"] <- list(NULL)
                                          }
                                          temp
                                        } else {
                                          lapply(cur$options[[1]], rep, length.out = nrows)
                                        }
                                        for(i in 1:nrows) {
                                          opts[[i]] <- lapply(iter.list, function(o) o[[i]])
                                        }
                                        return(opts)
                                      },
                                      update = function() {
                                        private$simulator$gen_data()
                                        private$cost.setup()
                                        # private$kernel.setup()
                                        for(mm in private$method) {
                                          private$estimate(mm) #updates private$estimate(mm)
                                          private$temp.output[[mm]] <- private$output.dt[!is.na(private$output.dt$estimate),] 
                                          data.table::set(private$output.dt, i = NULL, j = "estimate", value = NA_real_)
                                        }
                                      
                                        return(data.table::rbindlist(private$temp.output))
                                      },
                                      wass.calc = function(iter, estimand) {
                                        wass.df <- private$wass_df
                                        wass_iter <- 1L
                                        for(metric in private$metric) {
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
                                        data.table::set(private$output.dt, i = iter, j = "wasserstein" , value = list(wass.df))
                                    },
                                      weight.calc = function(cur, estimand, 
                                                             solver, delta,
                                                             cost = NULL, 
                                                             p = NULL,
                                                             grid.search = FALSE,
                                                             opt.hyperparam = TRUE,
                                                             opt.method = c("stan", "optim", "bayesian.optimization")) {
                                        method <- as.character(cur$method[[1]])
                                        if(grid.search & method == "SBW") delta <- private$standardized.difference.means
                                        if(grid.search & method == "RKHS") delta <- private$RKHS$lambdas
                                        if ((estimand != "ATE" & method != "RKHS" & method != "RKHS.dose") | 
                                            ((estimand == "ATE" | estimand == "feasible") & method == "RKHS") | 
                                            ((estimand == "ATE" | estimand == "feasible") & method == "RKHS.dose") | 
                                            method == "Logistic") {
                                          private$weights[[estimand]]<- 
                                            calc_weight(private$simulator,  constraint = delta,
                                                        estimate = estimand, method = method,
                                                        cost = cost, p = p[[1]],
                                                        transport.matrix = FALSE,
                                                        solver = solver,
                                                        opt.hyperparam = opt.hyperparam,
                                                        opt.method = opt.method,
                                                        grid.search = grid.search,
                                                        grid = delta,
                                                        theta = private$RKHS$theta,
                                                        gamma = private$RKHS$gamma,
                                                        iter = if(is.null(private$RKHS$iter)) 2000 else private$RKHS$iter,
                                                        maxit = if(is.null(private$RKHS$iter)) 2000 else private$RKHS$iter)
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
