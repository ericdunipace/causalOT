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
                             wass.dt <- data.table::rbindlist(lapply(output[,"wasserstein"], function(o) o), idcol = TRUE, fill = TRUE)
                             # dist <- unlist(lapply(output[,"wasserstein"], function(o) as.data.frame(o[4])), recursive = FALSE)
                             # wass.dt <- data.table::set(wass.dt, i = NULL, j = "dist", value = dist)
                             check.rep <- table(wass.dt[,".id"])
                             stopifnot(all(check.rep == check.rep[1]))
                             nrep <- check.rep[1]
                             
                             exclude <- c("wasserstein", "estimate", "ESS.frac","psis.k","psis.ess.frac")
                             exclude.cn <- (-which(colnames(output) %in% exclude))
                             wass.output <- output[rep(1:nrow(output), each = nrep), exclude.cn]
                             data.table::setDT(wass.output)
                             cn <- colnames(wass.dt)
                             cn <- cn[cn!= ".id"]
                             
                             for(j in cn) data.table::set(wass.output, i = NULL, j = j, value = wass.dt[,j])
                             
                             return(wass.output)
                           },
                           initialize = function(nsim = 100,
                                                 dataSim,
                                                 methods = NULL,
                                                 grid.search = TRUE,
                                                 truncations = NULL,
                                                 standardized.difference.means = NULL,
                                                 
                                                 estimands = c("ATT","ATC","cATE","ATE"),
                                                 RKHS = list(lambdas = NULL, theta = NULL, gamma = NULL, p = NULL,
                                                             metric = "malahanobis",
                                                             kernel = "RBF",
                                                             sigma_2 = NULL, opt = NULL, opt.method = NULL),
                                                 Wass = list(method = "networkflow",
                                                             niter = 0,
                                                             epsilon = 10, #0.05,
                                                             powers = 2,
                                                             # ground_powers = 2,
                                                             metrics = "Lp",
                                                             constrained.wasserstein.target = c("RKHS", "SBW"),
                                                             wasserstein.distance.constraints = NULL,
                                                             add.joint = TRUE,
                                                             add.margins = FALSE,
                                                             eval.method = "cross.validation",
                                                             cross.val.replicates = 10,
                                                             confidence.interval = FALSE,
                                                             add.divergence = FALSE),
                                                 outcome.model = "lm",
                                                 outcome.formula = list(none = NULL,
                                                                        augmentation = NULL),
                                                 propensity.formula = list(Logistic = NULL,
                                                                        SBW = NULL,
                                                                        CBPS = NULL,
                                                                        "Constrained Wasserstein" = NULL,
                                                                        "Wasserstein" = NULL),
                                                 model.augmentation = "both",
                                                 match = "both",
                                                 split.models = "both",
                                                 calculate.feasible = FALSE,
                                                 solver = "gurobi",
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
                             if(is.null(Wass)) {
                               Wass <- list(method = "sinkhorn_geom",
                                           niter = 1e3,
                                           epsilon = 10, #0.05,
                                           powers = 2,
                                           # ground_powers = 2,
                                           metrics = "mahalanobis",
                                           constrained.wasserstein.target = c("SBW"),
                                           wasserstein.distance.constraints = NULL,
                                           add.joint = TRUE,
                                           add.margins = FALSE,
                                           penalty = "L2",
                                           joint.mapping = FALSE,
                                           neg.weights = FALSE,
                                           eval.method = "bootstrap",
                                           cross.val.replicates = 10,
                                           confidence.interval = FALSE,
                                           add.divergence = FALSE)
                             }
                             private$wass.opt <- list()
                             if (!is.null(Wass$metrics)) {
                               private$metric <- match.arg(Wass$metrics,
                                                           dist.metrics(), several.ok = TRUE)
                             } else {
                               private$metric <- "mahalanobis"
                             }
                             
                             if (!is.null(Wass$wass_powers)) {
                               private$wass_powers <- as.numeric(Wass$wass_powers)
                             } else {
                               private$wass_powers <- 2.0
                             }
                             
                             if(!is.null(Wass$wass_powers)) {
                               private$ground_powers  <- private$wass_powers #as.numeric(Wass$ground_powers)
                             } else {
                               private$ground_powers <- 2.0
                             }
                             
                             if(!is.null(Wass$method)) {
                               private$wass.opt$method <- match.arg(Wass$method, choices = c(approxOT::transport_options(),"sinkhorn_geom"))
                             } else {
                               private$wass.opt$method <- "sinkhorn_geom"
                             }
                             
                             
                             if(!is.null(Wass$niter)) {
                               private$wass.opt$niter <- as.double(Wass$niter)
                             } else {
                               private$wass.opt$niter <- switch(private$wass.opt$method ,
                                                                 "networkflow" = 0,
                                                                 "exact" = 0,
                                                                 1000)
                             }
                             
                             if (!is.null(Wass$epsilon)) {
                               private$wass.opt$epsilon <- as.double(Wass$epsilon)
                             } else {
                               private$wass.opt$epsilon <-  10 #0.05
                             }
                             if (!is.null(Wass$constrained.wasserstein.target)) {
                               private$cwass.targ <- match.arg(Wass$constrained.wasserstein.target, c("RKHS", "SBW"))
                             } else {
                               private$cwass.targ <- "SBW"
                             }
                             if (!is.null(Wass$wasserstein.distance.constraints)) {
                               private$wasserstein.distance.constraints <- as.double(Wass$wasserstein.distance.constraints)
                             } else {
                               private$wasserstein.distance.constraints <- NA
                             }
                             if (!is.null(Wass$add.joint)) {
                               private$wass.opt$add.joint <- isTRUE(Wass$add.joint)
                             } else {
                               private$wass.opt$add.joint <- TRUE
                             }
                             
                             if (!is.null(Wass$add.margins)) {
                               private$wass.opt$add.margins <- sapply(Wass$add.margins, isTRUE)
                             } else {
                               private$wass.opt$add.margins <- FALSE
                             }
                             
                             if (!is.null(Wass$joint.mapping)) {
                               private$wass.opt$joint.mapping <- sapply(Wass$joint.mapping, isTRUE)
                             } else {
                               private$wass.opt$joint.mapping <- c(FALSE)
                             }
                             
                             if (!is.null(Wass$penalty)) {
                               private$wass.opt$penalty <- match.arg(Wass$penalty, c("L2", "variance","entropy", "none"), several.ok = TRUE)
                             } else {
                               private$wass.opt$penalty <- "L2"
                             }
                             if (!is.null(Wass$neg.weights)) {
                               private$wass.opt$neg.weights <- sapply(Wass$neg.weights, isTRUE)
                             } else {
                               private$wass.opt$neg.weights <- FALSE
                             }
                             
                             if (!is.null(Wass$eval.method)) {
                               private$wass.opt$eval.method <- match.arg(Wass$eval.method, c("bootstrap", "cross.validation"))
                             } else {
                               private$wass.opt$eval.method <- "bootstrap"
                             }
                             
                             if (!is.null(Wass$cross.val.replicates)) {
                               private$wass.opt$cross.val.replicates <- Wass$cross.val.replicates
                             } else {
                               private$wass.opt$cross.val.replicates <- 10
                             }
                             
                             if (!is.null(Wass$confidence.interval) ) {
                               private$wass.opt$confidence.interval <- match.arg(Wass$confidence.interval, c("asymptotic", "bootstrap"))
                             } else {
                               private$wass.opt$confidence.interval <- FALSE
                             }
                             
                             if (!is.null(Wass$add.divergence) ) {
                               private$wass.opt$add.divergence <- sapply(Wass$add.divergence, isTRUE)
                             } else {
                               private$wass.opt$add.divergence <- FALSE
                             }
                               
                             private$SBW.balconst <- list("ATT" = 0.1,
                                                          "ATC" = 0.1,
                                                          "cATE" = 0.1,
                                                          "ATE" = 0.1)
                             # private$metric <- match.arg(Wass$metrics,
                             #                             c("Lp", "mahalanobis", "RKHS"), several.ok = TRUE)
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
                               
                               if(is.null(RKHS$algorithm)) RKHS$algorithm <- "LBFGS"
                               private$RKHS$algorithm <- RKHS$algorithm
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
                               private$RKHS$algorithm <- "LBFGS"
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
                             
                             
                             
                             if(is.null(propensity.formula) || missing(propensity.formula)) {
                               propensity.formula <- list()
                             } else if (!is.list(propensity.formula) && is.character(propensity.formula)) {
                               propensity.formula <- list(logistic = propensity.formula,
                                    cbps = propensity.formula,
                                    sbw = propensity.formula,
                                    cwass = propensity.formula,
                                    wass = propensity.formula)
                             }
                             private$ps.formula <- list(logistic = unlist(propensity.formula$Logistic),
                                                        cbps = unlist(propensity.formula$CBPS),
                                                        sbw = unlist(propensity.formula$SBW),
                                                        cwass = unlist(propensity.formula$`Constrained Wasserstein`),
                                                        wass = unlist(propensity.formula$Wasserstein))
                             if(is.null(private$ps.formula$logistic)) private$ps.formula$logistic <- "z ~ ."
                             if(is.null(private$ps.formula$cbps)) private$ps.formula$cbps <- "z ~ ."
                             if(is.null(private$ps.formula$sbw)) private$ps.formula$sbw <- "~ . + 0"
                             if(is.null(private$ps.formula$cwass)) private$ps.formula$cwass <- NA #"z ~ . + 0"
                             if(is.null(private$ps.formula$wass)) private$ps.formula$wass <- NA #"z ~ . + 0"
                             
                             private$calculate.feasible <- FALSE #isTRUE(calculate.feasible)
                             options <- get_holder_options()
                             if(!is.null(estimands)) {
                               estimands <- match.arg(estimands, several.ok = TRUE)
                             } else {
                               estimands <- options$estimates
                             }
                             private$estimand <- options$estimates[options$estimates %in% estimands]
                             if(!private$calculate.feasible) private$estimand <- private$estimand[private$estimand != "feasible"]
                             #removing "wasserstein"
                             pot.methods <- c(options$weights[options$weights != "RKHS.dose"],
                                              "gp")
                             if(!is.null(methods)) {
                               private$method <- match.arg(methods, pot.methods,
                                                           several.ok = TRUE)
                             } else {
                               private$method <- pot.methods
                             }
                             if(any(private$method == "COT")) {
                               private$method[private$method == "COT"] <- "Wasserstein"
                               private$method <- unique(private$method)
                             }
                             private$weight_based_methods <- private$method[private$method != "gp"]
                             private$nmethod <- length(private$method)
                             private$model.augmentation <- match.arg(model.augmentation, c("both", "yes", "no"))
                             private$match <- match.arg(match, c("both", "yes", "no"))
                             private$split <- match.arg(split.models, c("both", "yes", "no"))
                             # private$cost.setup()
                             
                             private$solver <- match.arg(solver, c("mosek","gurobi","cplex","lbfgs","quadprog"))
                             # if(is.character(solver)) {
                             #   private$solver <- match.arg(solver, c("mosek","gurobi","cplex","lbfgs","quadprog"))
                             # } else if (is.list(solver)) {
                             #   solver <- lapply(solver, function(sname) match.arg(sname, c("mosek","gurobi","cplex","lbfgs","quadprog")))
                             #   # solver.find <- f.call.list.no.eval("switch", list.args = solver)
                             #   private$solver <- lapply(private$method ,  function(mm) f.call.list(fun = "switch", list.args = c(EXPR = mm, solver.find$envir)))
                             #   names(private$solver) <- private$method
                             # } else {
                             #   stop("argument `solver` must be a character specifying same solver for all methods or a named list specifying the solver for each method")
                             # }
                             
                             private$temp.output <- vector("list", length = private$nmethod)
                             names(private$temp.output) <- private$method
                             
                             n_g <- length(private$ground_powers)
                             n_m <- length(private$metric)
                             n_w <- length(private$wass_powers)
                             nrows <- n_g * n_m * n_w
                             nrows <- n_m * n_w
                             
                             private$wass_df <- list(metric = rep(private$metric, each = n_g * n_w))
                             # private$wass_df$ground_p <- rep_len(rep(private$ground_powers, each = n_w), length.out = nrows)
                             private$wass_df$wass_p <- rep_len(private$wass_powers, length.out = nrows)
                             private$wass_df$dist <- vector("list",nrows)
                             
                             private$method.setup()
                             
                             private$output.dt <- data.table::data.table(
                               method = rep(private$method[1], private$max.conditions),
                               estimand = character(private$max.conditions),
                               model = character(private$max.conditions),
                               model.augmentation = logical(private$max.conditions),
                               match = logical(private$max.conditions),
                               split.model = logical(private$max.conditions),
                               solver = character(private$max.conditions),
                               delta = numeric(private$max.conditions),
                               options = lapply(1:private$max.conditions, function(i) list(NA)), #vector("list", private$max.conditions),
                               wasserstein = lapply(1:private$max.conditions, function(i) private$wass_df),
                               ess.frac = lapply(1:private$max.conditions, function(i) private$wass_df),
                               psis.ess.frac = lapply(1:private$max.conditions, function(i) private$wass_df),
                               psis.k = lapply(1:private$max.conditions, function(i) private$wass_df),
                               confidence.interval = lapply(1:private$max.conditions, function(i) c(NA_real_, NA_real_)),
                               estimate = rep(NA_real_, private$max.conditions),
                               E_Y1 = rep(NA_real_, private$max.conditions),
                               E_Y0 = rep(NA_real_, private$max.conditions)
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
                                          ps.formula = "list",
                                          RKHS = "list",
                                          RKHS.opt = "list",
                                          SBW.balconst = "list",
                                          simulator = "DataSim",
                                          solver = "character",
                                          split = "character",
                                          standardized.difference.means = "vector",
                                          temp.output = "list",
                                          truncations = "vector",
                                          verbose = "logical",
                                          wass_df = "list",
                                          wass_powers = "vector",
                                          wass.opt = "list",
                                          wasserstein.distance.constraints = "vector",
                                          weights = "list",
                                          weight_based_methods = "character",
                                          #   cost.fun = function(X, Y, ground_p, direction, metric) {
                                          #     return(switch(metric,
                                          #             "Lp" = cost_calc_lp(X, Y, ground_p, direction),
                                          # "mahalanobis" = cost_mahalanobis(X, Y, ground_p, direction)))
                                          #     },
                                          cost.setup = function() {
                                            n_g <- length(private$ground_powers)
                                            n_g <- length(private$wass_powers)
                                            n_m <- length(private$metric)
                                            nrows <- n_g * n_m
                                            nrows <- n_m
                                            ns <- private$simulator$get_n()
                                            # private$costs <- vector("list", n_g)
                                            # names(private$costs) <- private$metric
                                            # list(metric = rep(private$metric, each = n_g))
                                            private$costs <- sapply(private$metric, 
                                                                    function(i) 
                                                                      # sapply(as.character(private$ground_powers), 
                                                                      sapply(as.character(private$wass_powers), 
                                                                                       function(j)
                                                                                         sapply(unique(c("ATT", private$estimand)),
                                                                                                  function(k) 
                                                                                                    matrix(0, nrow = ns["n0"], 
                                                                                                           ncol = ns["n1"] ), simplify = FALSE),
                                                                                       simplify = FALSE), simplify = FALSE)
                                            # private$costs$cost <- lapply(1:nrows, function(i) matrix(0, nrow = ns["n0"], ncol = ns["n1"] ))   
                                            x0 <- private$simulator$get_x0()
                                            x1 <- private$simulator$get_x1()
                                            x  <- private$simulator$get_x()
                                            z  <- private$simulator$get_z()
                                            for(i in as.character(private$metric)) {
                                              if(i == "RKHS") next
                                              # for (j in as.character(private$ground_powers)) {
                                              for (j in as.character(private$wass_powers)) {
                                                # function(X, Y, ground_p, direction, metric)
                                                for(k in as.character(private$estimand)) {
                                                  private$costs[[i]][[j]][[k]] <- cost_fun(x = x, z = z,
                                                                                    power = as.numeric(j), 
                                                                                    metric = i,
                                                                                    estimand = k)
                                                }
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
                                                                 estimand = est,
                                                                 algorithm = private$RKHS$algorithm,
                                                                 verbose = private$verbose)
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
                                                private$costs[["RKHS"]][[i]] <- cost_fun(x=x, z= z,
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
                                            if (isTRUE(all(is.na(weights$w0))) | isTRUE(length(weights$w0) == 0) | isTRUE(is.null(weights$w0))) {
                                              skip <- TRUE
                                            }
                                            if (isTRUE(all(is.na(weights$w1))) | isTRUE(length(weights$w1) == 0) | isTRUE(is.null(weights$w1))) {
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
                                            esteff <- ci.out <- NULL
                                            for (solver in cur$solver) {
                                              for (o in cur$options[[1]]) {
                                                if (method == "Wasserstein") {
                                                  if ( isTRUE(o$add.divergence) ) {
                                                    if (isTRUE(o$add.margins)) next
                                                    if (isTRUE(o$joint.mapping)) next
                                                    if (isTRUE(o$penalty != "entropy")) next
                                                    # if (isTRUE(!is.null(o$formula[[1]]) && !is.na(o$formula[[1]]))) next
                                                    # delta <- list(penalty = 1e4) #check how gridsearch handles this
                                                  } else {
                                                    if (isTRUE(o$penalty == "entropy")) solver <- "lbfgs"
                                                  }
                                                }
                                                if (private$verbose && (method == "Wasserstein" | method == "Constrained Wasserstein"  | method == "SCM") ) print(o)
                                                # if (method == "Wasserstein" || method == "Constrained Wasserstein") {
                                                  # if (isTRUE(o$penalty == "entropy") && isTRUE(o$joint.mapping == TRUE) ) next
                                                # }
                                                for (est in cur$estimand[[1]]) {
                                                  delta <- private$get_delta(o, est, method)
                                                  if ( isTRUE(is.null(delta)) ) next
                                                  if ( isTRUE(method == "RKHS" | method == "RKHS.dose" | method == "Wasserstein") & isTRUE(est == "feasible") ) next
                                                  if ( isTRUE(method == "RKHS.dose") & isFALSE(est == "ATE") ) next
                                                  if ( isTRUE(method == "Constrained Wasserstein" | method == "Wasserstein") & isTRUE(est == "feasible")) {
                                                    # if(o$metric == "RKHS") 
                                                      next
                                                  }
                                                  if ( isTRUE(method == "NNM") & isTRUE(est == "feasible")) next
                                                  if ( isTRUE(method == "Constrained Wasserstein") & isTRUE(o$penalty == "none")) next
                                                  if ( isTRUE(o$neg.weights) && isTRUE(o$penalty == "entropy")) next
                                                  
                                                  private$weight.calc(cur = cur, 
                                                                      estimand = est, 
                                                                      solver = solver,
                                                                      delta = delta,
                                                                      # cost = private$get_cost(o, method, est), 
                                                                      p = private$get_power(o),
                                                                      grid.search = isTRUE(o$grid.search),
                                                                      opt.hyperparam = isTRUE(o$opt),
                                                                      opt.method = o$opt.method,
                                                                      metric = o$metric,
                                                                      formula = o$formula[[1]],
                                                                      add.margins = isTRUE(o$add.margins),
                                                                      add.divergence = isTRUE(o$add.divergence),
                                                                      penalty = o$penalty,
                                                                      joint.mapping = isTRUE(o$joint.mapping)
                                                  )
                                                  if (private$check.skip(private$weights[[est]])) next
                                                  if ((o$formula[[1]] == "~. + 0" | o$formula[[1]] == "~ . + 0" |
                                                      o$formula[[1]] == "~.+0" | o$formula[[1]] == "~ .+ 0") && method == "SBW") private$SBW.balconst[[est]] <- private$weights[[est]]$args$constraint 
                                                  ess.frac <- list(ESS(private$weights[[est]])/ns)
                                                  psis.output <- PSIS_diag(private$weights[[est]])
                                                  psis.ess.frac <- list(sapply(psis.output, function(w) w$n_eff)/ns)
                                                  psis.k <- list(lapply(psis.output, function(w) w$pareto_k))
                                                  # if(method == "Constrained Wasserstein" & est == "ATE") if(o$metric=="RKHS") browser()
                                                  for (mods in cur$outcome.model[[1]]) {
                                                    for (aug in cur$model.aug[[1]]) {
                                                      # opt.dist <- list(private$simulator$opt_weight_dist(weight = private$weights[[est]], 
                                                      #                                                    estimand = est, augment = aug, 
                                                      #                                                    solver = private$solver[1]))
                                                      # opt.dist <- NA
                                                      if (mods == "ot_imputer" & !aug) next
                                                      for (match in cur$match[[1]]) {
                                                        for (split in cur$split[[1]]) {
                                                          # if (split) {
                                                          #   mn <- which(cur$match[[1]] == match)
                                                          #   an <- which(cur$model.aug[[1]] == aug)
                                                          #   if (mn > 1 & an > 1) next
                                                          # }
                                                          if (!split) {
                                                            if (match || !aug) next
                                                          }
                                                          data.table::set(private$output.dt, i = iter, j = "estimand" , value = est)
                                                          data.table::set(private$output.dt, i = iter, j = "model" , value = mods)
                                                          data.table::set(private$output.dt, i = iter, j = "model.augmentation" , value = aug)
                                                          data.table::set(private$output.dt, i = iter, j = "match" , value = match)
                                                          data.table::set(private$output.dt, i = iter, j = "split.model" , value = split)
                                                          
                                                          data.table::set(private$output.dt, i = iter, j = "solver" , value = solver)
                                                          data.table::set(private$output.dt, i = iter, j = "delta" , value = if (!is.null(o$delta)){o$delta} else {NA_real_})
                                                          data.table::set(private$output.dt, i = iter, j = "options" , value = list(list(o)))
                                                          # if (!private$grid.search) private$wass.calc(iter, est)
                                                          private$wass.calc(iter, est)
                                                          data.table::set(private$output.dt, i = iter, j = "ess.frac", value = ess.frac)
                                                          #PSIS
                                                          data.table::set(private$output.dt, i = iter, j = "psis.ess.frac", value = psis.ess.frac)
                                                          data.table::set(private$output.dt, i = iter, j = "psis.k", value = psis.k)
                                                          #opt wt dist
                                                          # data.table::set(private$output.dt, i = iter, j = "opt.dist", value = opt.dist)
                                                          
                                                          esteff <- estimate_effect(private$simulator, 
                                                                                    formula = cur$outcome.formula[[1]][[aug + 1]],
                                                                                    weights = private$weights[[est]],
                                                                                    hajek = TRUE,
                                                                                    p = 1,
                                                                                    doubly.robust = aug,
                                                                                    matched = match,
                                                                                    split.model = split,
                                                                                    estimand = est,
                                                                                    model = match.fun(mods))
                                                          if ( method == "Wasserstein" && !isFALSE(private$wass.opt$confidence.interval )) {
                                                            ci.out <- confint(esteff, method = private$wass.opt$confidence.interval)
                                                            data.table::set(private$output.dt, i = iter, j = "confidence.interval", 
                                                                            value = list(c(ci.out$CI)))
                                                          }
                                                          data.table::set(private$output.dt, i = iter, j = "estimate", 
                                                                          value = esteff$estimate
                                                          )
                                                          data.table::set(private$output.dt, i = iter, j = "E_Y1", 
                                                                          value = esteff$variance.components$E_Y1
                                                          )
                                                          data.table::set(private$output.dt, i = iter, j = "E_Y0", 
                                                                          value = esteff$variance.components$E_Y0
                                                          )
                                                          iter <- iter + 1L
                                                        }
                                                        
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
                                                  # data.table::set(private$output.dt, i = iter, j = "opt.dist", value = NA_real_)
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
                                          get_balconst = function(opts, estimand, method) {
                                            if (method == "Constrained Wasserstein" | method == "Wasserstein" ) {
                                                tf.est <- private$temp.output[["SBW"]]$estimand == estimand
                                                tf.delta <- private$temp.output[["SBW"]]$formula == opts$formula
                                                if (all(is.na(tf.delta))) tf.delta <- rep(TRUE, length(tf.delta))
                                                select.rows <- which(tf.delta & tf.est)[1]
                                                temp.delt <- private$temp.output[["SBW"]][select.rows, "delta"][[1]]
                                                # delta <- temp.wass$dist[[idx]]
                                             
                                              return(temp.delt)
                                              
                                            } else {
                                              return(NA)
                                            } 
                                          },
                                          get_delta = function(opts, estimand, method) {
                                            if(method == "Constrained Wasserstein" & isFALSE(opts$add.margins) & !private$grid.search) {
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
                                              if(estimand == "ATE") {
                                                out <- lapply(delta, function(w) list(joint = w))
                                              } else {
                                                out <- list(joint = delta)
                                              }
                                              # if(is.null(delta)) delta <- NA
                                              return(out)
                                            # } else if (method == "Wasserstein" | method == "RKHS") {
                                            # } else if ( (method == "Constrained Wasserstein" |
                                            #            method == "Wasserstein" ) & 
                                            #            isTRUE(opts$add.margins) &
                                            #            !private$grid.search) {
                                            #     if (private$cwass.targ == "SBW") {
                                            #       tf.est <- private$temp.output[["SBW"]]$estimand == estimand
                                            #       tf.delta <- private$temp.output[["SBW"]]$delta == opts$delta
                                            #       if (all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                            #       select.rows <- which(tf.delta & tf.est)[1]
                                            #       temp.wass <- private$temp.output[["SBW"]][select.rows, "wasserstein"][[1]]
                                            #       # idx <- which(temp.wass$metric == opts$metric & 
                                            #       #                temp.wass$ground_p == opts$ground_p &
                                            #       #                temp.wass$wass_p == opts$wass_p )
                                            #       # delta <- temp.wass$dist[[idx]]
                                            #     } else if (private$cwass.targ == "RKHS") {
                                            #       tf.est <- private$temp.output[["RKHS"]]$estimand == estimand
                                            #       # tf.delta <- private$temp.output[["RKHS"]]$delta == opts$delta
                                            #       # if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                            #       select.rows <- which(tf.est)[1]
                                            #       temp.wass <- private$temp.output[["RKHS"]][select.rows, "wasserstein"][[1]]
                                            #       # idx <- which(temp.wass$metric == opts$metric & 
                                            #       #                temp.wass$ground_p == opts$ground_p &
                                            #       #                temp.wass$wass_p == opts$wass_p )
                                            #       # delta <- temp.wass$dist[[idx]]
                                            #     } 
                                            #     # else if (private$cwass.targ == "RKHS.dose") {
                                            #     #   tf.est <- private$temp.output[["RKHS.dose"]]$estimand == estimand
                                            #     #   tf.delta <- private$temp.output[["RKHS.dose"]]$delta == opts$delta
                                            #     #   if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                            #     #   select.rows <- which(tf.delta & tf.est)[1]
                                            #     #   temp.wass <- private$temp.output[["RKHS.dose"]][select.rows, "wasserstein"][[1]]
                                            #     #   idx <- which(temp.wass$metric == opts$metric & 
                                            #     #                  temp.wass$ground_p == opts$ground_p &
                                            #     #                  temp.wass$wass_p == opts$wass_p )
                                            #     # }
                                                # nc <- private$simulator$get_p()
                                                # min_x <- matrix(apply(private$simulator$get_x(),2,min), 1,nc, byrow = TRUE)
                                                # max_x <- matrix(apply(private$simulator$get_x(),2,max), 1,nc, byrow = TRUE)
                                                # delta <- 0.05 * cost_calc_lp(min_x, max_x, opts$wass_p)^opts$wass_p
                                            #     # delta.joint <- temp.wass #temp.wass$dist[[idx]]
                                            #     
                                            #     delta.marg <- if (opts$metric == "Lp" ) {
                                            #       scale <-  matrixStats::colSds(private$simulator$get_x())
                                            #      sapply(delta.joint, function(dd) dd^opts$wass_p * scale^opts$wass_p/sum(scale^opts$wass_p)) 
                                            #     } else {
                                            #       lapply(delta.joint, function(dd) rep((dd^opts$wass_p/nc)^(1/opts$wass_p), nc))
                                            #     }
                                            #     if (private$wass.opt$add.joint) {
                                            #       if (estimand == "ATE") {
                                            #         delta <- list(list(margins = delta.marg[[1]], joint = delta.joint[[1]]),
                                            #              list(margins = delta.marg[[2]], joint = delta.joint[[2]]))
                                            #       } else {
                                            #         delta <- list(margins = delta.marg, joint = delta.joint)
                                            #       }
                                            #     } else {
                                            #       delta <- list(margins = delta.marg)
                                            #     }
                                            #     return(delta)
                                            } else if (method == "Wasserstein" & !private$grid.search) {
                                              # if (private$cwass.targ == "SBW") {
                                              #   tf.est <- private$temp.output[["SBW"]]$estimand == estimand
                                              #   tf.delta <- private$temp.output[["SBW"]]$delta == opts$delta
                                              #   if (all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                              #   select.rows <- which(tf.delta & tf.est)[1]
                                              #   temp.wass <- private$temp.output[["SBW"]][select.rows, "wasserstein"][[1]]
                                              #   # idx <- which(temp.wass$metric == opts$metric & 
                                              #   #                temp.wass$ground_p == opts$ground_p &
                                              #   #                temp.wass$wass_p == opts$wass_p )
                                              #   # delta <- temp.wass$dist[[idx]]
                                              # } else if (private$cwass.targ == "RKHS") {
                                              #   tf.est <- private$temp.output[["RKHS"]]$estimand == estimand
                                              #   # tf.delta <- private$temp.output[["RKHS"]]$delta == opts$delta
                                              #   # if(all(is.na(tf.delta)) & isTRUE(private$grid.search)) tf.delta <- rep(TRUE, length(tf.delta))
                                              #   select.rows <- which(tf.est)[1]
                                              #   temp.wass <- private$temp.output[["RKHS"]][select.rows, "wasserstein"][[1]]
                                              #   # idx <- which(temp.wass$metric == opts$metric & 
                                              #   #                temp.wass$ground_p == opts$ground_p &
                                              #   #                temp.wass$wass_p == opts$wass_p )
                                              #   # delta <- temp.wass$dist[[idx]]
                                              # } 
                                              if(is.na(opts$delta) || is.null(opts$delta)) {
                                                cost_temp <- private$costs[[opts$metric]][[as.character(opts$wass_p)]][estimand]
                                                if(is.null(cost_temp)) {
                                                  nc <- private$simulator$get_p()
                                                  
                                                  if(opts$metric == "sdLp") {
                                                    x  <- scale(private$simulator$get_x())
                                                  } else if (opts$metric == "mahalanobis") {
                                                    U  <- inv_sqrt_mat(private$simulator$get_x())
                                                    x  <- scale(private$simulator$get_x(), scale = FALSE) %*% U
                                                  } else {
                                                    x <- private$simulator$get_x()
                                                  }
                                                  min_x <- matrix(apply(x,2,min), 1, nc, byrow = TRUE)
                                                  max_x <- matrix(apply(x,2,max), 1, nc, byrow = TRUE)
                                                  # delta <- c(0.05 * cost_calc_lp(min_x, max_x, opts$wass_p)^opts$wass_p)
                                                  delta <- rep(c(cost_calc_lp(min_x, max_x, opts$wass_p)^opts$wass_p), 2)
                                                  
                                                } else {
                                                  if(estimand == "ATE") {
                                                    delta <- sapply(cost_temp[[1]], function(cc) max(cc)^opts$wass_p)
                                                  } else {
                                                    delta <- max(cost_temp[[1]])^opts$wass_p
                                                  }
                                                  
                                                }
                                                
                                              } else {
                                                delta <- rep(opts$delta, 2)
                                              }
                                              
                                              if (estimand == "ATE") {
                                                out <- list(list(penalty = delta[1]), list(penalty = delta[2]))
                                              } else {
                                                out <- list(penalty = delta[1])
                                              }
                                              
                                              return(out)
                                            } else if (method == "RKHS") {
                                              return(NA_real_)
                                            } else if (method == "SCM") {
                                              if (is.null(opts$delta) ) {
                                                return(list(penalty = 1))
                                              } else {
                                                list(penalty = opts$delta)
                                              }
                                              
                                            } else {
                                              return(opts$delta)
                                            }
                                          },
                                          get_cost = function(opt, method, estimand) {
                                            if (grepl("Wasserstein", method) | method == "NNM") {
                                              cost <- if (is.null(opt$metric) | is.null(opt$ground_p)) {
                                                NULL
                                                # } else if (!is.null(opt$theta) & !is.null(opt$gamma)) {
                                                #   private$kernel
                                              } else {
                                                if (opt$metric != "RKHS") {
                                                  private$costs[[opt$metric]][[as.character(opt$ground_p)]][[estimand]]
                                                } else {
                                                  private$costs[["RKHS"]][[estimand]]
                                                }
                                              }
                                            } else if (method == "RKHS") {
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
                                                                                                                           None = "calc_weight",
                                                                                                                           Logistic = "calc_weight",
                                                                                                                           NNM = "calc_weight",
                                                                                                                           SBW = "calc_weight",
                                                                                                                           RKHS = "calc_weight",
                                                                                                                           SCM = "calc_weight",
                                                                                                                           'Constrained Wasserstein' = "calc_weight",    
                                                                                                                           Wasserstein = "calc_weight",
                                                                                                                           NA_character_
                                            ))
                                            private$method.lookup$estimand <- lapply(1:nrows, function(i) private$estimand)
                                            private$method.lookup$outcome.model <- lapply(1:nrows, function(i) private$outcome.model)
                                              
                                            #   if(length(private$outcome.model) == 1) {
                                            #   lapply(1:nrows, function(i) private$outcome.model)
                                            # } else {
                                            #   lapply(private$outcome.model, function(i) i)
                                            # }
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
                                            private$method.lookup$split <- lapply(1:nrows, function(i) switch(private$split, 
                                                                                                              "both" = c(FALSE, TRUE),
                                                                                                              "yes" = TRUE,
                                                                                                              "no" = FALSE))
                                            private$method.lookup$solver <- sapply(private$method, function(mm) switch(mm,
                                                                                                                       Logistic = "glm",
                                                                                                                       NNM = NA_character_,
                                                                                                                       SBW = private$solver,
                                                                                                                       SCM = private$solver,
                                                                                                                       RKHS = private$solver,
                                                                                                                       RKHS.dose = private$solver,
                                                                                                                       'Constrained Wasserstein' = private$solver,    
                                                                                                                       Wasserstein = private$solver,
                                                                                                                       None = NA_character_,
                                                                                                                       NA_character_
                                            ))
                                            sdm <- private$standardized.difference.means
                                            wdc <- private$wasserstein.distance.constraints
                                            lambdas <- private$RKHS$lambdas
                                            theta = list(private$RKHS$theta)
                                            gamma = list(private$RKHS$gamma)
                                            rkhs_p = private$RKHS$p
                                            sigma_2 = private$RKHS$sigma_2
                                            kernel = private$RKHS$kernel
                                            grid.search = private$grid.search
                                            if(isTRUE(private$grid.search)) {
                                              sdm <- NA
                                              lambdas <- NA
                                              wdc <- NA
                                            }
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
                                            nnm_list <- list(metric = private$metric,
                                                                        # ground_p = private$ground_powers,
                                                                        wass_p = private$wass_powers,
                                                                        delta = NA
                                            )
                                            Cwass_list <- list(metric = private$metric,
                                                               # ground_p = private$ground_powers,
                                                               wass_p = private$wass_powers,
                                                               std_diff = switch(private$cwass.targ,
                                                                                 "SBW" = sdm, 
                                                                                 "RKHS" = NA,
                                                                                 "RKHS.dose" = lambdas),
                                                               # RKHS.metric = private$RKHS$metric,
                                                               delta = wdc,
                                                               grid.search = private$grid.search,
                                                               formula = private$ps.formula$cwass,
                                                               add.margins = private$wass.opt$add.margins,
                                                               joint.mapping = private$wass.opt$joint.mapping,
                                                               penalty = private$wass.opt$penalty,
                                                               neg.weights = private$wass.opt$neg.weights
                                            )
                                            wass_list <- list( metric = private$metric,
                                                              # ground_p = private$ground_powers,
                                                              wass_p = private$wass_powers,
                                                              std_diff = switch(private$cwass.targ,
                                                                                "SBW" = sdm, 
                                                                                "RKHS" = NA,
                                                                                "RKHS.dose" = lambdas),
                                                              # RKHS.metric = private$RKHS$metric,
                                                              delta = wdc,
                                                              grid.search = private$grid.search,
                                                              formula = private$ps.formula$wass,
                                                              add.margins = private$wass.opt$add.margins,
                                                              joint.mapping = private$wass.opt$joint.mapping,
                                                              penalty = private$wass.opt$penalty,
                                                              neg.weights = private$wass.opt$neg.weights,
                                                              add.divergence = private$wass.opt$add.divergence
                                                              
                                            )
                                            scm_list <- list(  delta = NA,
                                                               grid.search = FALSE,#private$grid.search,
                                                               # add.margins = private$wass.opt$add.margins,
                                                               joint.mapping = TRUE,
                                                               penalty = "none",#private$wass.opt$penalty,
                                                               neg.weights = private$wass.opt$neg.weights
                                                               
                                            )
                                            # if(!any(private$metric == "RKHS")) wass_list$RKHS.metric <- Cwass_list$RKHS.metric <- NULL
                                            wass_list$std_diff <- NA
                                            private$method.lookup$options <- sapply(private$method, function(mm) switch(mm,
                                                                  None = list(delta = NA_real_, other = NA_character_),
                                                                  Logistic = list(delta = private$truncations,
                                                                                  formula = private$ps.formula$logistic),
                                                                  Probit = list(delta = private$truncations,
                                                                                formula = private$ps.formula$logistic),
                                                                  CBPS = list(delta = NA_real_,
                                                                              formula = private$ps.formula$cbps),
                                                                  NNM = nnm_list,
                                                                  SBW = list(grid.search = private$grid.search,
                                                                             delta = sdm,
                                                                             formula = private$ps.formula$sbw),   
                                                                  SCM = scm_list,
                                                                  RKHS = RKHS_list,
                                                                  RKHS.dose = RKHS.dose_list,
                                                                  'Constrained Wasserstein' = Cwass_list,    
                                                                  Wasserstein = wass_list,
                                                                  COT = wass_list,
                                                                  gp = list(NA)), simplify = FALSE)
                                            if("Logistic" %in% private$method) private$method.lookup$estimand[private$method.lookup$method == "Logistic"][[1]] <- private$estimand[private$estimand != "feasible"]
                                            private$max.conditions <- private$max.cond.calc()
                                            for(i in private$method) private$method.lookup$options[[i]] <- private$set.opts(private$method.lookup[i == private$method.lookup$method,])
                                          },
                                          max.cond.calc = function() {
                                            dims <- sapply(private$method, function(mm){
                                              cur <- private$method.lookup[private$method.lookup$method == mm,]
                                              lens <- sapply(cur$options[[1]], length)
                                              noptions <- prod(lens[lens != 0])
                                              cur$options <- NULL
                                              dim <- prod( sapply(unlist(cur, recursive = FALSE), length) ) * noptions
                                              return(dim)
                                            })
                                            return(max(dims, na.rm = TRUE))
                                          },
                                          set.opts = function(cur) { #sets options for run in estimate function
                                            if (cur$method == "Constrained Wasserstein" & isFALSE(cur$options[[1]]$grid.search)) {
                                              cur$options[[1]]$delta <- cur$options[[1]]$std_diff
                                            }
                                            
                                            df.opts <- do.call("expand.grid", c(cur$options,
                                                               stringsAsFactors = FALSE))
                                            
                                            opts <- lapply(1:nrow(df.opts), function(i) df.opts[i,])
                                            
                                            return(opts)
                                          },
                                          update = function() {
                                            private$simulator$gen_data()
                                            private$cost.setup()
                                            # private$kernel.setup()
                                            if (private$verbose) message("  Method: ",appendLF = FALSE)
                                            for (mm in private$method) {
                                              if (private$verbose) {
                                                if (mm == private$method[length(private$method)]) {
                                                  message(mm)
                                                } else {
                                                  message(mm,", ",appendLF = FALSE)
                                                }
                                              }
                                              if (mm %in% private$weight_based_methods) {
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
                                            n <- sum(private$simulator$get_n())
                                            b <- rep(1/n, n)
                                            if (estimand == "ATE" || estimand == "cATE") {
                                              wass.output <- data.frame(w0 = sinkhorn_geom(x = private$simulator$get_x0(), 
                                                                                     y = private$simulator$get_x(),
                                                                                     a = private$weights[[estimand]]$w0,
                                                                                     b = b, scaling = 0.8,
                                                                                     power = 2,
                                                                                     blur = 100, metric = "Lp")$loss,
                                                                  w1 = sinkhorn_geom(x = private$simulator$get_x1(), 
                                                                                          y = private$simulator$get_x(),
                                                                                          a = private$weights[[estimand]]$w1,
                                                                                          b = b, scaling = 0.8,
                                                                                          power = 2,
                                                                                          blur = 10, metric = "Lp")$loss)
                                            } else if (estimand == "ATC") {
                                              wass.output <- data.frame(w0 = sinkhorn_geom(x = private$simulator$get_x0(), 
                                                                                     y = private$simulator$get_x1(),
                                                                                     a = private$weights[[estimand]]$w0,
                                                                                     b = private$weights[[estimand]]$w1,
                                                                                     power = 2, scaling = 0.8,
                                                                                     blur = 10, metric = "Lp")$loss,
                                                                  w1 = NA)
                                            } else if (estimand == "ATT") {
                                              wass.output <- data.frame(w0 = NA, 
                                                                  w1 = sinkhorn_geom(x = private$simulator$get_x0(), 
                                                                                          y = private$simulator$get_x1(),
                                                                                          a = private$weights[[estimand]]$w0,
                                                                                          b = private$weights[[estimand]]$w1,
                                                                                          power = 2, scaling = 0.8,
                                                                                          blur = 10, metric = "Lp")$loss)
                                            } else {
                                              stop("estimand not found in wass.calc")
                                            }
                                            
                                            data.table::set(private$output.dt, i = iter, j = "wasserstein" , value = list(wass.output))
                                          },
                                          weight.calc = function(cur, estimand, 
                                                                 solver, delta,
                                                                 cost = NULL, 
                                                                 p = NULL,
                                                                 grid.search = FALSE,
                                                                 opt.hyperparam = TRUE,
                                                                 opt.method = c("stan", "optim", "bayesian.optimization"),
                                                                 metric = metric,
                                                                 formula,
                                                                 add.margins = FALSE,
                                                                 add.divergence = FALSE,
                                                                 penalty,
                                                                 joint.mapping) {
                                            method <- as.character(cur$method[[1]])
                                            if (grid.search & method == "SBW") delta <- private$standardized.difference.means
                                            if (grid.search & method == "RKHS.dose") delta <- private$RKHS$lambdas
                                            # if (grid.search & method == "Wasserstein") delta <- private$standardized.difference.means
                                            if (grid.search & method == "Constrained Wasserstein") delta <- private$wasserstein.distance.constraints
                                            if (method == "RKHS") {
                                              if (opt.hyperparam & !is.null(private$RKHS.opt[[estimand]]) ) {
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
                                              if(grid.search) lambda <- 1e2
                                              power <- p[[1]]
                                              # kernel <- private$RKHS$kernel
                                            }
                                            if ( isTRUE(joint.mapping) && isTRUE(add.margins) && isTRUE(penalty != "none"))  {
                                              grid.length <- 4
                                            } else {
                                              grid.length <- 8
                                            }
                                            if (estimand != "cATE") {
                                              private$weights[[estimand]] <- 
                                                tryCatch(
                                                  calc_weight(data = private$simulator,  
                                                            constraint = delta,
                                                            formula = formula,
                                                            estimand = estimand, 
                                                            method = method,
                                                            cost = cost, p = power,
                                                            transport.matrix = TRUE,
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
                                                            iter = if (is.null(private$RKHS$iter)) 2000 else private$RKHS$iter,
                                                            maxit = if (is.null(private$RKHS$iter)) 2000 else private$RKHS$iter,
                                                            metric = metric,
                                                            # balance.constraints = (2 * private$SBW.balconst[[estimand]] + 0.001),
                                                            balance.constraints = private$SBW.balconst[[estimand]],
                                                            add.joint = TRUE, #private$wass.opt$add.joint,
                                                            add.margins = isTRUE(add.margins),
                                                            add.divergence = isTRUE(add.divergence),
                                                            unbiased = TRUE,
                                                            penalty = penalty,
                                                            joint.mapping = isTRUE(joint.mapping),
                                                            grid.length = grid.length,
                                                            wass.method = private$wass.opt$method,
                                                            wass.niter = private$wass.opt$niter,
                                                            epsilon = private$wass.opt$epsilon,
                                                            verbose = isTRUE(private$verbose),
                                                            eval.method = private$wass.opt$eval.method,
                                                            n.boot = if (method == "SBW") {1000}else{100},
                                                            K = 10,
                                                            R = private$wass.opt$cross.val.replicates),
                                                         error = function(e) {
                                                           warning("Error in weight method ", method, " with estimand ",estimand, ". ", e$message)
                                                           ns <- private$simulator$get_n()
                                                           e.out <- list(w0 = rep(NA_real_, ns[1]),
                                                                         w1 = rep(NA_real_, ns[2]),
                                                                         gamma = NULL,
                                                                         estimand = estimand,
                                                                         method = method, args = list(constraint = delta))
                                                           class(e.out) <- "causalWeights"
                                                           return(e.out)
                                                         }
                                                )
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
