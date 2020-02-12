sim.function2 <- function(dataGen, nsims = 100, ground_p = 2, p = 1, 
                         standardized.mean.difference = 0.1,
                         distance = c("Lp", "mahalanobis"),
                         solver = c("cplex","gurobi"),
                         parallel = FALSE, 
                         seed = NULL,
                         run.feasible = TRUE) {
  
  nsims <- as.integer(nsims)
  standardized.mean.difference <- as.numeric(standardized.mean.difference)
  ground_p <- as.numeric(ground_p)
  p <- as.numeric(p)
  
  stopifnot(nsims > 0)
  stopifnot(inherits(dataGen, "DataSim"))
  stopifnot(all(standardized.mean.difference > 0))
  
  #### dist mat ####
  dist <- match.arg(distance, several.ok = TRUE)
  
  #### set up cluster if applicable ####
  if(inherits(parallel,"cluster")) {
    cl <- parallel
    parallel <- TRUE
    stopcl <- FALSE
  } else if (isTRUE(parallel)) {
    # cl <- parallel::makeCluster(parallel::detectCores()-1)
    # stopcl <- TRUE
  } else {
    parallel <- FALSE
    stopcl <- FALSE
  }
  
  # if(parallel) {
  #   doParallel::registerDoParallel(cl)
  # } else {
  #   foreach::registerDoSEQ()
  # }
  
  if(parallel) {
    warning("parallel not yet supported")
  }
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  #### iterate through metrics and such ####
  addl.terms <- list(standardized.mean.difference = paste0("std_diff_",standardized.mean.difference),
                     wasserstein.power = paste0("wpower_",p), 
                     ground.power = paste0("gpower_",ground_p),
                     metric = paste0("metric_",dist)
  )
  names(standardized.mean.difference) <- addl.terms$standardized.mean.difference
  names(p)        <- addl.terms$wasserstein.power
  names(ground_p) <- addl.terms$ground.power
  names(dist)   <- addl.terms$metric
  
  term.values  <- list(standardized.mean.difference = standardized.mean.difference,
                       wasserstein.power = p, 
                       ground.power = ground_p,
                       metric = dist)
  
  fill.df <- list(labels = addl.terms,
                  values = term.values)
  
  #### run simulations ####
  SimHolder(nsim = 100,
            dataSim,
            grid.search = TRUE,
            truncations = NULL,
            standardized.difference.means = NULL,
            outcome.model = list("lm"),
            outcome.fomula = list(none = NULL,
                                  augmentation = NULL),
            model.augmentation = "both",
            match = TRUE,
            solver = "gurobi",
            wass_powers = 2,
            ground_powers = 2,
            metric = "Lp")
  if(stopcl) parallel::stopCluster(cl)
  
  outcome <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$outcome, outcome = TRUE, fill.df) ))
  
  # ATT <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATT)))
  # ATC <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATC)))
  # ATE <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATE)))
  # feasible <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$feasible)))
  pop_frac <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$pop_frac, outcome = FALSE, fill.df)))
  pop_frac$Population <- ifelse(grepl("Treated", rownames(pop_frac)), "Treated","Control")
  
  W2 <-  do.call("rbind", lapply(simulations, function(l) convert_holder(l$W2, outcome = FALSE, fill.df)))
  WassP <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$wass, outcome = FALSE, fill.df)))
  
  fac <- as.character(W2$estimate)
  fac[is.na(fac)] <- "Original"
  W2$estimate <- factor(fac, levels = c("Original","ATT","ATC","ATE"))
  
  fac <- as.character(WassP$estimate)
  fac[is.na(fac)] <- "Original"
  WassP$estimate <- factor(fac, levels = c("Original","ATT","ATC","ATE"))
  
  rng <- attr(simulations, "rng")
  
  
  output <- list(outcome,
                 pop_frac,
                 W2,
                 WassP,
                 rng)
  names(output) <- c("outcome","ESS/N",
                     "2-Wasserstein", 
                     "Wasserstein",
                     "RNG")
  return(output)
}

setGeneric("%dorng%", doRNG::`%dorng%`)
setGeneric("%dopar%", foreach::`%dopar%`)
