sim.function <- function(dataGen, nsims = 100, ground_p = 2, p = 1, 
                         standardized.mean.difference = 0.1,
                         distance = c("Lp", "mahalanobis"),
                         solver = c("cplex","gurobi"),
                         parallel = FALSE, seed = NULL) {
  
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
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    stopcl <- TRUE
  } else {
    parallel <- FALSE
    stopcl <- FALSE
  }
  
  if(parallel) {
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
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
  for(nn in 1:nsims) {
  # simulations <- foreach::foreach(sim = 1:nsims) %dorng% {

    #### Gen Data ####
    target <- dataGen$clone(deep = TRUE)
    target$gen_data()
    ns <- target$get_n()
    n1 <- ns["n1"]
    n0 <- ns["n0"]
    original_mass <- list(w0=rep(1/n0,n0),
                          w1=rep(1/n1,n1))
    
    #### setup holder data ####
    weights <- generate_holder(outcome = FALSE)
    pop_frac <- generate_holder(outcome = FALSE,
                                smd = addl.terms$standardized.mean.difference,
                                power = addl.terms$wasserstein.power, 
                                ground_p = addl.terms$ground.power,
                                dist = addl.terms$metric)
    wass <- generate_holder(outcome = FALSE,
                            smd = addl.terms$standardized.mean.difference,
                            power = addl.terms$wasserstein.power, 
                            ground_p = addl.terms$ground.power,
                            dist = addl.terms$metric)
    W2 <- generate_holder(outcome = FALSE,
                          smd = addl.terms$standardized.mean.difference,
                          power = addl.terms$wasserstein.power, 
                          ground_p = addl.terms$ground.power,
                          dist = addl.terms$metric)
    outcome <- generate_holder(outcome = TRUE,
                               smd = addl.terms$standardized.mean.difference,
                               power = addl.terms$wasserstein.power, 
                               ground_p = addl.terms$ground.power,
                               dist = addl.terms$metric)
    
    #names to iterate through
    options <- get_holder_options()
    estimates <- options$estimates
    wn <- options$weights
    DR <- options$dr
    mch<- options$matched
    # powers <- addl.terms$power
    # mean.diffs <- as.character(standardized.mean.difference)
    
    # use lists to avoid copy
    cost.calc <- delta <- std.mean.diff <- pp <- gp <- cost <- list(NULL)
    
    #### Naive Outcome
    outcome$Naive <- outcome_model(data = target, weights = original_mass,
                                   doubly.robust = FALSE, hajek = FALSE,
                                   matched = FALSE,
                                   target = "ATE")
    
    #### fill lists ####
    for(dist.name in  addl.terms$metric) {
      cost.calc[[1]] <- switch(dist[dist.name], "Lp" = causalOT::cost_calc_lp,
                          "mahalanobis" = causalOT::cost_mahalanobis)
      for (gpowers in addl.terms$ground.power) {
        gp[[1]] <- ground_p[[gpowers]]
        cost[[1]] <- cost.calc[[1]](X = target$get_x0(), Y = target$get_x1(), 
                               ground_p = gp[[1]], direction = "rowwise")
        for(powers in addl.terms$wasserstein.power) {
          pp[[1]] <- p[powers]
          for(diffs in addl.terms$standardized.mean.difference) {
            std.mean.diff[[1]] <- standardized.mean.difference[diffs]
            # original wass #
            wass[[dist.name]][[gpowers]][[powers]][[diffs]]$original <- 
              wasserstein_p(original_mass$w0, original_mass$w1, p = pp[[1]],
                            cost = cost[[1]])
            W2[[dist.name]][[gpowers]][[powers]][[diffs]]$`pre-match`<- 
              wasserstein_p(original_mass$w0, 
                             original_mass$w1, 
                             p = 2,
                             cost = cost[[1]])
            for(i in estimates) {
              for(j in wn ){
                delta[[1]] <- if(j == "SBW" | j == "Logistic") {
                  std.mean.diff[[1]]
                } else {
                  wass[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][["SBW"]]
                }
                if(i != "ATE" | j == "Logistic"){
                  if(j == "Logistic" & i == "feasible") next
                  weights[[i]][[j]] <- 
                    calc_weight(target,  constraint = delta[[1]],
                                                   estimate = i, method = j,
                                                   cost = cost[[1]], p = pp[[1]],
                                                   transport.matrix = TRUE,
                                                   solver = solver)
                } else {
                  weights[[i]][[j]] <- 
                    convert_ATE(weights[["ATT"]][[j]], 
                                 weights[["ATC"]][[j]],
                                 transport.matrix = TRUE,
                                 cost = cost[[1]], p = pp[[1]])
                }
                pop_frac[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][[j]] <- 
                  ESS(weights[[i]][[j]])/c(n0,n1)
                
                wass[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][[j]] <- 
                  wasserstein_p(a = weights[[i]][[j]], b = NULL,
                                                p = pp[[1]], tplan = NULL, cost = cost[[1]])
                W2[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][[j]] <- 
                if( pp[[1]] == 2) {
                  wass[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][[j]]
                } else {
                  wasserstein_p(a = weights[[i]][[j]], b = NULL,
                                p = 2, tplan = NULL, cost = cost[[1]])
                }
                for (k in DR) {
                  dr <- grepl("DR", k)
                  # for(l in mch[1:(dr+1)]){
                  for(l in mch){
                    matched <- grepl("Matched", l, fixed = TRUE)
                    outcome[[dist.name]][[gpowers]][[powers]][[diffs]][[i]][[k]][[l]][[j]] <- 
                      outcome_model(data = target, weights = weights[[i]][[j]],
                                                                 hajek = TRUE,
                                                                 doubly.robust = dr,
                                                                 matched = matched,
                                                                 target = i)
                  }
                  
                }
              }
            }
          }
        }
      }
    }
    
    out <- list(outcome = outcome,
                W2 = W2,
                wass = wass,
                pop_frac = pop_frac
    )
    return(out) 
  }
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
                     paste(c(p,"Wasserstein"), collapse="-"),
                     "RNG")
  return(output)
}

setGeneric("%dorng%", doRNG::`%dorng%`)
setGeneric("%dopar%", foreach::`%dopar%`)
