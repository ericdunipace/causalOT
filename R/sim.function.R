sim.function <- function(dataGen, nsims = 100, ground_p = 2, p = 1, 
                         std_mean_diff = 0.1,
                         distance = c("Lp", "mahalanobis"),
                         parallel = TRUE) {
  
  nsims <- as.integer(nsims)
  stopifnot(nsims > 0)
  stopifnot(inherits(dataGen, "DataSim"))
  
  #### dist mat ####
  dist <- match.arg(distance)
  cost.calc <- switch(dist, "Lp" = causalOT::cost_calc_lp,
                      "mahalanobis" = causalOT::cost_mahalanobis)
  
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
  
  #### run simulations ####
  simulations <- foreach::foreach(sim = 1:nsims) %dorng% {

    #### Gen Data ####
    target <- dataGen$clone()
    target$gen_data()
    ns <- target$get_n()
    n1 <- ns["n1"]
    n0 <- ns["n0"]
    original_mass <- list(w0=rep(1/n0,n0),
                          w1=rep(1/n1,n1))
    cost <- cost.calc(X = target$get_x0(), Y = target$get_x1(), 
                      ground_p = ground_p, direction = "rowwise")
    
    #### setup holder data ####
    weights <- generate_holder(outcome = FALSE)
    pop_frac <- generate_holder(outcome = FALSE)
    wass <- generate_holder(outcome = FALSE)
    W2 <- generate_holder(outcome = FALSE)
    outcome <- generate_holder(outcome = TRUE)
    
    #names to iterate through
    options <- get_holder_options()
    estimates <- options$estimates[1:3]
    wn <- options$weights
    DR <- options$dr
    
    delta <- list(NULL)

    #### original wass ####
    wass$original <- wasserstein_p(original_mass$w0, original_mass$w1, p = p,
                                   cost = cost)
    W2$`pre-match`<- wasserstein_p(original_mass$w0, original_mass$w1, p = 2,
                                   cost = cost)
    
    #### Naive Outcome
    outcome$Naive <- outcome_model(data = target, weights = original_mass,
                                   doubly.robust = FALSE, hajek = FALSE,
                                   target = "ATE")
    
    #### fill lists ####
    for(i in estimates) {
      for(j in wn ){
        delta <- if(j == "SBW") {
          std_mean_diff
        } else {
          wass[[i]][["SBW"]]
        }
        if(i != "ATE"){
          weights[[i]][[j]] <- calc_weight(target,  constraint = delta,
                                           estimate = i, method = j,
                                           cost = cost, p = p)
        } else {
          weights[[i]][[j]] <- convert_ATE(weights[["ATT"]][[j]], weights[["ATC"]][[j]])
        }
        pop_frac[[i]][[j]] <- ESS(weights[[i]][[j]])/c(n0,n1)
        
        wass[[i]][[j]] <- wasserstein_p(a = weights[[i]][[j]], b = NULL,
                                        p = p, tplan = NULL, cost = cost)
        W2[[i]][[j]] <- if( p == 2) {
          wass[[i]][[j]]
        } else {
          wasserstein_p(a = weights[[i]][[j]], b = NULL,
                        p = 2, tplan = NULL, cost = cost)
        }
        for (k in DR) {
          dr <- grepl("DR", k)
          outcome[[i]][[k]][[j]] <- outcome_model(data = target, weights = weights[[i]][[j]],
                                                  doubly.robust = dr, hajek = TRUE,
                                                  target = i)
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
  
  outcome <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$outcome, outcome = TRUE) ))
  
  # ATT <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATT)))
  # ATC <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATC)))
  # ATE <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATE)))
  # feasible <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$feasible)))
  pop_frac <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$pop_frac, outcome = FALSE)))
  pop_frac$Population <- ifelse(grepl("Treated", rownames(pop_frac)), "Treated","Control")

  W2 <-  do.call("rbind", lapply(simulations, function(l) convert_holder(l$W2, outcome = FALSE)))
  WassP <- do.call("rbind", lapply(simulations, function(l) convert_holder(l$wass, outcome = FALSE)))
  
  fac <- as.character(W2$estimate)
  fac[is.na(fac)] <- "Original"
  W2$estimate <- factor(fac, levels = c("Original","ATT","ATC","ATE"))
  
  fac <- as.character(WassP$estimate)
  fac[is.na(fac)] <- "Original"
  WassP$estimate <- factor(fac, levels = c("Original","ATT","ATC","ATE"))
  
  
  output <- list(outcome,
                 pop_frac,
                 W2,
                 WassP)
  names(output) <- c("outcome","ESS/N","2-Wasserstein", paste(c(p,"Wasserstein"),collapse="-"))
  return(output)
}

setGeneric("%dorng%", doRNG::`%dorng%`)
