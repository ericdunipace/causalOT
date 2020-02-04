compare_sbw_mine <- function(dataGen, nsims = 100, 
                           standardized.mean.difference = 0.1,
                           solver = c("cplex","gurobi"),
                           parallel = FALSE, 
                           seed = NULL) {
    
    nsims <- as.integer(nsims)
    standardized.mean.difference <- as.numeric(standardized.mean.difference)

    stopifnot(nsims > 0)
    stopifnot(inherits(dataGen, "DataSim"))
    stopifnot(all(standardized.mean.difference > 0))
    
    #### dist mat ####
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
    
    # if(parallel) {
    #   warning("parallel not yet supported")
    # }
    
    if(!is.null(seed)) {
      set.seed(seed)
    }
    
    #### iterate through metrics and such ####
    addl.terms <- list(standardized.mean.difference = paste0("std_diff_",standardized.mean.difference))
    
    names(standardized.mean.difference) <- addl.terms$standardized.mean.difference
    
    term.values  <- list(standardized.mean.difference = standardized.mean.difference)
    
    fill.df <- list(labels = addl.terms,
                    values = term.values)
    
    #### run simulations ####
    # for(nn in 1:nsims) {
    # simulations <- lapply(1:nsims, function(sim) {
    simulations <- foreach::foreach(sim = 1:nsims, .errorhandling = "pass") %dorng% {
      
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
                                  smd = addl.terms$standardized.mean.difference)
      outcome <- generate_holder(outcome = TRUE,
                                 smd = addl.terms$standardized.mean.difference)
      
      #names to iterate through
      options <- get_holder_options()
      estimates <- options$estimates
      wn <- c("SBW", "SBW_j")
      DR <- options$dr
      mch<- options$matched
      estimates <- estimates[estimates != "feasible"]
      
      for(diffs in addl.terms$standardized.mean.difference) {
        for(i in estimates) {
          for (k in DR) {
            for(l in mch) {
              outcome[[diffs]][[i]][[k]][[l]] <- sapply(wn, function(jlkj) NULL, simplify = FALSE)
            }
          }
        }
      }
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
      for(diffs in addl.terms$standardized.mean.difference) {
              std.mean.diff[[1]] <- standardized.mean.difference[diffs]
              for(i in estimates) {
                for(j in wn ){
                  delta[[1]] <- std.mean.diff[[1]]
                  if(j == "SBW") {
                    if(i != "ATE" ) {
                      weights[[i]][[j]] <- 
                        calc_weight(target,  constraint = delta[[1]],
                                    estimate = i, method = j,
                                    transport.matrix = TRUE,
                                    solver = solver)
                    } else {
                      weights[[i]][[j]] <- 
                        convert_ATE(weights[["ATT"]][[j]], 
                                    weights[["ATC"]][[j]])
                    }
                  } else {
                    if(i != "ATE" ) {
                      weights[[i]][[j]] <- sbw_jose(target,  constraint = delta[[1]],
                                                  estimate = i,
                                                  solver = solver)
                    } else {
                      weights[[i]][[j]] <- 
                        convert_ATE(weights[["ATT"]][[j]], 
                                    weights[["ATC"]][[j]])
                    }
                  }
                  pop_frac[[diffs]][[i]][[j]] <- 
                    ESS(weights[[i]][[j]])/c(n0,n1)
                  
                  for (k in DR) {
                    dr <- grepl("DR", k)
                    # for(l in mch[1:(dr+1)]){
                    for(l in mch){
                      matched <- grepl("Matched", l, fixed = TRUE)
                      outcome[[diffs]][[i]][[k]][[l]][[j]] <- 
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
      
      out <- list(outcome = outcome,
                  pop_frac = pop_frac
      )
      return(out) 
    # })
    }
    if(stopcl) parallel::stopCluster(cl)
    
    outcome <- do.call("rbind", lapply(simulations, function(l) convert_holder_sbw(l$outcome, outcome = TRUE, fill.df) ))
    
    # ATT <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATT)))
    # ATC <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATC)))
    # ATE <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$ATE)))
    # feasible <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$outcome$feasible)))
    pop_frac <- do.call("rbind", lapply(simulations, function(l) convert_holder_sbw(l$pop_frac, outcome = FALSE, fill.df)))
    pop_frac$Population <- ifelse(grepl("Treated", rownames(pop_frac)), "Treated","Control")
    
    
    
    output <- list(outcome,
                   pop_frac)
    names(output) <- c("outcome","ESS/N")
    return(output)
  
}