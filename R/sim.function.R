sim.function <- function(dataGen, nsims = 100L, ground_p = 2, p = 1, 
                          grid.search = TRUE,
                          match = c("both", "yes", "no"),
                          augmentation = c("both", "yes", "no"),
                         standardized.mean.difference = 0.1,
                         truncations = 0.1,
                         distance = c("Lp", "mahalanobis"),
                         solver = c("cplex","gurobi", "mosek"),
                         seed = NULL, ...) {
  
  nsims <- as.integer(nsims)
  grid.search <- isTRUE(grid.search)
  standardized.mean.difference <- as.numeric(standardized.mean.difference)
  truncations <- as.numeric(truncations)
  ground_p <- as.numeric(ground_p)
  p <- as.numeric(p)
  solver <- match.arg(solver)
  match <- match.arg(match)
  augmentation <- match.arg(augmentation)
  dots <- list(...)
  
  stopifnot(nsims > 0)
  stopifnot(inherits(dataGen, "DataSim"))
  stopifnot(all(standardized.mean.difference > 0))
  
  #### dist mat ####
  dist <- match.arg(distance, several.ok = TRUE)
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  #### run simulations ####
  sh <- SimHolder$new(nsim = nsims,
                dataSim = original,
                grid.search = grid.search,
                truncations = truncations,
                standardized.difference.means = standardized.mean.difference,
                outcome.model = dots$outcome.model,
                outcome.formula = list(none = dots$outcome.formula$none,
                                       augmentation = dots$outcome.formula$augmentation),
                model.augmentation = augmentation,
                match = match,
                solver = solver,
                wass_powers = p,
                ground_powers = ground_p,
                metrics = distance)
  sh$run()
  
  #### Get outputs ####
  sh.output <- sh$get.output()
  outcome <- sh$get.outcome(sh.output)
  wasserstein <- sh$get.wass(sh.output)
  # ESS.frac <- sh$get.ESSfrac(sh.output)
  
  #### create output list ####
  output <- list(outcome,
                 # ESS.frac,
                 wasserstein)
  names(output) <- c("outcome",
                     # "ESS/N",
                     "Wasserstein")
  class(output) <- "simOutput"
  
  return(output)
}

# setGeneric("%dorng%", doRNG::`%dorng%`)
# setGeneric("%dopar%", foreach::`%dopar%`)
