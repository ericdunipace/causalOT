sim.function <- function(dataGen, nsims = 100L, ground_p = 2, p = 1, 
                          methods = c("Logistic", "SBW", "RKHS", "NNM", "Wasserstein","Constrained Wasserstein", "gp"),
                          grid.search = TRUE,
                          RKHS = list(opt = TRUE, opt.method = "stan"),
                          match = c("both", "yes", "no"),
                          augmentation = c("both", "yes", "no"),
                          standardized.mean.difference = 0.1,
                          wasserstein.distance.constraints = NA,
                          truncations = 0.1,
                          distance = dist.metrics(),
                          wass.method = "networkflow",
                          wass.niter = 0,
                          wass.epsilon = 0.05,
                          solver = c("cplex","gurobi", "mosek"),
                          calculate.feasible = FALSE,
                          seed = NULL, ...) 
{
  
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
  stopifnot(all(standardized.mean.difference >= 0))
  stopifnot(all(truncations >= 0))
  stopifnot(all(truncations <= 1))
  
  
  #### dist mat ####
  dist <- match.arg(distance, several.ok = TRUE)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #### run simulations ####
  sh <- SimHolder$new(nsim = nsims,
                dataSim = dataGen,
                methods = methods,
                grid.search = grid.search,
                RKHS = RKHS,
                truncations = truncations,
                standardized.difference.means = standardized.mean.difference,
                estimands = dots$estimands,
                outcome.model = dots$outcome.model,
                outcome.formula = list(none = dots$outcome.formula$none,
                                       augmentation = dots$outcome.formula$augmentation),
                model.augmentation = augmentation,
                match = match,
                calculate.feasible = calculate.feasible,
                propensity.formula = dots$propensity.formula,
                solver = solver,
                Wass = list(wass_powers = p,
                          # ground_powers = ground_p,
                          metrics = dist,
                          method = wass.method,
                          niter = wass.niter,
                          epsilon = wass.epsilon,
                          wasserstein.distance.constraints = wasserstein.distance.constraints,
                          add.joint = dots$add.joint,
                          add.margins = dots$add.margins,
                          joint.mapping = dots$joint.mapping,
                          penalty = dots$penalty,
                          neg.weights = dots$neg.weights,
                          confidence.interval = dots$confidence.interval,
                          add.divergence = dots$add.divergence),
                verbose = isTRUE(dots$verbose))
  sh$run()
  
  #### Get outputs ####
  sh.output <- sh$get.output()
  outcome <- sh$get.outcome(sh.output)
  # if (isFALSE(grid.search)) {
  wasserstein <- sh$get.wass(sh.output)
  # } else {
  #   wasserstein <- NULL
  # }
  ESS.frac <- sh$get.ESS.frac(sh.output)
  psis <- sh$get.psis(sh.output)
  
  #### create output list ####
  output <- list(outcome,
                 ESS.frac,
                 psis,
                 wasserstein)
  names(output) <- c("outcome",
                     "ESS/N",
                     "PSIS",
                     "Wasserstein")
  class(output) <- "simOutput"
  
  return(output)
}

# setGeneric("%dorng%", doRNG::`%dorng%`)
# setGeneric("%dopar%", foreach::`%dopar%`)
