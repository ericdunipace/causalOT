setClass("causalWeights", slots = c(w0 = "numeric", w1="numeric", gamma = "NULL",estimate = "character"))
calc_weight.DataSim <- function(data, constraint,  estimate = c("ATT", "ATC","feasible"), 
                                method = c("SBW","Wasserstein", "Constrained Wasserstein"),
                                ...) {
  method <- match.arg(method)
  estimate <- match.arg(estimate)
  op <- quadprog(data, constraint,  estimate, 
           method,
           ...)
  dots <- list(...)
  if(is.null(dots$control)) {
    control <- list(trace = 0L, round = 1L)
  }
  # skip_cplex <- FALSE
  # if (method == "Constrained Wasserstein") {
  #   check <- check_wass_const(op)
  #   skip_cplex <- check$skip_cplex
  # }
  # if (skip_cplex) {
  #   res <- list(xopt = check$res, status = 1)
  # } else {
    res <- Rcplex::Rcplex(cvec = c(op$obj$L), Amat = op$LC$A, 
                        bvec = op$LC$vals, Qmat = op$obj$Q,
                        lb = 0, ub = Inf, control=control,
                        objsense = "min", sense = op$LC$dir,
                        vtype = "C", n = 1)
    Rcplex::Rcplex.close()
  # }
  if(res$status != 1) warning("Algorithm did not converge!!!")
  sol <- renormalize(res$xopt) # normalize to have closer to sum 1
  output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
  gamma <- NULL
  ns <- data$get_n()
  if(method == "Wasserstein" | method == "Constrained Wasserstein") {
    gamma <- matrix(sol, ns["n0"], ns["n1"])
    if(estimate == "ATT") {
      sol <- rowSums(gamma)
    } else if (estimate == "ATC") {
      sol <- colSums(gamma)
    } else if (estimate == "feasible") {
      sol <- gamma
    }
    output$gamma <- gamma
  }
  if(estimate == "ATT") {
    output$w0 <- sol
    output$w1 <- rep.int(1/ns["n1"],ns["n1"])
  } else if (estimate == "ATC") {
    output$w0 <- rep.int(1/ns["n0"],ns["n0"])
    output$w1 <- sol
  } else if (estimate == "feasible") {
    if(method == "Wasserstein" | method == "Constrained Wasserstein"){
      output$w0 <- rowSums(gamma)
      output$w1 <- colSums(gamma)
    } else if (method == "SBW") {
      output$w0 <- renormalize(res$xopt[1:ns["n0"]])
      output$w1 <- renormalize(res$xopt[ns["n0"] + 1:ns["n1"]])
    }
  }
  output$estimate <- estimate
  class(output) <- "causalWeights"
  return(output)
}

convert_ATE <- function(weight1, weight2) {
  list_weight <- list(weight1, weight2)
  check.vals <- sapply(list_weight, function(w) w$estimate)
  ATT.pres <- check.vals %in% "ATT"
  ATC.pres <- check.vals %in% "ATC"
  both <- (sum(c(ATT.pres, ATC.pres)) == 2)
  if(!both) stop("One set of weights must be from an ATT estimate and one must be from an ATC to combine")
  
  ATC.weight <- list_weight[[which(ATC.pres)]]
  ATT.weight <- list_weight[[which(ATT.pres)]]
  
  output <- list(w0 = NULL, w1 = NULL, gamma= NULL, estimate = NULL)
  class(output) <- "causalWeights"
  
  output$w0 <- ATT.weight$w0
  output$w1 <- ATC.weight$w1
  output$estimate <- "ATE"
  return(output)
}

setOldClass("DataSim")
setGeneric("calc_weight", function(data, ...) UseMethod("calc_weight"))
setMethod("calc_weight", "DataSim", calc_weight.DataSim)