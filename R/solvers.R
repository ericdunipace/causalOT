cplex_solver <- function(qp, ...) {
  dots <- list(...)
  
  num_param <- length(c(qp$obj$L))
  
  if(is.null(dots$control)) {
    control <- list(trace = 0L, round = 1L)
  }
  res <- Rcplex::Rcplex(cvec = c(qp$obj$L), Amat = qp$LC$A, 
                        bvec = qp$LC$vals, Qmat = qp$obj$Q,
                        lb = 0, ub = Inf, control=control,
                        objsense = "min", sense = qp$LC$dir,
                        vtype = "C", n = 1)
  Rcplex::Rcplex.close()
  
  if(res$status != 1) warning("Algorithm did not converge!!!")
  
  return(res$xopt[1:num_param])
}

gurobi_solver <- function(qp, ...) {
  num_param <- length(c(qp$ob$L))
  model <- list()
  model$Q <- qp$obj$Q
  model$modelsense <- 'min'
  model$obj <- c(qp$obj$L)
  model$A <- qp$LC$A
  model$sense <- rep(NA, length(qp$LC$dir))
  model$sense[qp$LC$dir=="E"] <- '='
  model$sense[qp$LC$dir=="L"] <- '<='
  model$sense[qp$LC$dir=="G"] <- '>='
  model$rhs <- qp$LC$vals
  model$vtype <- rep("C", num_param)
  model$lb <- rep(0, num_param)
  model$ub <- rep(Inf, num_param)
  params <- list(OutputFlag = 1)
  
  # dots <- list(...)
  # model$pstart <- dots$init.sol
  
  res <- gurobi::gurobi(model, params)
  if(res$status != "OPTIMAL") {
    # browser()
    warning("Algorithm did not converge!!!")
  }
  sol <- (res$x)[1:num_param]
  sol[sol<0] <- 0
  
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  # dual_vars <- out$pi[-length(out$pi)]
  return(sol)
}

mosek_solver <- function(qp, ...) {
  num_param <- length(c(qp$obj$L))
  model <- list()
  model$sense <- 'min'
  model$qobj <- qp$obj$Q
  model$c <- c(qp$obj$L)
  model$A <- qp$LC$A
  model$bx <- rbind(blx = rep(0, num_param), bux = rep(Inf, num_param))
  
  #constraint bounds
  blc <- rep(-Inf, num_param)
  buc <- rep(Inf, num_param)
  blc[qp$LC$dir=="E"] <- buc[qp$LC$dir=="E"] <- qp$LC$vals[qp$LC$dir=="E"]
  buc[qp$LC$dir=="L"] <- qp$LC$vals[qp$LC$dir=="L"]
  blc[qp$LC$dir=="G"] <- qp$LC$vals[qp$LC$dir=="G"]
  
  model$bc <- rbind(blc = blc,
                    buc = buc)
  
  opts <- list(verbose = 0)
  
  # dots <- list(...)
  # model$sol <- c(dots$init.sol)
  
  res <- Rmosek(problem = model, opts = opts)
  if(res$response$code != "OPTIMAL") {
    # browser()
    warning("Algorithm did not converge!!!")
  }
  sol <- (res$sol)[1:num_param]
  sol[sol<0] <- 0
  
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  # dual_vars <- out$pi[-length(out$pi)]
  return(sol)
}

QPsolver <- function(qp, solver = c("cplex","gurobi"), ...) {
  solver <- match.arg(solver)
  
  solve.fun <- switch(solver,
                      "cplex" = "cplex_solver",
                      "gurobi" = "gurobi_solver",
                      "mosek" = "mosek_solver")
  sol <- do.call(solve.fun, list(qp, ...))
  return(renormalize(sol))
}
