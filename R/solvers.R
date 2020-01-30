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
  model$Q <- qp$Q
  model$modelsense <- 'min'
  model$obj <- c(qp$obj$L)
  model$A <- qp$LC$A
  model$sense <- rep(NA, length(qp$LC$dir))
  model$sense[qp$LC$dir=="E"] <- '='
  model$sense[qp$LC$dir=="L"] <- '<='
  model$sense[qp$LC$dir=="G"] <- '>='
  model$rhs <- qp$LC$vals
  model$vtype <- rep("C", num_param)
  model$ub <- rep(Inf, num_param)
  params <- list(OutputFlag = 0)
  
  res <- gurobi::gurobi(model, params)
  
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  # dual_vars <- out$pi[-length(out$pi)]
  return((res$x)[1:num_param])
}

QPsolver <- function(qp, solver = c("cplex","gurobi"), ...) {
  solver <- match.arg(solver)
  
  solve.fun <- switch(solver,
                      "cplex" = cplex_solver,
                      "gurobi" = gurobi_solver)
  sol <- solve.fun(qp, ...)
  return(renormalize(sol))
}
