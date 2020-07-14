cplex_solver <- function(qp, ...) {
  dots <- list(...)
  
  num_param <- length(c(qp$obj$L))
  
  if(is.null(dots$control)) {
    control <- list(trace = 0L, round = 1L)
  }
  
  if(!is.null(qp$obj$Q)) {
    if(inherits(qp$obj$Q, "ddiMatrix")) qp$obj$Q <- as(as(qp$obj$Q, "dsparseMatrix"), "dtCMatrix")
    if(!(inherits(qp$obj$Q, "dsCMatrix") | inherits(qp$obj$Q, "dtCMatrix"))) qp$obj$Q <- as(qp$obj$Q, "dsCMatrix")
    qp$obj$Q <-qp$obj$Q * 2
  }
  
  fn.capt <- tempfile(pattern = "cplex_capture", fileext = ".txt")
  invisible(capture.output( res <-
                              Rcplex::Rcplex(cvec = c(qp$obj$L), Amat = qp$LC$A, 
                                             bvec = qp$LC$vals, Qmat = qp$obj$Q,
                                             lb = 0, ub = Inf, control=control,
                                             objsense = "min", sense = qp$LC$dir,
                                             vtype = "C", n = 1) , type = "message", file = fn.capt))
  invisible(capture.output(Rcplex::Rcplex.close(), type = "message", append = TRUE, file = fn.capt))
  
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
  params <- list(OutputFlag = 0)
  
  # dots <- list(...)
  # model$pstart <- dots$init.sol
  
  res <- gurobi::gurobi(model, params)
  if(res$status != "OPTIMAL") {
    # browser()
    warning("Algorithm did not converge!!!")
  }
  sol <- (res$x)[1:num_param]
  sol[sol<0] <- 0
  if(all(sol == 0)) stop("All weights are 0!")
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
  # Quad <- if(!inherits(qp$obj$Q, 'dgTMatrix')) {
  #   as(qp$obj$Q, "dgtTMatrix")
  # } else {
  #   qp$obj$Q
  # }
  # ind <- which(Quad@i > Quad@j)
  # Quad_lower <- Matrix::sparseMatrix(i = Quad@i[ind], 
  #                                j = Quad@j[ind],
  #                                x = Quad@x[ind] ,
  #                                dims = dim(Quad),
  #                                giveCsparse = F, index1 = FALSE)
  # model$qobj <- list(i = Quad_lower@i+1, j = Quad_lower@j+1, v = Quad_lower@x)
  if(!is.null(qp$obj$Q) ){
    if(Matrix::isDiagonal(qp$obj$Q)) {
      if("diag" %in% slotNames(qp$obj$Q)) {
        if(qp$obj$Q@diag == "U") {
          vals <- 1
        } else {
          vals <-  mean(qp$obj$Q@x)
        }
      } else {
        vals <- mean(qp$obj$Q@x)
      }
      model$qobj <- list(i = 1:num_param, j =  1:num_param, v = rep(2*vals,num_param))
    } else {
      if(!inherits(qp$obj$Q,"dgTMatrix")) {
        qp$obj$Q <- as(qp$obj$Q, "dgTMatrix")
      }
      trimat <- Matrix::tril(qp$obj$Q)
      model$qobj <- list(i = trimat@i+1, j = trimat@j+1, v = trimat@x * 2)
    }
  }
  model$c <- c(qp$obj$L)
  model$A <- qp$LC$A
  model$bx <- rbind(blx = rep(0, num_param), bux = rep(Inf, num_param))
  
  #constraint bounds
  num_bounds <- sum(qp$LC$dir == "E") + sum(qp$LC$dir == "L") + sum(qp$LC$dir == "G")
  if(num_bounds > 0) {
    blc <- rep(-Inf, num_bounds)
    buc <- rep(Inf, num_bounds)
    blc[qp$LC$dir=="E"] <- buc[qp$LC$dir=="E"] <- qp$LC$vals[qp$LC$dir=="E"]
    buc[qp$LC$dir=="L"] <- qp$LC$vals[qp$LC$dir=="L"]
    blc[qp$LC$dir=="G"] <- qp$LC$vals[qp$LC$dir=="G"]
    
    model$bc <- rbind(blc = blc,
                      buc = buc)
  }
  
  opts <- list(verbose = 0L)
  
  # dots <- list(...)
  # model$sol <- c(dots$init.sol)
  
  # model$dparam <- list(ANA_SOL_INFEAS_TOL = 1e-6)
  
  res <- Rmosek::mosek(problem = model, opts = opts)
  if(res$response$code != 0) {
    # browser()
    warning("Algorithm did not converge!!! Mosek solver message: ", res$response$msg)
  }
  sol <- (res$sol$itr$xx)[1:num_param]
  sol[sol<0] <- 0
  
  if(all(sol == 0)) stop("All weights are 0!")
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  # dual_vars <- out$pi[-length(out$pi)]
  return(sol)
}

QPsolver <- function(qp, solver = c("cplex","gurobi","mosek"), ...) {
  solver <- match.arg(solver)
  
  # solve.fun <- switch(solver,
  #                     "cplex" = "cplex_solver",
  #                     "gurobi" = "gurobi_solver",
  #                     "mosek" = "mosek_solver")
  # sol <- do.call(solve.fun, list(qp, ...))
  sol <- switch(solver,
                "cplex" = cplex_solver(qp, ...),
                "gurobi" = gurobi_solver(qp, ...),
                "mosek" = mosek_solver(qp, ...))
  return(renormalize(sol))
}
