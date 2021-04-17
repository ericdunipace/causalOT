cplex_solver <- function(qp, neg.weights = FALSE, get.dual = FALSE, ...) {
  dots <- list(...)
  neg.wt <- as.numeric(isTRUE(neg.weights)) + 1
  get.dual <- isTRUE(get.dual)
  
  qp <- convert_cones(qp)
  qp <- bc_to_dir_const(qp)
  num_param <- length(c(as.numeric(qp$obj$L)))
  
  if (is.null(dots$control)) {
    control <- list(trace = 0L, round = 1L)
  }
  
  if (!is.null(qp$obj$Q)) {
    if (inherits(qp$obj$Q, "ddiMatrix")) qp$obj$Q <- as(as(qp$obj$Q, "dsparseMatrix"), "dtCMatrix")
    if (!(inherits(qp$obj$Q, "dsCMatrix") | inherits(qp$obj$Q, "dtCMatrix"))) qp$obj$Q <- as(qp$obj$Q, "symmetricMatrix")
    qp$obj$Q <- qp$obj$Q * 2
  }
  lb <- qp$bounds$lb
  ub <- qp$bounds$ub
  # lb <- switch(neg.wt,
  #              0,
  #              -Inf)
  fn.capt <- tempfile(pattern = "cplex_capture", fileext = ".txt")
  invisible(capture.output( res <-
                              Rcplex::Rcplex(cvec = c(as.numeric(qp$obj$L)), Amat = qp$LC$A, 
                                             bvec = qp$LC$vals, Qmat = qp$obj$Q,
                                             lb = lb, ub = ub, control = control,
                                             objsense = "min", sense = qp$LC$dir,
                                             vtype = "C", n = 1) , type = "message", file = fn.capt))
  invisible(capture.output(Rcplex::Rcplex.close(), type = "message", append = TRUE, file = fn.capt))
  
  if (res$status != 1) warning("Algorithm did not converge!!!")
  
  sol <- res$xopt[1:num_param]
  sol <- switch(neg.wt,
                sol * as.numeric(sol > 0),
                sol)
 if (get.dual) {
    dual_vars <- res$extra$lambda
    names(dual_vars) <- names(qp$LC$vals)
  } else {
    dual_vars <- NULL
  }
  return(list(sol = sol, dual = dual_vars))
  return(sol[1:qp$nvar])
}

gurobi_solver <- function(qp, neg.weights = FALSE, get.dual = FALSE, ...) {
  neg.wt <- as.numeric(isTRUE(neg.weights)) + 1
  get.dual <- isTRUE(get.dual)
  
  #convert from mosek format to gurobi
  qp <- convert_cones(qp)
  qp <- bc_to_dir_const(qp)
  
  
  num_param <- length(c(as.numeric(qp$obj$L)))
  model <- list()
  model$Q <- qp$obj$Q
  model$modelsense <- 'min'
  model$obj <- c(as.numeric(qp$obj$L))
  model$A <- qp$LC$A
  model$sense <- rep(NA, length(qp$LC$dir))
  model$sense[qp$LC$dir=="E"] <- '='
  model$sense[qp$LC$dir=="L"] <- '<='
  model$sense[qp$LC$dir=="G"] <- '>='
  model$rhs <- qp$LC$vals
  model$vtype <- rep("C", num_param)
  # model$lb <- switch(neg.wt,
  #                    rep(0, num_param),
  #                    rep(-Inf, num_param))
  # model$ub <- rep(Inf, num_param)
  model$lb <- qp$bounds$lb
  model$ub <- qp$bounds$ub
  params <- list(OutputFlag = 0)
  
  
  # dots <- list(...)
  # model$pstart <- dots$init.sol
  
  res <- gurobi::gurobi(model, params)
  if (res$status != "OPTIMAL") {
    # browser()
    warning("Algorithm did not converge!!!")
  }
  sol <- (res$x)[1:num_param]
  sol <- switch(neg.wt,
                sol * as.numeric(sol > 0),
                sol)
  if (all(sol == 0)) stop("All weights are 0!")
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  if (get.dual) {
    dual_vars <- res$pi
    names(dual_vars) <- names(qp$LC$vals)
  } else {
    dual_vars <- NULL
  }
  
  return(list(sol = sol, dual = dual_vars))
  
  if (dots$save.solution) {
    return(list(result = sol, res = res))
  } else {{
    return(list(result = sol))
  }}
}

mosek_solver <- function(qp, neg.weights = FALSE, get.dual = FALSE, ...) {
  neg.wt <- as.numeric(isTRUE(neg.weights)) + 1
  get.dual <- isTRUE(get.dual)
  
  num_param <- length(c(as.numeric(qp$obj$L)))
  
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
          vals <- rep(1, num_param)
        } else {
          vals <-  qp$obj$Q@x
        }
      } else {
        vals <- rep(mean(qp$obj$Q@x), num_param)
      }
      model$qobj <- list(i = 1:num_param, j =  1:num_param, v = 2 * vals)
    } else {
      # if(!inherits(qp$obj$Q, "dsCMatrix")) {
      #   qp$obj$Q <- Matrix::Matrix(qp$obj$Q, sparse = TRUE)
      # }
      not.pos.def <- check_pos_sdef(qp$obj$Q, symmetric = TRUE)
      if (not.pos.def) {
        
        eigs <- RSpectra::eigs_sym(qp$obj$Q, 
                                   k = 2, 
                                   which = "LA", 
                                   sigma = -100,
                       opts = list(retvec = FALSE,
                                   maxitr = 2000,
                                   tol = 1e-7))$values
        if (sum(eigs < 0) == 1) {
          
          decomp <- eigen(qp$obj$Q, symmetric = TRUE)
          qp$obj$Q <- NULL
          d <- sqrt(abs(decomp$values))
          nd <- length(d)
          qvec <- decomp$vec[,c(nd, 1:(nd-1))]
          d <- d[c(nd, 1:(nd-1))]
          if (!is.null(qp$cones)) {
            qp$cones <- list(F = rbind(qp$cones$F,  d * t(qvec)),
                             g = c(qp$cones$g, rep(0, nd)),
                             cones = cbind(qp$cones$cones, matrix(list("QUAD",  nd, NULL),
                                            nrow = 3, ncol = 1)))
          } else {
            qp$cones <- list(F = d * t(qvec),
                             g = rep(0, nd),
                             cones = matrix(list("QUAD",  nd, NULL),
                                            nrow = 3, ncol = 1))
          }
          qp$obj$L <- qp$obj$L + c(d[1] *qvec[,1])
        } else {
          qp$obj$Q <- Matrix::Matrix(pos_sdef(qp$obj$Q, symmetric = TRUE),
                                     sparse = TRUE)
          if (!inherits(qp$obj$Q,"dgTMatrix")) {
            qp$obj$Q <- as(qp$obj$Q, "dgTMatrix")
          }
          trimat <- Matrix::tril(qp$obj$Q)
          model$qobj <- list(i = trimat@i + 1, j = trimat@j + 1, v = trimat@x * 2)
        }
      } else {
        if (!inherits(qp$obj$Q,"dgTMatrix")) {
          qp$obj$Q <- as(qp$obj$Q, "dgTMatrix")
        }
        trimat <- Matrix::tril(qp$obj$Q)
        model$qobj <- list(i = trimat@i + 1, j = trimat@j + 1, v = trimat@x * 2)
      }
      # M <- qp$obj$Q + Matrix::t(qp$obj$Q)
      # not.pos.def <- tryCatch(isFALSE(is.matrix(chol(diag(10000, 10,10)))), error = function(e) TRUE)
      # if (not.pos.def) {
      # qp$obj$Q<- as(qp$obj$Q, "dsCMatrix")
      
      # }
      
    }
  }
  model$c <- c(as.numeric(qp$obj$L))
  if (is.null(model$c)) model$c <- rep(0, num_param)
  model$A <- qp$LC$A
  model$bx <- rbind(qp$bounds$lb, qp$bounds$ub)
  
  #   switch(neg.wt,
  #   rbind(blx = rep(0, num_param), bux = rep(Inf, num_param)),
  #   rbind(blx = rep(-Inf, num_param), bux = rep(Inf, num_param))
  # )
  
  #constraint bounds
  # num_bounds <- sum(qp$LC$dir == "E") + sum(qp$LC$dir == "L") + sum(qp$LC$dir == "G")
  # if (num_bounds > 0) {
  #   blc <- rep(-Inf, num_bounds)
  #   buc <- rep(Inf, num_bounds)
  #   blc[qp$LC$dir=="E"] <- buc[qp$LC$dir=="E"] <- qp$LC$vals[qp$LC$dir=="E"]
  #   buc[qp$LC$dir=="L"] <- qp$LC$vals[qp$LC$dir=="L"]
  #   blc[qp$LC$dir=="G"] <- qp$LC$vals[qp$LC$dir=="G"]
  #   
  #   model$bc <- rbind(blc = blc,
  #                     buc = buc)
  # }
  
  model$bc <- rbind(blc = qp$LC$lc,
                    buc = qp$LC$uc)
  
  model$cones <- qp$cones$cones
  model$F <- qp$cones$F
  model$g <- qp$cones$g
  
  dots <- list(...)
  # 
  # model$sol <- dots$sol
  # model$iparam <- list(OPTIMIZER = dots$OPTIMIZER)
  # 
  # if(!is.null(dots$sol)) model$iparam$OPTIMIZER <- "MSK_OPTIMIZER_FREE_SIMPLEX"#"FREE_SIMPLEX" # "MSK_OPTIMIZER_FREE_SIMPLEX"
  # 
  opts <- list()
  if (is.null(dots$verbose)) {
    opts$verbose <- 0L
  } else {
    opts$verbose <- as.integer(dots$verbose)
  }
  
  opts$usesol <- dots$usesol
  opts$useparam <- dots$useparam
  opts$soldetail <- dots$soldetail
  opts$getinfo <- dots$getinfo
  opts$writebefore <- dots$writebefore
  opts$writeafter <- dots$writeafter
  
  
  # dots <- list(...)
  # model$sol <- c(dots$init.sol)
  
  # model$dparam <- list(ANA_SOL_INFEAS_TOL = 1e-6)
  
  res <- Rmosek::mosek(problem = model, opts = opts)
  if (is.nan(res$response$code) || res$response$code != 0 || res$sol$itr$solsta != "OPTIMAL") {
    # browser()
    warning("Algorithm did not converge!!! Mosek solver message: ", res$response$msg)
  }
  if (res$sol$itr$solsta == "PRIMAL_INFEASIBLE_CER" || res$sol$itr$prosta == "PRIMAL_INFEASIBLE") stop("Problem infeasible")
  sol <- (res$sol$itr$xx)[1:qp$nvar]
  sol <- switch(neg.wt,
                sol * as.numeric(sol > 0),
                sol)
  
  if (all(sol == 0)) stop("All weights are 0!")
  # obj_total <- out$obj
  # 
  # status <- out$status
  # 
  if (get.dual) {
    dual_vars <- cbind(res$sol$itr$slc, res$sol$itr$suc)
    rownames(dual_vars) <- names(qp$LC$vals)
  } else {
    dual_vars <-  NULL
  }
  
  return(list(sol = sol, dual = dual_vars))
  
  
  if (dots$save.solution) {
    return(list(result = sol, res = res))
  } else {{
    return(list(result = sol))
  }}
}

QPsolver <- function(qp, solver = c("mosek","gurobi","cplex"), ...) {
  solver <- match.arg(solver)
  # neg.weights <- isTRUE(list(...)$neg.weights)
  
  # solve.fun <- switch(solver,
  #                     "cplex" = "cplex_solver",
  #                     "gurobi" = "gurobi_solver",
  #                     "mosek" = "mosek_solver")
  # sol <- do.call(solve.fun, list(qp, ...))
  res <- switch(solver,
                "cplex" = cplex_solver(qp, ...),
                "gurobi" = gurobi_solver(qp, ...),
                "mosek" = mosek_solver(qp, ...))
  
  # res$sol <- renormalize(res$sol)
  return(res)
  
  # sol$result <- renormalize(sol$result)
  # return(sol)
}


convert_cones <- function(qp) {
  if (!is.null(qp$cones)) {
    nvar <- qp$nvar
    cones <- qp$cones$cones
    varnums <- 1:sum(unlist(cones[2,]))
    idx     <- unlist(sapply(1:length(cones[2,]), function(i) rep(i, cones[2,i][[1]])))
    cumcone <- cones
    cumcone[2,] <- lapply(1:ncol(cones), function(i) varnums[idx == i])
    which.cones <- sapply(cumcone[1,], 
                          function(cc) which(cc == "RQUAD"))
    vars <- cumcone[2, which.cones]
    # Fmat <- qp$cones$F
    Qc.list <- lapply(vars, function(i) qp$cones$F[i,])
    constraint.list <- lapply(vars, function(i) qp$cones$g[i])
    
    qc     <- vector("list", length(Qc.list))
    idx.constraint <- rep(NA_integer_, length(Qc.list))
    lambda.list <- vector("list", length(Qc.list))
    for (i in seq_along(Qc.list) ) {
      idx.constraint[[i]] <- which(diff(Qc.list[[i]][1,, drop = FALSE]@p) == 1)
      
      qc[[i]] <- Qc.list[[i]][-c(1:2),-idx.constraint[[i]],drop = FALSE]#,
      lambda.list[[i]] <- rep(qp$obj$L[idx.constraint[[i]] ], nrow(qc[[i]]))
      constraint.list[[i]] <- constraint.list[[i]][-c(1:2)]
    }
    qp$obj$L <- qp$obj$L[-idx.constraint]
    qp$LC$A <- qp$LC$A[ , -idx.constraint, drop = FALSE]
    QQ      <- do.call("rbind", qc)
    QQ      <- cbind(QQ, Matrix::Diagonal(n = nrow(QQ), x = -1))
    addl.varnum <- nrow(QQ)
    qp$LC$A <- rbind(cbind(qp$LC$A, Matrix::sparseMatrix(i = integer(0),
                                                         j = integer(0),
                                                         x = 0,
                                                         dims = c(nrow(qp$LC$A),
                                                                  addl.varnum))),
                     QQ)
    if (is.null(qp$obj$Q)) {
      qp$obj$Q <- Matrix::bdiag(Matrix::Diagonal(n = length(qp$obj$L), x = 0.0) ,
                                Matrix::Diagonal(n = addl.varnum, x = 0.5 * unlist(lambda.list)))
    } else {
      qp$obj$Q <- Matrix::bdiag(qp$obj$Q[-idx.constraint, -idx.constraint], 
                                Matrix::Diagonal(n = nrow(QQ), x = 0.5 * unlist(lambda.list)))
    }
    qp$obj$L <- c(qp$obj$L, rep(0, addl.varnum))
    
    qp$LC$vals <- c(qp$LC$vals, unlist(constraint.list))
    qp$LC$dir <- c(qp$LC$dir, rep("E", addl.varnum))
    
    qp$bounds$lb <- c(qp$bounds$lb[-idx.constraint], rep(0, addl.varnum))
    qp$bounds$ub <- c(qp$bounds$ub[-idx.constraint], rep(Inf, addl.varnum))
  }
  return(qp)
}
