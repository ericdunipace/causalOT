pos_sdef <- function(X, symmetric = FALSE) {
  p <- ncol(X)
   if (inherits(X, "dsTMatrix")) {
     X <- as(as(X, "dsCMatrix"),"dgCMatrix")
     symmetric <- TRUE
   }
   if (inherits(X, "dgTMatrix")) {
     X <- as(X, "dgCMatrix")
     symmetric <- TRUE
   }
  if (symmetric) {
    # emax <- RSpectra::eigs_sym(X, k = 1, which = "LM")$values
    emin <- RSpectra::eigs_sym(X, k = 1, which = "LA", sigma = -100,
                               opts = list(retvec = FALSE,
                                           maxitr = 2000,
                                           tol = 1e-7))$values
  } else {
    # emax <- RSpectra::eigs(X, k = 1, which = "LM")$values
    emin <- RSpectra::eigs(X, k = 1, which = "LA", sigma = -100,
                           opts = list(retvec = FALSE,
                                       maxitr = 2000,
                                       tol = 1e-7))$values
  }
  if (emin < 0) {
    adjust <- abs(emin) + 1e-6
  } else {
    return(X)
  }
  return(X + Matrix::Diagonal(n = p, x = adjust))
}

check_pos_sdef <- function(X, symmetric = FALSE) {
  p <- ncol(X)
  if (inherits(X, "dsTMatrix")) {
    X <- as(as(X, "dsCMatrix"),"dgCMatrix")
    symmetric <- TRUE
  }
  if (inherits(X, "dgTMatrix")) {
    X <- as(X, "dgCMatrix")
    symmetric <- TRUE
  }
  if (symmetric) {
    # emax <- RSpectra::eigs_sym(X, k = 1, which = "LM")$values
    emin <- RSpectra::eigs_sym(X, k = 1, which = "LA", sigma = -100,
                               opts = list(retvec = FALSE,
                                           maxitr = 2000,
                                           tol = 1e-7))$values
  } else {
    # emax <- RSpectra::eigs(X, k = 1, which = "LM")$values
    emin <- RSpectra::eigs(X, k = 1, which = "LA", sigma = -100,
                           opts = list(retvec = FALSE,
                                       maxitr = 2000,
                                       tol = 1e-7))$values
  }
  return(emin < 0)
}

robust_sqrt_mat <- function(X) {
  X <- pos_sdef(X, symmetric = TRUE)
  return(Matrix::Matrix(chol(X), sparse = TRUE))
}

# round_pi <- function(f,g, cost, lambda, a, b) {
#   n <- length(a)
#   m <- length(b)
#   
#   f_prime <- log(a) - f/lambda + row_log_sum_exp((matrix(g, n, m, byrow = TRUE) - cost)/lambda)
#   f_prime <- f_prime * (f_prime < 0)
#   
#   g_prime <- log(b) - g/lambda + col_log_sum_exp((matrix(f, n, m) - cost)/lambda)
#   g_prime <- g_prime * (g_prime < 0)
#   
#   pi_prime <- exp((matrix(f_prime, n, m) + matrix(g_prime, n, m, byrow = TRUE) - cost)/lambda)
#   err_row <- a - rowSums(pi_prime)
#   err_col <- b - colSums(pi_prime)
#   
#   return(pi_prime + matrix(err_row, n, m) * matrix(err_col,n,m,byrow=TRUE)/ sum(abs(err_row)))
#   
# }

round_pi <- function(raw_pi, a, b) {
  n <- length(a)
  m <- length(b)
  
  x <- a/rowSums(raw_pi)
  x[x > 1] <- 1
  
  y <- b/colSums(raw_pi)
  y[y > 1] <- 1
  
  X <- diag(x)
  Y <- diag(y)
  
  pi_prime <-  matrix(x, n, m) *  raw_pi
  pi_2_prime <- pi_prime * matrix(y, n, m, byrow = TRUE)
  err_row <- a - rowSums(pi_2_prime)
  err_col <- b - colSums(pi_2_prime)
  
  return(pi_2_prime + matrix(err_row, n, m) * matrix(err_col,n,m,byrow=TRUE)/ sum(abs(err_row)))
  
}

sqrt_mat <- function(X, symmetric = FALSE) {
  p <- ncol(X)
  decomp <- eigen(X, symmetric = symmetric)
  return(tcrossprod(decomp$vectors %*% diag(sqrt(abs(decomp$values)), p, p), decomp$vectors))
}

inv_sqrt_mat <- function(X, symmetric = FALSE) {
  p <- ncol(X)
  decomp <- eigen(as.matrix(X), symmetric = symmetric)
  return(tcrossprod(decomp$vectors %*% diag(1/sqrt(abs(decomp$values)), p, p), decomp$vectors))
}

inv_mat <- function(X, symmetric = FALSE) {
  p <- ncol(X)
  decomp <- eigen(as.matrix(X), symmetric = symmetric)
  return(tcrossprod(decomp$vectors %*% diag(1/abs(decomp$values), p, p), decomp$vectors))
}

mahal_transform <- function(X, Y, symmetric = FALSE) {
  p <- ncol(X)
  decomp <- eigen(as.matrix(0.5 * cov(X) + 0.5 * cov(Y)), symmetric = symmetric)
  L_inv <- tcrossprod(decomp$vectors %*% diag(1/sqrt(abs(decomp$values)), p, p), decomp$vectors)
  
  return(list(X %*% L_inv, Y %*% L_inv))
}

vec_to_row_constraints <- function(rows, cols) {
  # #linear algebra formulation, slower
  # ones_cols <- Matrix::Matrix(data = 1, nrow = 1, ncol = cols)
  # diag_rows <- Matrix::Diagonal(rows, x = 1)
  # return(Matrix::kronecker(ones_cols, diag_rows))
  
  col_idx <- seq(1,rows*cols,rows)
  
  return(Matrix::sparseMatrix(i = rep(1:rows, each = cols),
                       j = c(sapply(0:(rows - 1), function(i) i + col_idx)),
                       x = rep.int(1,rows * cols),
                       dims = c(rows, rows * cols), repr = "T"))
}

vec_to_col_constraints <- function(rows, cols) {
  # ones_rows <- Matrix::Matrix(data = 1, nrow = 1, ncol = rows)
  # diag_cols <- Matrix::Diagonal(cols, x = 1)
  # return(Matrix::kronecker(ones_rows, diag_cols))
  
  
  return( Matrix::sparseMatrix(i = rep(1:cols, each = rows),
                                         j = c(sapply(0:(cols - 1), function(i) i * rows + 1:rows)),
                                         x = rep(1,rows * cols),
                                         dims = c(cols, rows * cols), repr = "T")
  
  )
}

vec_to_col_constraints_csparse <- function(rows, cols) {
  # ones_rows <- Matrix::Matrix(data = 1, nrow = 1, ncol = rows)
  # diag_cols <- Matrix::Diagonal(cols, x = 1)
  # return(Matrix::kronecker(ones_rows, diag_cols))
  
  
  return( Matrix::sparseMatrix(i = rep(1:cols, each = rows),
                               j = c(sapply(0:(cols - 1), function(i) i * rows + 1:rows)),
                               x = rep(1,rows * cols),
                               dims = c(cols, rows * cols), repr = "C")
          
  )
}

zero_mat_sp <- function(rows, cols) {
  return(Matrix::sparseMatrix(i = integer(0),
                              j = integer(0),
                              x = 0,
                              dims = c(rows, cols), repr = "T"))
}

form_all_squares <- function(form, data.names) {
  if (is.character(form)) {
    split.form   <- strsplit(form, "~")[[1]]
    form.temp    <- split.form[2]
    form.outcome <- split.form[1]
    
  } else if (inherits(form,"formula")) {
    nform        <- length(form)
    form.temp    <- as.character(form[nform])
    form.outcome <- as.character(form[nform - 1])
    if (form.outcome == "~") form.outcome <- NULL
  }
  form.terms <- strsplit(form.temp, "\\+")[[1]]
  is.square  <- grepl("I\\(\\s*\\.\\^2\\s*\\)", form.terms)
  form.terms <- form.terms[!is.square]
  form.nsq   <- paste0(form.terms, collapse = "+")
  square.terms <- NULL
  if ( any(is.square) ) {
    square.terms <- paste0("I(",data.names, "^2)", collapse = " + ")
  }
  form <- as.formula(paste0(form.outcome,
                            "~",
                            paste0(c(form.nsq, square.terms), 
                                   collapse = " + "))
  )
  return(form)
}

cot.model.matrix <- function(formula, object) {
  model.matrix( formula, data=object, contrasts.arg = 
                  lapply(data.frame(object[,sapply(data.frame(object), is.factor)]),
                         contrasts, contrasts = FALSE))
}

f.call.list <- function(fun, list.args) {
  fun <- as.name(as.character(fun))
  
  list.args <- list.args[!duplicated(names(list.args))]
  list.args <- list.args[!sapply(list.args, is.null)]
  name.args <- lapply(names(list.args), as.name)
  names(name.args) <- names(list.args)
  f.call <- as.call(c(list(fun), name.args))
  
  return(eval(expr = f.call, envir = list.args))
}

f.call.list.no.eval <- function(fun, list.args) {
  fun <- as.name(as.character(fun))
  
  list.args <- list.args[!duplicated(names(list.args))]
  list.args <- list.args[!sapply(list.args, is.null)]
  name.args <- lapply(names(list.args), as.name)
  names(name.args) <- names(list.args)
  f.call <- as.call(c(list(fun), name.args))
  
  return(list(expr = f.call, envir = list.args))
}

#from PropCis package
z2stat <- function(p1x, nx, p1y, ny, dif) 
{
  diff = p1x - p1y - dif
  if (abs(diff) == 0) {
    fmdiff = 0
  }
  else {
    t = ny/nx
    a = 1 + t
    b = -(1 + t + p1x + t * p1y + dif * (t + 2))
    c = dif * dif + dif * (2 * p1x + t + 1) + p1x + t * p1y
    d = -p1x * dif * (1 + dif)
    v = (b/a/3)^3 - b * c/(6 * a * a) + d/a/2
    s = sqrt((b/a/3)^2 - c/a/3)
    if (v > 0) {
      u = s
    }
    else {
      u = -s
    }
    w = (3.141592654 + acos(v/u^3))/3
    p1d = 2 * u * cos(w) - b/a/3
    p2d = p1d - dif
    nxy = nx + ny
    var = (p1d * (1 - p1d)/nx + p2d * (1 - p2d)/ny) * nxy/(nxy - 
                                                             1)
    fmdiff = diff^2/var
  }
  return(fmdiff)
}

#from PropCIs package
scoreci <- function(x, n, conf.level) 
{
  zalpha <- abs(qnorm((1 - conf.level)/2))
  phat <- x/n
  bound <- (zalpha * ((phat * (1 - phat) + (zalpha^2)/(4 * 
                                                         n))/n)^(1/2))/(1 + (zalpha^2)/n)
  midpnt <- (phat + (zalpha^2)/(2 * n))/(1 + (zalpha^2)/n)
  uplim <- round(midpnt + bound, digits = 4)
  lowlim <- round(midpnt - bound, digits = 4)
  cint <- c(lowlim, uplim)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

#from PropCis package
diffpropci <- function(x1, n1, x2, n2, conf.level) 
{
  px = x1/n1
  py = x2/n2
  z = qchisq(conf.level, 1)
  proot = px - py
  dp = 1 - proot
  niter = 1
  while (niter <= 50) {
    dp = 0.5 * dp
    up2 = proot + dp
    score = z2stat(px, n1, py, n2, up2)
    if (score < z) {
      proot = up2
    }
    niter = niter + 1
    if ((dp < 1e-07) || (abs(z - score) < 1e-06)) {
      niter = 51
      ul = up2
    }
  }
  proot = px - py
  dp = 1 + proot
  niter = 1
  while (niter <= 50) {
    dp = 0.5 * dp
    low2 = proot - dp
    score = z2stat(px, n1, py, n2, low2)
    if (score < z) {
      proot = low2
    }
    niter = niter + 1
    if ((dp < 1e-07) || (abs(z - score) < 1e-06)) {
      ll = low2
      niter = 51
    }
  }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

entropy <- function(x) {
  x_pos <- x[x > 0]
  return(sum(-x_pos * log(x_pos)))
}

# log sum exp function
log_sum_exp <- function(x) {
  # if(is.vector(x)) {
  if(all(is.infinite(x))) return(x[1])
  mx <- max(x)
  x_temp <- x - mx
  return(log(sum(exp(x_temp)))+ mx)
  # } else if (is.matrix(x)) {
  #   mx <- apply(x, 1, max)
  #   x_temp <- x - mx
  #   return(log(rowSums(exp(x_temp)))+ mx)
  # }
}

# log sum exp for two vectors
log_sum_exp2 <- function(x,y) {
  mx <- pmax(x,y)
  # if(is.infinite(mx)) return(mx)
  
  temp <- cbind(x,y) - mx
  temp[mx == -Inf,] <- -Inf
  return(log(rowSums(exp(temp))) + mx)
}

# log sum exp function by column
col_log_sum_exp <- function(x) {
  # if(is.vector(x)) {
  if(all(is.infinite(x))) return(x[1])
  mx <- apply(x,2,max)
  mx_mat <- matrix(mx,nrow(x),ncol(x),byrow=TRUE)
  x_temp <- x - mx_mat
  return(log(colSums(exp(x_temp))) + mx)
  # } else if (is.matrix(x)) {
  #   mx <- apply(x, 1, max)
  #   x_temp <- x - mx
  #   return(log(rowSums(exp(x_temp)))+ mx)
  # }
}

# log sum exp function by row
row_log_sum_exp <- function(x) {
  # if(is.vector(x)) {
  if(all(is.infinite(x))) return(x[1])
  mx <- apply(x,1,max)
  mx_mat <- matrix(mx,nrow(x),ncol(x))
  x_temp <- x - mx_mat
  return(log(rowSums(exp(x_temp))) + mx)
  # } else if (is.matrix(x)) {
  #   mx <- apply(x, 1, max)
  #   x_temp <- x - mx
  #   return(log(rowSums(exp(x_temp)))+ mx)
  # }
}

# make vector sum to 1
renormalize <- function(x) {
  if (all(is.na(x))) return(x)
  
  if (isTRUE(any(x < 0))) {
    # warning("Negative weights found! Normalizing to sum to 1 with less accurate function. Make sure negative weights make sense for your problem")
    return(x/sum(x, na.rm = TRUE))
  }
  if (isTRUE(all(x == 0)) ) return(rep(0, length(x)))
  l_x <- log(x)
  return(exp(l_x - log_sum_exp(l_x)))
}

# project weights onto simplex
simplex_proj <- function(y) { #simplex projection of Condat 2015
  N <- length(y)
  v <- v_tilde <- rep(NA_real_, N)
  v_count <- 1
  vt_count <- 0
  v[1] <- y[1]
  rho <- y[1] - 1
  
  for(n in 2:N) {
    if(y[n] > rho){
      rho <- rho + (y[n] - rho)/(v_count + 1)
      if(rho > y[n] - 1) {
        v[v_count + 1] <- y[n]
        v_count <- v_count + 1
      } else {
        v_tilde[(vt_count+1):(v_count + vt_count)] <- v[1:v_count]
        vt_count <- vt_count + v_count
        v[[1]] <- y[n]
        v[2:N] <- NA_real_
        rho <- y[n] - 1
        v_count <- 1
      }
    }
  }
  if(!all(is.na(v_tilde))) { #ie, output non-empty
    v_tilde <- v_tilde[!is.na(v_tilde)]
    for(x in v_tilde) {
      if(x > rho) {
        v[[v_count]] <- x
        v_count <- v_count + 1
        rho <- rho + (x - rho)/v_count
      }
    }
  }
  change <- 1
  v_count <- sum(!is.na(v))
  while(change == 1) {
    change <- 0
    v <- v[!is.na(v)]
    for(n in 1:length(v)) {
      x <- v[n]
      if(x <= rho) {
        v[[n]] <- NA_real_
        v_count <- v_count - 1
        rho <- rho + (rho - x)/v_count
        change <- 1
      }
    }
  }
  tau <- rho
  K <- sum(!is.na(v))
  x <- pmax(y - tau, 0)
  return(x)
}


#' Covert the 2-dimensional index to 1-dimensional index
#'
#' @param i Index of row
#' @param j Index of column
#' @param n Total number of rows
#' @param m Total number of columns
#'
#' @return a 1d index for easy matrix entry
#' 
#' @keywords internal
dist_2d_to_1d <- function (i, j, n, m) {
  valid <- (i >= 1) & (j >= 1) & (i <= n) & (j <= m)
  k <- (j - 1) * n + i
  k[!valid] <- NA_real_
  return(k)
}


