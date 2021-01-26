sqrt_mat <- function(X) {
  p <- ncol(X)
  decomp <- eigen(X)
  return(tcrossprod(decomp$vectors %*% diag(sqrt(decomp$values), p, p), decomp$vectors))
}

inv_sqrt_mat <- function(X) {
  p <- ncol(X)
  decomp <- eigen(X)
  return(tcrossprod(decomp$vectors %*% diag(1/sqrt(decomp$values), p, p), decomp$vectors))
}

mahal_transform <- function(X, Y) {
  p <- ncol(X)
  decomp <- eigen(0.5 * cov(X) + 0.5 * cov(Y))
  L_inv <- tcrossprod(decomp$vectors %*% diag(1/sqrt(decomp$values), p, p), decomp$vectors)
  
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
                       dims = c(rows, rows * cols), giveCsparse = FALSE))
}
vec_to_col_constraints <- function(rows, cols) {
  # ones_rows <- Matrix::Matrix(data = 1, nrow = 1, ncol = rows)
  # diag_cols <- Matrix::Diagonal(cols, x = 1)
  # return(Matrix::kronecker(ones_rows, diag_cols))
  
  
  return( Matrix::sparseMatrix(i = rep(1:cols, each = rows),
                                         j = c(sapply(0:(cols - 1), function(i) i * rows + 1:rows)),
                                         x = rep(1,rows * cols),
                                         dims = c(cols, rows * cols), giveCsparse = FALSE)
  
  )
}

# marg_constraint_to_transport_matrix_row <- function(constraints, rows, cols) {
#   
# }

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

#from PropCis package
z2stat <- function (p1x, nx, p1y, ny, dif) 
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