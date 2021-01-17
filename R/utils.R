sqrt_mat <- function(X) {
  p <- ncol(X)
  decomp <- eigen(X)
  return(tcrossprod(decomp$vectors %*% diag(sqrt(decomp$values), p, p), decomp$vectors))
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

