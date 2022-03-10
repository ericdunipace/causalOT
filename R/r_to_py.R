check_python_modules <- function(method = "auto", conda = "auto") {
  
  if(!reticulate::py_module_available("scipy")) {
    stop("Required python package not installed. Package missing: ","scipy")
    # reticulate::py_install("scipy", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("numpy")) {
    stop("Required python package not installed. Package missing: ","numpy")
    # reticulate::py_install("numpy", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("torch")) {
    stop("Required python package not installed. Package missing: ","torch")
    # reticulate::py_install("torch", method = method, conda = conda)
  }
  
  if(!reticulate::py_module_available("geomloss")) {
    stop("Required python package not installed. Package missing: ","geomloss")
    # reticulate::py_install("geomloss", method = method, conda = conda)
  }
  # if(!reticulate::py_module_available("utils")) {
  #   stop("Required python package not installed. Package missing: ","utils")
  #   # reticulate::py_install("utils", method = method, conda = conda)
  # }
  # if(!reticulate::py_module_available("logging")) {
  #   stop("Required python package not installed. Package missing: ","logging")
  #   # reticulate::py_install("logging", method = method, conda = conda)
  # }
}

skip_if_no_geomloss <- function() {
  have_geomloss <- reticulate::py_module_available("geomloss")
  if (!have_geomloss)
    testthat::skip("geomloss not available for testing")
}



#' Sinkhorn Loss
#' 
#' This function serves as an R wrapper to the Python function SamplesLoss in the GeomLoss package <http://www.kernel-operations.io/geomloss/api/pytorch-api.html?highlight=samplesloss#geomloss.SamplesLoss>
#'
#' @param x covariates for the first set of samples. Should be of class matrix.
#' @param y covariates for the second set of samples. Should be of class matrix.
#' @param a The empirical measure of the first set of samples.
#' @param b The empirical measure of the second set of samples.
#' @param power power of the optimal transport distance.
#' @param blur The finest level of detail that should be handled by the loss function to prevent overfitting on the samples/ locations.
#' @param reach specifies the typical scale associated to the constraint strength
#' @param diameter A rough indication of the maximum distance between points, which is used to tune the espilon-scaling descent and provide a default heuristic for clustering multiscale schemes. If None, a conservative estimate will be computed on-the-fly.
#' @param scaling specifies the ratio between successive values of sigma in the epsilon-scaling descent. This parameter allows you to specify the trade-off between speed (scaling < .4) and accuracy (scaling > .9).
#' @param truncate If backend is "multiscale", specifies the effective support of a Gaussian/Laplacian kernel as a multiple of its standard deviation
#' @param metric Set the metric. One of "Lp","sdLp", or "mahalanobis".
#' @param cluster_scale If backend is "multiscale", specifies the coarse scale at which cluster centroids will be computed. If NULL, a conservative estimate will be computed from diameter and the ambient space's dimension, making sure that memory overflows won't take place.
#' @param debias specifies if we should compute the unbiased Sinkhorn divergence instead of the classic, entropy-regularized "SoftAssign" loss.
#' @param verbose if backend is "multiscale", specifies whether information on the clustering and epsilon-scaling descent should be displayed in the standard output.
#' @param backend one of "auto", "tensorized", "online", or "multiscale"
#' @param ... not currently used. Used to absorb extra arguments passed 
#' by other functions without throwing an error.
#'
#' @return a list with slots "loss", "f", "g". "loss" is the Sinkhorn distance,
#' "f" is the potential corresponding to data `x`, and "g" is the potential
#' corresponding to data `y`.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # requires Python and GeomLoss package
#' x <- stats::rnorm(100, 100, 10)
#' a <- rep(1/100, 100)
#' 
#' y <- stats::rnorm(50, 50, 10)
#' b <- rep(1/50, 50)
#' 
#' sink <- sinkhorn(x = x, y = y, a = a, b = b, power = 2,
#'               metric = "Lp", debias = TRUE)
#' 
#' # sinkhorn distance, de-biased
#' print(sink$loss)
#' 
#' # potentials for first 5 obs in each group
#' print(sink$f[1:5])
#' print(sink$g[1:5])
#' }
sinkhorn <- function(x, y, a, b, power = 2, 
                          blur = 0.05, reach = NULL, diameter = NULL,
                          scaling = 0.5, truncate = 5,
                          metric = "Lp",
                          cluster_scale=NULL, 
                          debias=TRUE, 
                          verbose=FALSE, backend='auto', ... ) {
  
  # retrieve python packages
  np <- reticulate::import("numpy", convert = TRUE)
  torch <- reticulate::import("torch", convert = TRUE)
  geomloss <- reticulate::import("geomloss", convert = TRUE)
  # cmake <- reticulate::import("cmake", convert = TRUE)
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)
  
  # get data dimensions
  n <- nrow(x)
  m <- nrow(y)
  d <- ncol(x)
  
  # adjust data if metric is "sdLp" or "mahalanobis"
  if (metric == "mahalanobis") {
    total <- rbind(x,y)
    U <- inv_sqrt_mat(cov(total), symmetric = TRUE)
    
    update <- (total - matrix(colMeans(total), nrow = n+m,
                                       ncol = d, byrow = TRUE)) %*% U
    
    x <- update[1:n,,drop = FALSE]
    y <- update[-(1:n),,drop = FALSE]
    
  } else if (metric == "sdLp") {
    total <- rbind(x,y)
    update <- scale(total)
    x <- update[1:n,,drop = FALSE]
    y <- update[-(1:n),,drop = FALSE]
  }
  
  # set up appropriate cost function for backend
  if (backend == "tensorized" || (n*m <= 5000^2 && backend != "multiscale" && backend != "online")) {
    if (power == 2) {
      cost <- geomloss$utils$squared_distances
    } else if (power == 1) {
      # reticulate::source_python(file = lp_python_path)
      pycot <- reticulate::import_from_path("pycot", path = pycot_path, convert = TRUE)
      cost <- pycot$python_lp$l1_loss
    } else {
      # reticulate::source_python(file = lp_python_path)
      # cost <- lp_loss
      cost <- paste0("Sum(Pow(X-Y,", power,"))")
      backend <- "online"
      power <- 1L
    }
  } else if (n*m > 5000^2 || backend == "multiscale" || backend == "online") {
    if (power == 2) {
      cost <- "SqDist(X,Y)"
    } else if (power == 1) {
      cost <- "Sum(Abs(X - Y))"
    } else {
      cost <- paste0("Sum(Pow(X-Y,", power,"))")
      power <- 1L
    }
    pykeops <- reticulate::import("pykeops", convert = TRUE)
    # pykeops$clean_pykeops()
  } else {
    cost <- NULL
  }
  use_cuda <- torch$cuda$is_available()
  dtype <- if(use_cuda){torch$cuda$DoubleTensor} else {torch$DoubleTensor}
  
  # sets up python data types
  xt <- dtype(np$array(x))$contiguous()
  yt <- dtype(np$array(y))$contiguous()
  at <- dtype(a)$contiguous()
  bt <- dtype(b)$contiguous()
  
  
  # sets up python function
  Loss <- geomloss$SamplesLoss("sinkhorn", p = power, blur = blur, reach = reach,
                      diameter = diameter, 
                      scaling = scaling, 
                      cost = cost, kernel = NULL,
                      cluster_scale = cluster_scale,
                      debias = debias,
                      potentials = TRUE,
                      verbose = verbose,
                      backend = backend)
  
  # run python function from geomloss library
  potentials.torch <- Loss(at, xt, bt, yt)
  
  # get potentials
  f  <- c(potentials.torch[[1]]$cpu()$numpy())
  g  <- c(potentials.torch[[2]]$cpu()$numpy())
  
  # use potentials to calculate the sinkhorn loss
  loss <- c(torch$add(torch$dot(potentials.torch[[1]]$squeeze(), at), 
                            torch$dot(potentials.torch[[2]]$squeeze(), bt))$cpu()$numpy())
  
  # setup output matrix
  retVal <- list(f = f,
                 g = g,
                 loss = loss)
  
  return(retVal)
  
}

sinkhorn_geom <- sinkhorn
