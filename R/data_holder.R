generate_holder <- function(outcome = TRUE) {
  # type = c("list","data.frame")
  # type <- match.arg(type)
  
  options <- get_holder_options()
  
  estimates <- options$estimates
  dr <- options$dr
  weights <- options$weights
  
  weight.list <- sapply(weights, function(i){NULL}, simplify = FALSE)
  
  if(outcome) {
    dr.list <- sapply(dr, function(i){weight.list}, simplify=FALSE)
    output <- sapply(estimates, function(i) {dr.list}, simplify=FALSE)
  } else {
    output <- sapply(estimates, function(i) {weight.list}, simplify = FALSE)
  }
  class(output) <- "outputHolder"
  return(output)
}

get_holder_options <- function() {
  estimates <- c("ATT", "ATC", "ATE", "feasible")
  dr <- c("Hajek", "DR Hajek")
  weights <- c("SBW", "Constrained Wasserstein", "Wasserstein")
  return(list(estimates = estimates, dr = dr, weights = weights))
}

convert_holder <- function(x, outcome = TRUE) {
  stopifnot(inherits(x, "outputHolder"))
  options <- get_holder_options()
  
  find.names <- function(x, lookup) {
    out <- rep(NA, length(x))
    for(i in lookup) {
      idx <- which(grepl(i, x))
      out[idx] <- i
    }
    return(out)
  }
  
  estimates <- c("Naive", options$estimates)
  dr <- options$dr
  weights <- options$weights
  
  un <- unlist(x)
  poss.names <- names(un)
  # dr.vals <- grep("DR", dr, value = TRUE)
  # hj.vals <- grep("Hajek", dr, value = TRUE)
  
  dr.tf <- hj.tf <- NULL
  est.vec <- factor(find.names(poss.names, estimates), levels = estimates)
  weight.vec <- factor(find.names(poss.names, sort(weights, decreasing = TRUE)), levels = weights)
  
  if(outcome) {
    dr.tf <- grepl("DR", poss.names)
    hj.tf <- grepl("Hajek", poss.names)
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec, 
                     hajek = hj.tf,  doubly.robust = dr.tf)
  } else {
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec)
  }
  fac <- df$weighting.method
  df$weighting.method <- factor(ifelse(is.na(fac), "None", as.character(fac)), levels = c(levels(fac), "None"))
  
  return(df)
}

data_sim_holder <- function(design, overlap) {
  des.list <- sapply(design, function(i){NULL}, simplify=FALSE)
  output <- sapply(overlap, function(i) {des.list}, simplify = FALSE)
  return(output)
}

sim_holder <- function(design, overlap, distance, power, std_mean_diff) {
  smd.list <- sapply(as.character(std_mean_diff), function(i){NULL}, simplify=FALSE)
  plist    <- sapply(power, function(i){smd.list}, simplify=FALSE)
  dlist    <- sapply(distance, function(i){plist}, simplify=FALSE)
  des.list <- sapply(design, function(i){dlist}, simplify=FALSE)
  output <- sapply(overlap, function(i) {des.list}, simplify = FALSE)
  return(output)
}
