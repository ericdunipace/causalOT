supported.methods <- function() {
  return(c("Logistic","SBW", "RKHS", "RKHS.dose", "NNM", 
           "Wasserstein",
           "Constrained Wasserstein",
           "None"))
}

ot.methods <- function() {
  return(c("NNM", 
           "Wasserstein",
           "Constrained Wasserstein"))
}

dist.metrics <- function() {
  c("Lp", "mahalanobis", "RKHS", "sdLp")
}