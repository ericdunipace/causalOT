supported.methods <- function() {
  return(c("Logistic","SBW", "SCM",
           "RKHS", "RKHS.dose", "NNM", 
           "Wasserstein",
           "Constrained Wasserstein",
           "None"))
}

ot.methods <- function() {
  return(c("NNM", 
           "Wasserstein",
           "Constrained Wasserstein",
           "SCM"))
}

dist.metrics <- function() {
  c("Lp", "mahalanobis", "RKHS", "sdLp")
}