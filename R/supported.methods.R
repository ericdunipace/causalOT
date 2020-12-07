supported.methods <- function() {
  return(c("Logistic","SBW", "RKHS", "RKHS.dose", "NNM", 
           "Wasserstein",
           "Constrained Wasserstein"))
}

ot.methods <- function() {
  return(c("NNM", 
           "Wasserstein",
           "Constrained Wasserstein"))
}