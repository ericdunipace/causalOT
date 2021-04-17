supported.methods <- function() {
  return(c("Logistic","SBW", "SCM", 
           # "EBW",
           "CBPS",
           "RKHS", "RKHS.dose", "NNM", 
           "Wasserstein",
           "Constrained Wasserstein",
           "None"))
}

ot.methods <- function() {
  return(c("NNM"
           ,"Wasserstein"
           ,"Constrained Wasserstein"
           ,"SCM"
           # ,"EBW"
           ))
}

dist.metrics <- function() {
  c("sdLp", "mahalanobis",  "Lp", "RKHS" )
}