seed.gen <- function(design, overlap, niter, seed) {
  
  nd <- length(design)
  no <- length(overlap)
  ni <- as.integer(niter)
  
  num.seeds <- nd*no*ni
  
  set.seed(seed)
  
  seeds.out <- sample.int(.Machine$integer.max, num.seeds, replace = FALSE)
  
  return(seeds.out)
}