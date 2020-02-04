rm(list=ls())
library(causalOT)

overlap <- c("low","high")
design <- c("A","B")

no    <- length(overlap)
nd    <- length(design)
nsims <- 1000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 190954522) #seed from random.org

seed_array <- array(seeds, dim=c(nd, no, nsims),
                    dimnames = list(design = design,
                                    overlap = overlap,
                                    nsims = 1:nsims)
)
dump("seed_array", file="inst/seed_generation/seeds.Rdmped")

q("no")