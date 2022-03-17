#### Hainmueller ####

rm(list=ls())
library(causalOT)

overlap <- c("low","medium","high")
design <- c("A","B")

no    <- length(overlap)
nd    <- length(design)
nsims <- 1000

seeds <- causalOT:::seed.gen(design = design, overlap = overlap, niter = nsims, seed = 190954522) #seed from random.org

seed_array <- array(seeds, dim=c(nd, no, nsims),
                    dimnames = list(design = design,
                                    overlap = overlap,
                                    nsims = 1:nsims)
)
dump("seed_array", file="code/original_sim/seeds/hainmueller_seeds.Rdmped")


#### convergence ####
rm(list=ls())
library(causalOT)

overlap <- "none"

nsims <- 1000

seeds <- causalOT:::seed.gen(design = overlap, overlap = overlap, niter = nsims, seed = 555069822) #seed from random.org


seed_array <- seeds
dump("seed_array", file="code/original_sim/convergence_seeds.Rdmped")


#### algorithm ####
rm(list=ls())
library(causalOT)

overlap <- "none"

nsims <- 1000

seeds <- causalOT:::seed.gen(design = overlap, overlap = overlap, niter = nsims, seed = 878560115) #seed from random.org


seed_array <- seeds
dump("seed_array", file="code/original_sim/algorithm.Rdmped")

