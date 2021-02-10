#### Hainmueller ####

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
dump("seed_array", file="inst/seed_generation/hainmueller_seeds.Rdmped")


#### Sonabend ####
rm(list=ls())
library(causalOT)

overlap <- c("low","high")
design <- c("A","B")

no    <- length(overlap)
nd    <- length(design)
nsims <- 1000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 366637009) #seed from random.org

seed_array <- array(seeds, dim=c(nd, no, nsims),
                    dimnames = list(design = design,
                                    overlap = overlap,
                                    nsims = 1:nsims)
)
dump("seed_array", file="inst/seed_generation/sonabed_seeds.Rdmped")

#### Kang and Schafer ####
rm(list=ls())
library(causalOT)

overlap <- c("low","high")
design <- c("A","B")

no    <- length(overlap)
nd    <- length(design)
nsims <- 1000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 386100177) #seed from random.org

seed_array <- array(seeds, dim=c(nd, no, nsims),
                    dimnames = list(design = design,
                                    overlap = overlap,
                                    nsims = 1:nsims)
)
dump("seed_array", file="inst/seed_generation/kangschafer_seeds.Rdmped")


#### Fan et Al ####
rm(list=ls())
library(causalOT)

overlap <- c("none")
design <- c("A","B")

no    <- length(overlap)
nd    <- length(design)
nsims <- 1000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 963017530) #seed from random.org

seed_array <- array(seeds, dim=c(nd, nsims),
                    dimnames = list(design = design,
                                    nsims = 1:nsims)
)
dump("seed_array", file="inst/seed_generation/fanetal_seeds.Rdmped")


#### LaLonde CI ####
rm(list = ls())
library(causalOT)

overlap <- c("none")
design <- c("none")

no    <- length(overlap)
nd    <- length(design)
nsims <- 10000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 517691423) #seed from random.org

seed_array <- array(seeds, dim = c(nsims),
                    dimnames = list(nsims = 1:nsims)
)
dump("seed_array", file = "inst/seed_generation/lalonde_seeds_ci.Rdmped")

#### LaLonde CI ####
rm(list = ls())
library(causalOT)

overlap <- c("none")
design <- c("none")

no    <- length(overlap)
nd    <- length(design)
nsims <- 10000

seeds <- seed.gen(design = design, overlap = overlap, niter = nsims, seed = 517691423) #seed from random.org

seed_array <- array(seeds, dim = c(nsims),
                    dimnames = list(nsims = 1:nsims)
)
dump("seed_array", file = "inst/seed_generation/lalonde_seeds_est.Rdmped")

# q("no")