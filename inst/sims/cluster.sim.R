rm(list=ls())

#### Load Packages ####
library(causalOT)
require(gurobi)

#### Environmental Parameters ####
arraynum <- Sys.getenv('SLURM_ARRAY_TASK_ID')
jobid <- Sys.getenv('SLURM_ARRAY_JOB_ID')
design <-  Sys.getenv('DESIGN')
overlap <- Sys.getenv('OVERLAP')

#### Set Seed ####
seed.file <- file.path("seeds.Rdmped")
source(seed.file)
seed <- seed_array[design, overlap, arraynum]

#### Sim param ####
n <- 2^9
p <- 6
nsims <- 10
distance <- c("Lp", "mahalanobis")
wass_power <- c(1,2)
ground_power <- 1:2
std_mean_diff <- seq(0,0.2, length.out = 100)
trunc <- c(0, 0.01, 0.05, 0.1, 0.2)
solver <- "gurobi"
augmentation <- match <- "both"
grid.search <- TRUE
# overlap don't set b/c environmental parameter
# design don't set b/c environmental parameter

#### get simulation functions ####
dataGen <- Hainmueller$new(n = n, p = p, 
                           design = design, overlap = overlap)

#### Simulations ####
times <- proc.time()
output <- sim.function(dataGen = dataGen,
                       nsims = nsims, 
                       ground_p = ground_power, 
                       p = wass_power,
                       grid.search = grid.search,
                       match = match,
                       augmentation = augmentation,
                       standardized.mean.difference = std_mean_diff,
                       truncations = trunc,
                       distance = distance, 
                       solver = solver,
                       seed = seed)


run.time <- (proc.time() - times)
print(run.time)

date <- gsub(" ", "_", as.name(as.character(Sys.time())))
date <- gsub(":", "=", date)
term <- paste0(c(date, ".rds"), collapse="")
otfl <- paste0(c("causalOT",design,overlap,n,p,jobid,arraynum,term),collapse="_")
otdr <- file.path("Output", design,overlap,n,p)
otfn <- file.path(otdr, otfl)
if(!dir.exists(otdr)) {
  drfl <- strsplit(otdr, "/")[[1]]
  for(fn in seq_along(drfl)) {
    curfile <- paste0(drfl[1:fn], collapse="/")
    if(dir.exists(curfile)) next
    dir.create(paste0(drfl[1:fn], collapse="/"))
  }
}
saveRDS(output, file=otfn)

warnings()

q("no")
