rm(list=ls())

node <- Sys.getenv('SLURM_JOB_NODELIST')
message(node)

#### Load Packages ####
library(causalOT)

#### Environmental Parameters ####
arraynum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
jobid <- Sys.getenv('SLURM_ARRAY_JOB_ID')
arrayset <- as.numeric(Sys.getenv('ARRAY_SET'))

arraylookup_master <- readRDS(file = "sim_arraylookup.rds")
n <- as.numeric(Sys.getenv('NOBS'))

arraylookup <- arraylookup_master[arraylookup_master$n == n & arraylookup_master$arrayset == arrayset, , drop = FALSE]

data <- as.character(arraylookup$data[arraynum])
overlap <- as.character(arraylookup$overlap[arraynum])
design <- as.character(arraylookup$design[arraynum])
p <- as.numeric(arraylookup$p[arraynum])
expernum <- as.numeric(arraylookup$experiment.number[arraynum])
penalty <- as.character(arraylookup$penalty[arraynum])
metric <- as.character(arraylookup$metric[arraynum])
formula <- as.character(arraylookup$formula[arraynum])
methods <- as.character(arraylookup$methods[[arraynum]][[1]])

methods <- c("Logistic", "CBPS", "SBW", "NNM", "Wasserstein", "SCM")
penalty <- c("entropy","L2")
metric <- "sdLp"

#### Set Seed ####
  seed.file <- file.path("seeds/hainmueller_seeds.Rdmped")
  source(seed.file)
  seed <- seed_array[design, overlap, expernum]

#### Sim param ####
nsims <- 4
wass_power <- c(1,2)
std_mean_diff <- seq(1e-04, p^(-1/2), length.out = 10)
std_mean_diff <- 1e-04
trunc <- 0
solver <- "mosek"
augmentation <- match <- "both"
grid.search <- TRUE
grid.search <- FALSE
RKHS <- list(opt = FALSE, opt.method = "stan",
             kernel = "polynomial")
calc.feasible <- FALSE #FALSE
cwass.targ <- c("SBW")

psform  <- if (choose(p,2) < n) {
  list(Logistic = list("z~."),
       SBW = list("~. + 0"),
       Wasserstein = list(NA, "~. + 0"),
       "Constrained Wasserstein" = list(NA, "~. + 0"))
} else {
  list(Logistic = list("z~."),
       SBW = list("~. + 0"),
       Wasserstein = list(NA, "~. + 0"),
       "Constrained Wasserstein" = list(NA, "~. + 0"))
}
outcome.models <- c("lm")
# overlap don't set b/c environmental parameter
# design don't set b/c environmental parameter
reticulate::use_condaenv("COT", require = TRUE)

#### get simulation functions ####
  dataGen <- Hainmueller$new(n = n, p = p,
                             design = design, overlap = overlap)

message("Dataset: ", data)
message("N: ", n)
message("P: ", p)
message("overlap: ", overlap)
message("design: ", design)
message("formula: ", formula)
message("methods: ", paste(methods, collapse = ", "))
message("penalty: ", paste(penalty, collapse = ", "))

#### Simulations ####
times <- proc.time()
output <- causalOT:::sim.function(dataGen = dataGen,
                       methods = methods,
                       nsims = nsims,
                       p = wass_power,
                       grid.search = grid.search,
                       RKHS = RKHS,
                       match = match,
                       propensity.formula = psform,
                       outcome.model = outcome.models,
                       estimands = c("ATE"),
                       augmentation = augmentation,
                       standardized.mean.difference = std_mean_diff,
                       truncations = trunc,
                       distance = metric, #distance,
                       solver = solver,
                       calculate.feasible = calc.feasible,
                       constrained.wasserstein.target = cwass.targ,
                       seed = seed, wass.method = "sinkhorn",
                       wass.niter = 1e4,
                       add.joint = TRUE, add.margins = c(FALSE),
                       joint.mapping = c(FALSE),
                       add.divergence = c(TRUE, FALSE),
                       neg.weights = FALSE, #c(TRUE, FALSE),
                       penalty = penalty,
                       confidence.interval = "asymptotic",
                       verbose = TRUE)


run.time <- (proc.time() - times)
print(run.time)

#### Save data ####
date <- gsub(" ", "_", as.name(as.character(Sys.time())))
date <- gsub(":", "=", date)
term <- paste0(c(date, ".rds"), collapse="")
otfl <- paste0(c("causalOT",data,design,overlap,n,p,jobid,expernum,arraynum,term),collapse="_")
otdr <- file.path("Output", data, design, overlap, n, p)
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
