rm(list=ls())

set.seed(NULL)

ns <- 2^(7:13)
ds <- 2:10
arraynum <- Sys.getenv('SLURM_ARRAY_TASK_ID')
jobid <- Sys.getenv('SLURM_ARRAY_JOB_ID')

dat.gen <- function(n, d) {
  x0 <- matrix(rnorm(n*d),d,n)
  x1 <- matrix(rnorm(n*d),d,n)
  mass <- rep(1/n, n)
  cost <- causalOT::cost_calc_lp(x0, x1, 2, "colwise")
  wass <-  transport::wasserstein(mass, mass, costm= cost, p = 2)
  hilb <- approxOT::wasserstein(X = x0, Y = x1, p = 2, ground_p = 2, observation.orientation = "colwise",
                                method = "hilbert")
  sink <- approxOT::wasserstein(X = x0, Y = x1, p = 2, ground_p = 2, observation.orientation = "colwise",
                                method = "greenkhorn", epsilon = 0.01, niter = 1000)
  return(data.frame(n = n, d= d, value = c(wass, hilb, sink), 
                    method = c("wasserstein","hilbert","sinkhorn")))
}

run.fun <- function() {
  return(data.table::rbindlist(mapply(dat.gen, n = rep(ns, each = length(ds)), d = ds, SIMPLIFY=FALSE)))
}

output <- run.fun()

date <- gsub(" ", "_", as.name(as.character(Sys.time())))
date <- gsub(":", "=", date)
term <- paste0(c(date, ".rds"), collapse="")
otfl <- paste0(c("approxOT",jobid,arraynum,term),collapse="_")
otdr <- file.path("Output", "approx")
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

