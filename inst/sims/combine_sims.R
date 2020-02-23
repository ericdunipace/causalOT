library(causalOT)

#### Output settings ####
jobids      <- as.character(c(44893116))
n           <- "512"
p           <- "6"
design      <- "A"
overlap     <- "high"
simsett     <- "hainmueller"
curdate     <- as.character(Sys.Date())

#### Get Data ####
path        <- file.path("Output", design, overlap, n, p)
files       <- list.files(path = path)
sel.files   <- unlist(c(sapply(jobids, grep, files, value = TRUE)))
file.holder <- vector("list", length(sel.files))
for(i in seq_along (sel.files) ) {
  file.holder[[i]] <- readRDS(file.path(path, sel.files[i]))
}

#### Combine ####
output      <- causalOT:::combine_sims(file.holder, simulation_name = simsett)

outname     <- paste(c(simsett, "des", design, "overlap", overlap, "n", n, "date", curdate), collapse="_")
saveRDS(output, file = file.path("Output", paste0(outname,".rds")))


