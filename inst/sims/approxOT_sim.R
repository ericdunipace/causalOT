library(causalOT)

#### Output settings ####
jobids      <- as.character(45385874)
simsett     <- "approxOT"
curdate     <- as.character(Sys.Date())

#### Get Data ####
path        <- file.path("Output", "approx")
files       <- list.files(path = path)
sel.files   <- unlist(c(sapply(jobids, grep, files, value = TRUE)))
file.holder <- vector("list", length(sel.files))
for(i in seq_along (sel.files) ) {
  file.holder[[i]] <- readRDS(file.path(path, sel.files[i]))
}

output <- data.table::rbindlist(file.holder)

outname     <- paste(c(simsett, "date", curdate), collapse="_")
saveRDS(output, file = file.path("Output", paste0(outname,".rds")))
