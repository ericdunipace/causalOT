library(causalOT)

#### Output settings ####
jobids      <- as.character(62255860) #-Ahi : 59681680, -Alo: 59407198,59680922, -Bhi:59158432, -Blo:59134039, Blo update 62255860
n           <- "512"
p           <- "6"
design      <- "B"
overlap     <- "lown"
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
output      <- causalOT:::combine_sims(file.holder, simulation_name = paste(simsett, "design", design, "overlap", overlap, sep="_"))

outname     <- paste(c(simsett, "des", design, "overlap", overlap, "n", n, "date", curdate), collapse="_")
saveRDS(output, file = file.path("Output", paste0(outname,".rds")), compress = FALSE)


