# combines raw simulations into one data.frame for each simulation

library(causalOT)

#### load and edit ####
loadEdit <- function(files, data, design, overlap , n, p) {
  out <- do.call("rbind", lapply(files, function(f) {
      temp <- readRDS(f) 
      cbind(temp$outcome, temp$Wasserstein[,c("w0","w1")])
    }
    ))
  
  out$data <- data
  out$design <- design
  out$overlap <- overlap
  out$n <- n
  out$p <- p
  return(out)
}

loadEdit_conv <- function(files, data, design, overlap , n, p) {
  out <- do.call("rbind", lapply(files, function(f) {
    readRDS(f) 
  }
  ))
  
  return(out)
}
#adopted https://stackoverflow.com/questions/54840918/how-to-unlist-nested-lists-while-keeping-vectors
g <- function(L){
  out <- unlist(L, recursive = FALSE)
  if (!all(sapply(out, is.data.frame))) {
    out <- sapply(out, g)
  }
  out
}

#### Get Data Names ####
root <- "data/raw_results"
sim_data <- list.files(root, recursive = FALSE)

output <- list()

# gets all individual simulation saves
for (s in sim_data) {
  designs <- list.files(file.path(root,s), recursive = FALSE)
  for (d in designs) {
    overlap <- list.files(file.path(root,s, d), recursive = FALSE)
    for (o in overlap) {
      sampsize <- list.files(file.path(root, s, d, o), recursive = FALSE)
      for (ss in sampsize) {
        param <- list.files(file.path(root,s, d, o, ss), recursive = FALSE)
        for (p in param) {
          output[[s]][[d]][[o]][[ss]][[p]] <- if(!grepl("convergence",s) ) {
            loadEdit(list.files(file.path(root,s, d, o, ss,p),
                                full.names = TRUE),
                     s, d, o, ss, p)
          } else {
            loadEdit_conv(list.files(file.path(root,s, d, o, ss,p),
                                full.names = TRUE),
                     s, d, o, ss, p)
          }
        }
      }
    }
  }
}

#### combine all the lists for one dataset ####
save_data <- list()
for (s in sim_data) {
  save_data[[s]] <- do.call("rbind", g(output[[s]]))
}

#### Save data ####
root.save <- "data"

for (s in sim_data) {
  saveRDS(save_data[[s]], file = file.path(root.save, paste0(s, ".rds")))
}

