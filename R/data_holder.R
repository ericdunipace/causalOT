generate_holder <- function(outcome = TRUE, ...) {
  # type = c("list","data.frame")
  # type <- match.arg(type)
  dots <- list(...)
  options <- get_holder_options()
  
  estimates <- options$estimates
  dr <- options$dr
  matched <- options$matched
  weights <- options$weights
  
  weight.list <- sapply(weights, function(i){NULL}, simplify = FALSE)
  
  if(outcome) {
    m.list <- sapply(matched, function(i){weight.list}, simplify=FALSE)
    dr.list <- sapply(dr, function(i){m.list}, simplify=FALSE)
    output <- sapply(estimates, function(i) {dr.list}, simplify=FALSE)
  } else {
    output <- sapply(estimates, function(i) {weight.list}, simplify = FALSE)
  }
  if(length(dots)) {
    temp_vec <- sapply(dots, function(i) {as.character(i)}, simplify=FALSE)
    for(i in seq_along(dots)) {
      holder <- output
      output <- sapply(temp_vec[[i]], function(i) {holder}, simplify = FALSE)
    }
  } 
  class(output) <- "outputHolder"
  return(output)
}

get_holder_options <- function() {
  estimates <- c("ATT", "ATC", "cATE", "ATE", "feasible")
  dr <- c("Hajek", "DR Hajek")
  weights <- supported.methods()
  matched <- c("Unmatched", "Matched")
  return(list(estimates = estimates, dr = dr, weights = weights, matched = matched))
}

convert_holder <- function(x, outcome = TRUE, addl.terms = NULL) {
  stopifnot(inherits(x, "outputHolder"))
  options <- get_holder_options()
  
  find.names <- function(x, lookup) {
    out <- rep(NA, length(x))
    for(i in lookup) {
      idx <- which(grepl(i, x))
      out[idx] <- i
    }
    return(out)
  }
  
  estimates <- c("Naive", options$estimates)
  dr <- options$dr
  weights <- options$weights
  matched <- options$matched
  
  un <- unlist(x)
  poss.names <- names(un)
  # dr.vals <- grep("DR", dr, value = TRUE)
  # hj.vals <- grep("Hajek", dr, value = TRUE)
  
  dr.tf <- hj.tf <- NULL
  est.vec <- factor(find.names(poss.names, estimates), levels = estimates)
  weight.vec <- factor(find.names(poss.names, sort(weights, decreasing = TRUE)), levels = weights)
  
  if(outcome) {
    m.tf  <- !grepl("Unmatched", poss.names)
    dr.tf <- grepl("DR", poss.names)
    hj.tf <- grepl("Hajek", poss.names)
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec, 
                     hajek = hj.tf,  doubly.robust = dr.tf,
                     matched = m.tf)
  } else {
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec)
  }
  fac <- df$weighting.method
  df$weighting.method <- factor(ifelse(is.na(fac), "None", as.character(fac)), levels = c(levels(fac), "None"))
  
  if(!is.null(addl.terms)) {
    labels <- addl.terms$labels
    vals <- addl.terms$values
    nat <- names(labels)
    for(i in seq_along(labels)) {
      df[[nat[i]]] <- factor(find.names(poss.names, labels[[i]]), levels = labels[[i]], labels = vals[[i]])
    }
  }
  
  return(df)
}

convert_holder_sbw <- function(x, outcome = TRUE, addl.terms = NULL) {
  stopifnot(inherits(x, "outputHolder"))
  options <- get_holder_options()
  
  find.names <- function(x, lookup) {
    out <- rep(NA, length(x))
    for(i in lookup) {
      idx <- which(grepl(i, x))
      out[idx] <- i
    }
    return(out)
  }
  
  estimates <- c("Naive", options$estimates)
  dr <- options$dr
  weights <- c("SBW_j", "SBW")
  matched <- options$matched
  
  un <- unlist(x)
  poss.names <- names(un)
  # dr.vals <- grep("DR", dr, value = TRUE)
  # hj.vals <- grep("Hajek", dr, value = TRUE)
  
  dr.tf <- hj.tf <- NULL
  est.vec <- factor(find.names(poss.names, estimates), levels = estimates)
  weight.vec <- factor(find.names(poss.names, sort(weights, decreasing = FALSE)), levels = weights)
  
  if(outcome) {
    m.tf  <- !grepl("Unmatched", poss.names)
    dr.tf <- grepl("DR", poss.names)
    hj.tf <- grepl("Hajek", poss.names)
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec, 
                     hajek = hj.tf,  doubly.robust = dr.tf,
                     matched = m.tf)
  } else {
    df <- data.frame(values = un, estimate = est.vec, 
                     weighting.method = weight.vec)
  }
  fac <- df$weighting.method
  df$weighting.method <- factor(ifelse(is.na(fac), "None", as.character(fac)), levels = c(levels(fac), "None"))
  
  if(!is.null(addl.terms)) {
    labels <- addl.terms$labels
    vals <- addl.terms$values
    nat <- names(labels)
    for(i in seq_along(labels)) {
      df[[nat[i]]] <- factor(find.names(poss.names, labels[[i]]), levels = labels[[i]], labels = vals[[i]])
    }
  }
  
  return(df)
}


data_sim_holder <- function(design, overlap) {
  des.list <- sapply(design, function(i){NULL}, simplify=FALSE)
  output <- sapply(overlap, function(i) {des.list}, simplify = FALSE)
  return(output)
}

sim_holder <- function(design, overlap) {
  # smd.list <- sapply(as.character(std_mean_diff), function(i){NULL}, simplify=FALSE)
  # plist    <- sapply(power, function(i){smd.list}, simplify=FALSE)
  # dlist    <- sapply(distance, function(i){plist}, simplify=FALSE)
  des.list <- sapply(design, function(i){NULL}, simplify=FALSE)
  output <- sapply(overlap, function(i) {des.list}, simplify = FALSE)
  return(output)
}

check.names <- function(...) {
  names <- lapply(..., colnames)
  test.name <- names[[1]]
  stopifnot(all(sapply(names, all.equal, test.name)))
}

check.all.df <- function(...) {
  stopifnot(all(sapply(..., inherits, "data.frame")))
}

combine_holder_df <- function(..., simulation_name = NULL) {
  
  check.names(...)
  check.all.df(...)
  
  new_df <- do.call("rbind", ...)
  if(!is.null(simulation_name)) {
    simulation_name <- as.character(simulation_name)
    new_df$simulation.setting <- simulation_name
  }
  
  return(new_df)
}

combine_sims <- function(..., simulation_name = NULL) {
  
  if(is.list(...)) {
    dots <- (...)
  } else {
    dots <- list(...)
  }
  slot.names <- names(dots[[1]])
  length.slots <- length(slot.names)
  
  outlist <- vector("list", length.slots)
  names(outlist) <- slot.names
  temp_list <- list(NULL)
  
  for(i in 1:length.slots) {
    temp_list[[1]] <- lapply(dots, function(out) out[[slot.names[i]]])
    if(slot.names[i] != "RNG") {
      outlist[[i]] <- do.call("combine_holder_df", 
                            list(temp_list[[1]], simulation_name = simulation_name))
    } else {
      outlist[[i]] <- unlist(temp_list[[1]], recursive = FALSE)
    }
  }
  
  return(outlist)
}
