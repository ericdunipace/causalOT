log_sum_exp <- function(x) {
  # if(is.vector(x)) {
  if(all(is.infinite(x))) return(x[1])
  mx <- max(x)
  x_temp <- x - mx
  return(log(sum(exp(x_temp)))+ mx)
  # } else if (is.matrix(x)) {
  #   mx <- apply(x, 1, max)
  #   x_temp <- x - mx
  #   return(log(rowSums(exp(x_temp)))+ mx)
  # }
}

log_sum_exp2 <- function(x,y) {
  mx <- pmax(x,y)
  # if(is.infinite(mx)) return(mx)
  
  temp <- cbind(x,y) - mx
  temp[mx == -Inf,] <- -Inf
  return(log(rowSums(exp(temp))) + mx)
}

renormalize <- function(x) {
  if (all(is.na(x))) return(x)
  
  if (isTRUE(any(x < 0))) {
    # warning("Negative weights found! Normalizing to sum to 1 with less accurate function. Make sure negative weights make sense for your problem")
    return(x/sum(x, na.rm = TRUE))
  }
  if (isTRUE(all(x == 0)) ) return(rep(0, length(x)))
  l_x <- log(x)
  return(exp(l_x - log_sum_exp(l_x)))
}

simplex_proj <- function(y) { #simplex projection of Condat 2015
  N <- length(y)
  v <- v_tilde <- rep(NA_real_, N)
  v_count <- 1
  vt_count <- 0
  v[1] <- y[1]
  rho <- y[1] - 1
  
  for(n in 2:N) {
    if(y[n] > rho){
      rho <- rho + (y[n] - rho)/(v_count + 1)
      if(rho > y[n] - 1) {
        v[v_count + 1] <- y[n]
        v_count <- v_count + 1
      } else {
        v_tilde[(vt_count+1):(v_count + vt_count)] <- v[1:v_count]
        vt_count <- vt_count + v_count
        v[[1]] <- y[n]
        v[2:N] <- NA_real_
        rho <- y[n] - 1
        v_count <- 1
      }
    }
  }
  if(!all(is.na(v_tilde))) { #ie, output non-empty
    v_tilde <- v_tilde[!is.na(v_tilde)]
    for(x in v_tilde) {
      if(x > rho) {
        v[[v_count]] <- x
        v_count <- v_count + 1
        rho <- rho + (x - rho)/v_count
      }
    }
  }
  change <- 1
  v_count <- sum(!is.na(v))
  while(change == 1) {
    change <- 0
    v <- v[!is.na(v)]
    for(n in 1:length(v)) {
      x <- v[n]
      if(x <= rho) {
        v[[n]] <- NA_real_
        v_count <- v_count - 1
        rho <- rho + (rho - x)/v_count
        change <- 1
      }
    }
  }
  tau <- rho
  K <- sum(!is.na(v))
  x <- pmax(y - tau, 0)
  return(x)
}
