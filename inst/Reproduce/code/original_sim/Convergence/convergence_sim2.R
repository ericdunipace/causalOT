# setup torch
options(timeout=600)
kind <- "cu116" 
version <- "0.10.0"
options(repos = c(
  torch = sprintf("https://storage.googleapis.com/torch-lantern-builds/packages/%s/%s/", kind, version),
  CRAN = "https://cloud.r-project.org" # or any other from which you want to install the other R dependencies.
))

# misc packages and torch
if(!require("torch")) install.packages("torch")
if(!require("remotes")) install.packages(c("remotes"))
if(!require("reticulate")) install.packages(c("reticulate"))
if(!require("foreach")) install.packages(c("foreach"))
if(!require("doRNG")) install.packages(c("doRNG"))
if(!require("doParallel")) install.packages(c("doParallel"))

# install rkeops and check
# if(!require("rkeops")) install.packages("rkeops")
# rkeops::use_gpu()

# > 2.0
#system("export CUDA_ROOT=/usr/local/cuda-11.6")
#system("export CUDA_ROOT=/usr/local/cuda-11.6")
#system("export CUDA_PATH=/usr/local/cuda-11.6")
#Sys.setenv("CUDA_PATH"="/usr/local/cuda-11.6")
#Sys.setenv("CUDA_ROOT"="locate libcuda.so.1")
if(!any(grepl("rkeops", installed.packages()[,1]))) remotes::install_github("getkeops/keops", subdir = "rkeops", force = TRUE)
# reticulate::virtualenv_create("rkeops")
# reticulate::use_virtualenv(virtualenv = "rkeops", required = TRUE)
# reticulate::py_config()
# reticulate::py_install("pykeops", method = "virtualenv", pip = TRUE)
# pykeops <- reticulate::import("pykeops")
# pykeops$test_numpy_bindings()
# stopifnot("Pykeops can't recognize CUDA" = pykeops$pykeopsconfig$gpu_available)

# install COT
if(!require("causalOT")) remotes::install_github("ericdunipace/causalOT", ref = "s4class", quick = TRUE, force = TRUE)

# load packages
library(causalOT)
library(foreach)
library(doRNG)
library(doParallel)

#setup sim
set.seed(853333764) #random.org
estimand <- "ATT"
nsims <- 100
max_n <- switch(estimand, "ATE" = 2^13, "ATT" = 2^14,
                "ATC" = 2^14)
n_seq <- 2^(4:(log(max_n)/log(2)))
n_gen <- max_n*3
n_worker <- 7L
worker_iteration <- ceiling(nsims/n_worker)
d <- 6

# function to setup matrix for SBW
# from https://stackoverflow.com/questions/49538911/r-given-a-matrix-and-a-power-produce-multiple-matrices-containing-all-unique?rq=1
fun <- function(mat,p) {
  mat <- as.data.frame(mat)
  combs <- do.call(expand.grid,rep(list(seq(ncol(mat))),p)) # all combinations including permutations of same values
  combs <- combs[!apply(combs,1,is.unsorted),]              # "unique" permutations only
  rownames(combs) <- apply(combs,1,paste,collapse="-")      # Just for display of output, we keep info of combinations in rownames
  combs <- combs[order(rownames(combs)),]                   # sort to have desired column order on output
  apply(combs,1,function(x) Reduce(`*`,mat[,x]))            # multiply the relevant columns
}

n_cuts_sbw <- cumsum(sapply(2:12, function(i) choose(d + i - 1, i)))


methods <- c(
  "COT"
  , "NNM.2"
  , "NNM.4"
  , "SBW.5" # all possible covariate combinations upto 5
  , "SBW.2" # up to 2 moments
  , "SBW.1" # just first moments
  , "Logistic", "Probit"
)

registerDoParallel(n_worker)

# run sim
cat(c("Running ", nsims, " simulations.\nThe maximum sample size is ", max_n, ".\nThe estimand is ",estimand,".\n"))
start <- proc.time()
result <- 
  foreach(nw=1:n_worker, .errorhandling = "pass",.export = c("estimand", "worker_iteration")) %dorng% {
    if(utils::packageVersion("rkeops") >= "2.0") {
      vename <- paste0("rkeops", nw)
      #Sys.setenv("CUDA_PATH"="/usr/local/cuda-11.6")
      #Sys.setenv("CUDA_ROOT"="locate libcuda.so.1")
      if(reticulate::virtualenv_exists(vename)) reticulate::virtualenv_remove(vename, confirm = FALSE)
      reticulate::virtualenv_create(vename)
      reticulate::use_virtualenv(virtualenv = vename, required = TRUE)
      reticulate::py_config()
      reticulate::py_install("pykeops", method = "virtualenv")
      rkeops::rkeops_use_gpu()
      rkeops::rkeops_use_float32()
    } else {
      rkeops::compile4gpu()
      rkeops::compile4float32()
      rkeops::use_gpu()
    }
    worker_res <- foreach(i=1:worker_iteration, .errorhandling = "pass", .packages=c('causalOT'), .export = "estimand") %do% {
      cat(c("Worker:",nw,", Iteration number", i,"of",worker_iteration,"\n"))
      
      
      datasim <- Hainmueller$new(n = n_gen, overlap = "high")
      datasim$gen_data()
      
      zt <- datasim$get_z()
      
      pscore <- c(datasim$get_pscore())
      ps1_full <- pscore[zt==1][1:max_n]
      ps0_full <- pscore[zt==0][1:max_n]
      t_idx_first <- which(zt==1)[1:max_n]
      c_idx_first <- which(zt==0)[1:max_n]
      tot_idx_first <- c(t_idx_first,c_idx_first)
      
      x_full  <- datasim$get_x()[tot_idx_first,]
      y_full  <- datasim$get_y()[tot_idx_first]
      
      t_idx   <- 1:max_n
      c_idx   <- max_n + 1:max_n
      
      
      #iterate over each sample size
      output <- foreach::foreach(n = n_seq, .combine = rbind, .export = "estimand") %do% {
        # cat(c("  sample size ", n,"\n"))
        t_sub <- t_idx[1:n]
        c_sub <- c_idx[1:n]
        tot_sub <- c(t_sub,
                     c_sub)
        x1 <- x_full[t_sub,]
        x0 <- x_full[c_sub,]
        
        x <- rbind(x1,x0)
        z  <- c(rep(1,n),
                rep(0,n))
        y <- y_full[c(t_sub,c_sub)]
        ps1   <- ps1_full[1:n]
        ps0   <- ps0_full[1:n]
        h1  <- causalOT:::renormalize(1/ps1)
        h0  <- causalOT:::renormalize(1/(1-ps0))
        
        N <- n * 2
        n0 <- n1 <- n
        b <- switch(estimand, "ATE" = rep(1/N,N), rep(1/n,n))
        torch::torch_manual_seed(sample.int(.Machine$integer.max, 1))
        cost.online <- switch(isTRUE(n > 2^11) + 1L, "tensorized", "online")
        
        output <- foreach::foreach(meth = methods, .combine = rbind, .export = "estimand") %do% {
          
          # cat(c("    Method: ", meth))
          
          x_run <- x
          run.meth <- meth
          p <- 2
          if (meth == "COT") {
            opts <- cotOptions(lambda.bootstrap = Inf, p = as.integer(p),
                               cost.online = cost.online)
          } else if (meth == "NNM.4") {
            run.meth <- "NNM"
            opts <- cotOptions(p = as.integer(4), cost.online = cost.online)
          } else if (meth == "NNM.2") {
            run.meth <- "NNM"
            opts <- cotOptions(p = as.integer(2), cost.online = cost.online)
          } else if (meth == "SBW.5") {
            if( n >= 2^9) {
              # cat("\n")
              return (NULL)
            }
            opts <- sbwOptions(grid.length = 20L, verbose = FALSE, max_iter = 10000L)
            run.meth <- "SBW"
            for(j in seq_along(n_cuts_sbw)) {
              if(n < n_cuts_sbw[j]) break
              x_run <- cbind(x_run, fun(x, j + 1))
            }
          } else if (meth == "SBW.2") {
            opts <- sbwOptions(grid.length = 20L, verbose = FALSE, max_iter = 10000L)
            run.meth <- "SBW"
            for(j in seq_along(n_cuts_sbw[1])) {
              if(n < n_cuts_sbw[j]) break
              x_run <- cbind(x_run, fun(x, j + 1))
            }
          } else if (meth == "SBW.1") {
            opts <- sbwOptions(grid.length = 20L, verbose = FALSE, max_iter = 10000L)
            run.meth <- "SBW"
          } else {
            opts <- NULL
          }
          startmeth <- proc.time()
          weights <- calc_weight(x = x_run, z = z, estimand = estimand, method = run.meth, options = opts)
          est <- estimate_effect(causalWeights = weights, x = x, y = y)
          ci <- confint(est)
          estimate <- coef(est)
          rm(est)
          
          est.aug <- estimate_effect(causalWeights = weights, x = x, y = y, 
                                     model.function = "lm", augment.estimate = TRUE)
          ci.aug <- confint(est.aug)
          estimate.aug <- coef(est.aug)
          rm(est.aug)
          
          w0 <- causalOT:::renormalize(weights@w0)
          w1 <- causalOT:::renormalize(weights@w1)
          
          rm(weights)
          
          w_1 <- ot_distance(x1 = x1, x2 = x1, a = h1, b = w1, p = 2,
                             debias = TRUE, penalty = 10, online.cost = cost.online)
          w_0 <- ot_distance(x1 = x0, x2 = x0, a = h0, b = w0, p = 2,
                             debias = TRUE, penalty = 10, online.cost = cost.online)
          if(estimand == "ATE") {
            w_xy_1 <- ot_distance(x1 = x1, x2 = x, a = w1, b = b, p = 2,
                                  debias = TRUE, penalty = 10, online.cost = cost.online)
            w_xy_0 <- ot_distance(x1 = x0, x2 = x, a = w0, b = b, p = 2,
                                  debias = TRUE, penalty = 10, online.cost = cost.online)
          } else {
            w_xy_0 <- w_xy_1 <- ot_distance(x1 = x0, x2 = x1, a = w0, b = w1, p = 2,
                                            debias = TRUE, penalty = 10, online.cost = cost.online)
            
          }
          
          
          output.meth <- data.frame(n1 = n1, n0 = n0, n = N,
                                    estimand = estimand,
                                    w_b_c = w_xy_0, w_b_t = w_xy_1,
                                    w_c = w_0, w_t = w_1,
                                    and_c = mean((h0 - w0)^2/(h0 * (1-h0))), and_t = mean((h1 - w1)^2/(h1 * (1-h1))),
                                    l2_c = sum((h0 - w0)^2), l2_t = sum((h1 - w1)^2),
                                    ci.lwr = ci[1],
                                    ci.upr = ci[2],
                                    estimate = estimate,
                                    ci.lwr.aug = ci.aug[1],
                                    ci.upr.aug = ci.aug[2],
                                    estimate.aug = estimate.aug,
                                    method = meth)
          endmeth <- proc.time()
          # cat(c(", ", endmeth[3] - startmeth[3], "\n"))
          return(output.meth)
        }
        
        return(output)
      }
      
      return(output)
      
    }
    reticulate::virtualenv_remove(vename, confirm = FALSE)
    return(worker_res)
  }

end <- proc.time()

# check sim
length(result)
print(end - start)
nsims
stopImplicitCluster()

saveRDS(result, file = "convergence.rds")

# check first
#str(result[[1]])
#str(result[[8]])
cat(c("Pct completed: ", sum(sapply(unlist(result, recursive = FALSE),is.data.frame))/n_worker/worker_iteration))

# close
#rkeops::clean_rkeops()
#reticulate::virtualenv_remove(envname = "rkeops")

