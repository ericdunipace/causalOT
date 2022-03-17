library(causalOT)
library(doRNG)

#### cluster param ####
ARRAYID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
JOBID <- Sys.getenv('SLURM_ARRAY_JOB_ID')

outdir <- file.path("Output", "convergence","pen",JOBID)
DATE   <- Sys.Date()

#### set seed ####
seed.file <- file.path("seeds/convergence_seeds.Rdmped")
source(seed.file)
set.seed(seed_array[ARRAYID])


#### Setup Data ####
ns <- 2^(5:11)
max_n <- max(ns)
lambda <- 100
methods <- c(
	       "NNM", 
	       "SBW",
               "Probit",
               "Wass.1",
               "Wass.2",
               "Wass.3"
             , "Wass.4",
               "Wass.5",
               "Wass.6"
             )

dataGen <- causalOT::Hainmueller$new(n = max_n*3, p = 6, overlap = "high")
dataGen$gen_data()
zt <- dataGen$get_z()
pscore <- c(dataGen$get_pscore())
ps1_full   <- pscore[zt==1][1:max_n]
ps0_full   <- pscore[zt==0][1:max_n]
t_idx_first <- which(zt==1)[1:max_n]
c_idx_first <- which(zt==0)[1:max_n]
tot_idx_first <- c(t_idx_first,c_idx_first)

x_full  <- dataGen$get_x()[tot_idx_first,]
y_full  <- dataGen$get_y()[tot_idx_first]

t_idx   <- 1:max_n
c_idx   <- max_n + 1:max_n

sbw.bal <- NA_real_

reticulate::use_condaenv("COT")

scratch.dir <- file.path("scratch2", ARRAYID)
dir.create(scratch.dir, recursive = TRUE)

pykeops <- reticulate::import("pykeops", convert = TRUE)

pykeops$set_bin_folder(scratch.dir)

pykeops$clean_pykeops()

#### Run Simulations ####
start <- proc.time()
out <- foreach::foreach(n = ns, .combine = rbind) %do% {
  # internal.time.start <- proc.time()
  cat(c("sample size ", n,"\n"))
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
  cost <- cost_fun(x = x, z = z, p = 2, metric = "sdLp", estimand = "ATE")          
  max_cost <- lapply(10^(-4:2), function(mm) lapply(cost, function(cc) list(penalty = mm * max(cc)^2)))
  rm(cost)
  df <- data.frame(x, z = z, y = y)
  
  ps1   <- ps1_full[1:n]
  ps0   <- ps0_full[1:n]
  h1  <- causalOT::renormalize(1/ps1)
  h0  <- causalOT::renormalize(1/(1-ps0))
  
  n0  <- n1 <- n
  N   <- n * 2
  
  b <- rep(1/N,N)
  
  output <- foreach::foreach(meth = methods, .combine = rbind) %do% {
      cat(c("  Method: ", meth,"\n"))
      orig.meth <- meth
      if (meth == "Wass.1") {
        meth <- "Wasserstein"
        f <- NULL
        balance.constraints <- NULL
        add.divergence <- TRUE
        penalty <- "entropy"
        p <- 2
        grid.search <- TRUE
        constraints <- max_cost
        solver <- "mosek"
        grid.length <- 7
      } else if (meth == "Wass.2") {
        meth <- "Wasserstein"
        f <- NULL
        balance.constraints <- NULL
        add.divergence <- FALSE
        penalty <- "L2"
        p <- 4
        grid.search <- TRUE
        constraints <- NULL
        solver <- "mosek"
        grid.length <- 7
        if (n > 2^9) return(NULL)
      } else if (meth == "Wass.3") {
        meth <- "Wasserstein"
        f <- NULL
        balance.constraints <- NULL
        add.divergence <- FALSE
        penalty <- "entropy"
        p <- 2
        grid.search <- TRUE
        constraints <- NULL
        solver <- "lbfgs"
        grid.length <- 7
        if (n > 2^10) return(NULL)
      } else if (meth == "Wass.4") {
        meth <- "Wasserstein"
        f <- "~. + 0"
        balance.constraints <- sbw.bal
        add.divergence <- TRUE
        penalty <- "entropy"
        p <- 2
        grid.search <- TRUE
        constraints <- max_cost
        grid.length <- 7
        solver <- "mosek"
        if (n > 2^10) return(NULL)
      } else if (meth == "Wass.5") {
        meth <- "Wasserstein"
        f <- "~. + 0"
        balance.constraints <- sbw.bal
        add.divergence <- FALSE
        penalty <- "L2"
        p <- 4
        grid.search <- TRUE
        constraints <- NULL
        solver <- "mosek"
        grid.length <- 7
        if (n > 2^9) return(NULL)
      } else if (meth == "Wass.6") {
        meth <- "Wasserstein"
        f <- "~. + 0"
        balance.constraints <- sbw.bal
        add.divergence <- FALSE
        penalty <- "entropy"
        p <- 2
        grid.search <- TRUE
        constraints <- NULL
        solver <- "lbfgs"
        grid.length <- 7
        if (n > 2^10) return(NULL)
      } else if (meth == "NNM") {
        f <- NULL
        add.divergence <- FALSE
        penalty <- "none"
        p <- 4
        grid.search <- FALSE
        constraints <- NULL
        solver <- "mosek"
        balance.constraints <- NULL
        grid.length <- 7
      } else if (meth == "SBW") {
        f <- "~. + 0"
        add.divergence <- FALSE
        penalty <- "none"
        p <- 2
        grid.search <- TRUE
        constraints <- seq(1e-2, 6^(-0.5), length.out = 7)
        balance.constraints <- NULL
        solver <- "mosek"
        grid.length <- 7
      } else {
        f <- "z~ . "
        add.divergence <- FALSE
        penalty <- "none"
        p <- 2
        grid.search <- TRUE
        constraints <- NULL
        solver <- "mosek"
        balance.constraints <- NULL
        grid.length <- 7
      }
      
        
      weights <- calc_weight(data = df, balance.covariates = colnames(x), treatment.indicator = "z",
                                outcome = "y",
                                method = meth, add.divergence = add.divergence, grid.search = grid.search,
                                formula = f, balance.constraints = balance.constraints,
                                solver = solver, penalty = penalty, 
                                constraint = constraints,
                                grid.length = grid.length,
                                estimand = "ATE", metric = "sdLp", p = p, verbose = TRUE,
                                n.boot = 100, niter = 1e6, control = list(maxit = 1e6),
                                backend = "auto")
      
      if (meth == "SBW") sbw.bal <<- weights$args$constraint
      
      est.aug <- estimate_effect(data = df, weights = weights,
                             balance.covariates = colnames(x), treatment.indicator = "z",
                             outcome = "y",
                             doubly.robust = TRUE, split.model = TRUE, matched = FALSE)
      
      ci.aug  <- confint(est.aug, method = "asymptotic", model = "lm",
      					formula = list(treated = "y ~ .",
                                       control = "y ~ ."))
      estimate.aug <- est.aug$estimate
      rm(est.aug)
      
      est <- estimate_effect(data = df, weights = weights,
                             balance.covariates = colnames(x), treatment.indicator = "z",
                             outcome = "y",
                             doubly.robust = FALSE, split.model = TRUE, matched = FALSE)
      
      ci  <- confint(est, method = "asymptotic", model = "lm",
      					formula = list(treated = "y ~ .",
                                       control = "y ~ ."))
      estimate<- est$estimate
      rm(est)
        
      w0 <- renormalize(weights$w0)
      w1 <- renormalize(weights$w1)
      rm(weights)
      
      n0  <- length(w0)
      n1  <- length(w1)
      w_1 <- causalOT::sinkhorn(x = x1, y = x1, a = h1, b = w1, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 10,
                                     backend = "auto")$loss
      w_0 <- causalOT::sinkhorn(x = x0, y = x0, a = h0, b = w0, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 10,
                                     backend = "auto")$loss
      w_xy_1 <- causalOT::sinkhorn(x = x1, y = x, a = w1, b = b, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 10,
                                     backend = "auto")$loss
      w_xy_0 <- causalOT::sinkhorn(x = x0, y = x, a = w0, b = b, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 10,
                                     backend = "auto")$loss
      
      output.meth <- data.frame(n1 = n1, n0 = n0, n = N,
                           w_b_c = w_xy_0, w_b_t = w_xy_1,
                           w_c = w_0, w_t = w_1,
                           and_c = mean((h0 - w0)^2/(h0 * (1-h0))), and_t = mean((h1 - w1)^2/(h1 * (1-h1))),
                           l2_c = sum((h0 - w0)^2), l2_t = sum((h1 - w1)^2),
                           ci.lwr = ci$CI[1],
                           ci.upr = ci$CI[2],
                           estimate = estimate,
                           ci.lwr.aug = ci.aug$CI[1],
                           ci.upr.aug = ci.aug$CI[2],
                           estimate.aug = estimate.aug,
                           method = meth, penalty = penalty)
      return(output.meth)
  }
  
  return(output)
}
end <- proc.time()
print(end - start)

#### Save output ####
if(!dir.exists(outdir)) {
  dir.create(path = outdir, recursive = TRUE)
}
saveRDS(out, file = file.path(outdir, paste0("pen_conv_",ARRAYID,"_",DATE,".rds")))

print(warnings())

q("no")
