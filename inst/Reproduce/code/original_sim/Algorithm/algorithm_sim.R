library(causalOT)
library(doRNG)

#### cluster param ####
ARRAYID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
JOBID <- Sys.getenv('SLURM_ARRAY_JOB_ID')

outdir <- file.path("Output", "algorithm",JOBID)
DATE   <- Sys.Date()

#### set seed ####
seed.file <- file.path("seeds/algorithm_seeds.Rdmped")
source(seed.file)
set.seed(seed_array[ARRAYID])

#### function ####
bp_bs <- function(weights, x,z, cost, estimand = "ATE", n.boot = 100, lambda = .5) {
  
  eval_fun <- function(x1, x2, w, a, b, cost, lambda = 0.5) {
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    
    a_bs <- causalOT:::renormalize(w * c( rmultinom(1, n1, a) ) )
    b_bs <- causalOT:::renormalize(b * c( rmultinom(1, n2, b) ) )
    
    epsilon <- lambda/median(cost)
    
    tplan <- approxOT::transport_plan_given_C(mass_x = a_bs, mass_y = b_bs, p = 2, cost = cost, method = "sinkhorn",
                                     debias = FALSE, epsilon = epsilon)
    
    # create normalized version of gamma for matrix sums
    gamma_norm <- matrix(0, n1, n2)
    gamma_norm[causalOT:::dist_2d_to_1d(tplan$from, tplan$to,n1, n2)] <- tplan$mass * 1/b_bs[tplan$to]
    
    loss <- sum((crossprod(x1, gamma_norm) - t(x2))^2 %*% (b*(b_bs>0)))
    return(loss)
  }
  
  x1 <- x[z == 1, ]
  x0 <- x[z == 0,]
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  n  <- nrow(x)
  # lambda1 <- weights$args$constraint[[2]]$penalty
  # lambda0 <- weights$args$constraint[[1]]$penalty
  lambda1 <- lambda0 <- lambda
  w1 <- weights$w1
  w0 <- weights$w0
  
  sw <- causalOT:::get_sample_weight(NULL, z)
  
  # x1_bs <- lapply(1:nboot, function() {sw$a * rmultinom(1, n1, sw$a) })
  # x0_bs <- lapply(1:nboot, function() {sw$b * rmultinom(1, n0, sw$b) })
  # x0_bs <- lapply(1:nboot, function() {sw$tot * rmultinom(1, n, sw$tot) })
  # 
  x1_evals <- mean(replicate(n.boot, eval_fun(x1,x,w1, sw$b, sw$total, cost[[2]], lambda = lambda1 )))
  x0_evals <- mean(replicate(n.boot, eval_fun(x0,x,w0, sw$a, sw$total, cost[[1]], lambda = lambda0 )))
  return(list(x0 = x0_evals, x1 = x1_evals))
}

sink_bs <- function(weights, x,z, cost, estimand = "ATE", n.boot = 100) {
  
  eval_fun <- function(x1, x2, w, a, b, cost) {
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    
    a_bs <- causalOT:::renormalize(w * c( rmultinom(1, n1, a) ) )
    b_bs <- causalOT:::renormalize(b * c( rmultinom(1, n2, b) ) )
    
    # tplan <- approxOT::transport_plan_given_C(mass_x = a_bs, mass_y = b_bs, p = 2, cost = cost, method = "sinkhorn",
    #                                           unbiased = FALSE)
    # gamma <- matrix(0, n1, n2)
    # gamma[causalOT:::dist_2d_to_1d(tplan$from, tplan$to,n1, n2)] <- tplan$mass
    # 
    # blur <- 0.05*median(cost)
    blur <- 100
    loss <- sinkhorn(x = x1, y = x2, a = a_bs, b = b_bs, power = 2,blur = blur, debias = TRUE)
    
    # gamma_check <- exp((matrix(loss$f, n1,n2) + matrix(loss$g,n1,n2, byrow = TRUE) - cost^2)/blur) *tcrossprod(a_bs,b_bs)
    # 
    # print(sum(gamma * cost^2))
    # print(sum(gamma_check * cost^2))
    # print(loss$loss)
    
    return(loss$loss)
  }
  
  x1 <- x[z == 1, ]
  x0 <- x[z == 0,]
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  n  <- nrow(x)
  
  sw <- causalOT:::get_sample_weight(NULL, z)
  
  # x1_bs <- lapply(1:nboot, function() {sw$a * rmultinom(1, n1, sw$a) })
  # x0_bs <- lapply(1:nboot, function() {sw$b * rmultinom(1, n0, sw$b) })
  # x0_bs <- lapply(1:nboot, function() {sw$tot * rmultinom(1, n, sw$tot) })
  # 
  x1_evals <- mean(replicate(n.boot, eval_fun(x1,x,weights$w1, sw$b, sw$total, cost[[2]])))
  x0_evals <- mean(replicate(n.boot, eval_fun(x0,x,weights$w0, sw$a, sw$total, cost[[1]])))
  return(list(x0 = x0_evals, x1 = x1_evals))
}

#### Setup Data ####
ns <- 2^10 #2^(5:11)
max_n <- max(ns)
lambda <- 1e2
n.boot <- 100
methods <- c(
               "Wass.1"
               # ,
               # "Wass.2"
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

scratch.dir <- file.path("scratch_alg", ARRAYID)
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
  # max_cost <- lapply(10^(-4:2), function(mm) lapply(cost, function(cc) list(penalty = mm * max(cc)^2)))
  max_cost <- lapply(10^(-2:5), function(mm) lapply(cost, function(cc) list(penalty = mm )))
  # rm(cost)
  df <- data.frame(x, z = z, y = y)
  
  ps1   <- ps1_full[1:n]
  ps0   <- ps0_full[1:n]
  h1  <- causalOT:::renormalize(1/ps1)
  h0  <- causalOT:::renormalize(1/(1-ps0))
  
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
        grid.search <- FALSE
        constraints <- max_cost
        solver <- "mosek"
        grid.length <- 7
      } else if (meth == "Wass.2") {
        f <- "~. + 0"
        add.divergence <- FALSE
        penalty <- "none"
        p <- 2
        grid.search <- TRUE
        constraints <- seq(1e-2, 6^(-0.5), length.out = 7)
        balance.constraints <- NULL
        solver <- "mosek"
        grid.length <- 7
        sbwweights <- calc_weight(data = df, 
                                  balance.covariates = colnames(x),
                                  treatment.indicator = "z",
                               outcome = "y",
                               method = "SBW", add.divergence = add.divergence, grid.search = grid.search,
                               formula = f, balance.constraints = balance.constraints,
                               solver = solver, penalty = penalty, 
                               constraint = constraints,
                               grid.length = grid.length,
                               estimand = "ATE", 
                               metric = "Lp", p = p, verbose = FALSE,
                              niter = 1e6, control = list(maxit = 1e6),
                               backend = "auto")
        
        meth <- "Wasserstein"
        f <- "~. + 0"
        balance.constraints <- sbwweights$args$constraint
        add.divergence <- TRUE
        penalty <- "entropy"
        p <- 2
        grid.search <- FALSE
        constraints <- max_cost
        grid.length <- 7
        solver <- "mosek"
        if (n > 2^10) return(NULL)
      }
      
      
      weights <- lapply(max_cost, function(lambda)
        calc_weight(data = df, 
                    balance.covariates = colnames(x), 
                    treatment.indicator = "z",
                    outcome = "y",
                    method = meth, 
                    add.divergence = add.divergence, 
                    grid.search = grid.search,
                    cost = cost,
                   formula = f, 
                   balance.constraints = balance.constraints,
                   solver = solver, penalty = penalty, 
                   # cost = cost,
                   search = "LBFGS",
                   constraint = lambda,
                    estimand = "ATE", metric = "Lp", p = p, verbose = FALSE,
                    niter = 1e6, control = list(maxit = 1e6),
                    backend = "auto")
      )
      
      
      est.aug <- lapply(weights, function(w)
        estimate_effect(data = df, weights = w,
                             balance.covariates = colnames(x), treatment.indicator = "z",
                             outcome = "y",
                             doubly.robust = TRUE, 
                        split.model = TRUE, matched = FALSE)
        )
      
      
      est <- lapply(weights, function(w)
        estimate_effect(data = df, weights = w,
                             balance.covariates = colnames(x),
                        treatment.indicator = "z",
                             outcome = "y",
                             doubly.robust = FALSE, split.model = TRUE, matched = FALSE)
      )
      
      E_Y1 <- sapply(est, function(x) x$variance.components$E_Y1)
      E_Y0 <- sapply(est, function(x) x$variance.components$E_Y0)
      E_Y1.aug <- sapply(est.aug, function(x) x$variance.components$E_Y1)
      E_Y0.aug <- sapply(est.aug, function(x) x$variance.components$E_Y0)
      
      estimate <- sapply(est, function(x) x$estimate)
      estimate.aug <- sapply(est.aug, function(x) x$estimate)
      
      
      
      args <- list(data = df, 
                   x0 = x[z==0,],
                   x1 = x[z==1,],
                   x = x,
                   z = z,
                   grid = max_cost,  
                   n.boot = n.boot,
                   K = 10, 
                   R = 10,
                   eval.method = "bootstrap",
                   wass.method = "sinkhorn_geom",
                   wass.iter = 1e3,
                   lambda = lambda,
                   sample_weight = causalOT:::get_sample_weight(NULL, z),
                   estimand = "ATE", 
                   method = meth, 
                   solver = solver, 
                   metric = "Lp",
                   unbiased = TRUE,
                   p = 2, 
                   cost = cost, 
                   add.joint = TRUE,
                   add.margins = FALSE, 
                   add.divergence = add.divergence,
                   joint.mapping = FALSE,
                   neg.weights = FALSE,
                   cgd = FALSE, 
                   verbose = FALSE,
                   balance.covariates = colnames(x),
                   treatment.indicator = "z",
                   outcome = "y")
      balcheck <- causalOT:::eval_weights(weights, args)
      sel0 <- balcheck$sel$control[[1]]
      sel1 <- balcheck$sel$treated[[1]]
      # start <- proc.time()
      bpcheck <- lapply(weights, function(w) bp_bs(w, x, z, cost, estimand = "ATE", n.boot = 100))
      # end <- proc.time()
      # print(end - start)
      
      sel1 <- which.min(sapply(bpcheck, function(b) b[[2]]))[1]
      sel0 <- which.min(sapply(bpcheck, function(b) b[[1]]))[1]
      
      # start <- proc.time()
      # sinkcheck <- lapply(weights, function(w) sink_bs(w, x, z, cost, estimand = "ATE", n.boot = 100))
      # end <- proc.time()
      # print(end - start)
      
      w_1 <- sapply(weights, function(w)
        causalOT::sinkhorn(x = x1, y = x1, a = h1, b = w$w1, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 1e-1,
                                     backend = "auto")$loss)
      
      w_0 <- sapply(weights, function(w)
        causalOT::sinkhorn(x = x0, y = x0, a = h0, b = w$w0, power = 2,
                                     metric = "Lp", debias = TRUE, blur = 1e-1,
                                     backend = "auto")$loss)
      
      
      l2_c = sapply(weights, function(w) sum((h0 - w$w0)^2/h0^2))
      l2_t = sapply(weights, function(w) sum((h1 - w$w1)^2/h1^2))
      and_c = sapply(weights, function(w) mean((h0 - w$w0)^2/(h0 * (1-h0))))
      and_t = sapply(weights, function(w) mean((h1 - w$w1)^2/(h1 * (1-h1))))
      
      output.meth <- data.frame(n1 = n1, n0 = n0, n = N,
                          lambda_0 = sapply(max_cost,
                                            function(x) x[[1]]$penalty),
                          lambda_1 = sapply(max_cost,
                                            function(x) x[[2]]$penalty),
                           ot0 = balcheck$ot$control,
                           ot1 = balcheck$ot$treated,
                           sel0 = sel0,
                           sel1 = sel1,
                           E_Y1 = E_Y1,
                           E_Y0 = E_Y0,
                           E_Y1.aug = E_Y1.aug,
                           E_Y0.aug = E_Y0.aug,
                           estimate = estimate,
                           estimate.aug = estimate.aug,
                          l2_c = l2_c,
                          l2_t = l2_t,
                          and_c = and_c,
                          and_t = and_t,
                          w_0 = w_0,
                          w_1 = w_1,
                           method = meth, bf = switch(orig.meth,
                                                      "Wass.1" = FALSE,
                                                      "Wass.2" = TRUE),
                           arrayid = ARRAYID)
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

#q("no")
