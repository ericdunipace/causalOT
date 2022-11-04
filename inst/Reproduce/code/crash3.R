
node <- Sys.getenv('SLURM_JOB_NODELIST')
message(node)

#### load packages ####
library(causalOT)
library(dplyr)

#### Environmental params
arraynum <- as.numeric(Sys.getenv('SLURM_TASK_ARRAY_ID'))

#### set.seed ####
set.seed(838697155) #from random.org

#### evaluation functions ####

sample.eval <- function(data, 
                        power = power,
                        mapped = FALSE)
{
  
  ns <- data$get_n()
  n0 <- ns[1]
  n1 <- ns[2]
  
  weights <- list()
  weights$Naive <- calc_weight(data = data, 
                               estimand = "ATE", 
                               method = "None",
                               grid.search = TRUE)
  weights$IPW  = calc_weight(data = data, 
                             estimand = "ATE", 
                             method = "Logistic",
                             grid.search = TRUE,
                             # formula = "z ~ .",
                             formula = "z ~ . + .*. + I(.^2)")
  
  weights$CBPS  = calc_weight(data = data, 
                             estimand = "ATE", 
                             method = "CBPS",
                             grid.search = TRUE,
                             # formula = "z ~ .",
                             formula = "z ~ . + .*. + I(.^2)")
  
  weights$SBW  = tryCatch(calc_weight(data = data,
                                      estimand = "ATE",
                                      method = "SBW",
                                      grid.search = TRUE,
                                      solver = "mosek",
                                      grid.length = 20,
                                      grid = c(seq(0, 0.5, length.out = 20), 1, 2,10),
                                      formula = "z ~ . + .*. + I(.^2) + .*.*. + .*.*.*."), 
                          error = function(e) {warning(e$message); return(list(w0 = rep(NA_real_, n0), w1 = rep(NA_real_, n1))) })
  weights$NNM = calc_weight(data = data, estimand = "ATE", method = "NNM",
                            p = power,
                            metric = "sdLp")
  
  weights$COT =  tryCatch(calc_weight(data = data,
                                     estimand = "ATE",
                                     method = "Wasserstein",
                                     add.divergence = TRUE,
                                     grid.search = TRUE,
                                     p = power, grid.length = 7,
                                     metric = "sdLp",
                                     wass.method = "sinkhorn",
                                     n.boot = 100,
                                     backend = "tensorized",
                                     verbose = TRUE),
                         error = function(e) {
                           warning(e$message)
                           return(list(w0 = rep(NA_real_, n0),
                                       w1 = rep(NA_real_, n1)))
                         })
  
  
  
  hajek <- lapply(weights,
    function(w) {
        return(estimate_effect(data = data, weights = w, 
                               matched = FALSE, 
                               split.model = TRUE,
                               doubly.robust = FALSE,
                               estimand = "ATE",
                               treatment.indicator = "z", 
                               balance.covariates = balance.covariates,
                               outcome = "y"))
    })
  
  linear.est <- sapply(hajek, function(h) h$estimate)
  linear.confint <- lapply(hajek, confint)
  linear.ci <- list(lwr = sapply(linear.confint, function(m) m$CI[1]),
                    upr = sapply(linear.confint, function(m) m$CI[2]))
  
  noutcome <- length(linear1)
  
  if (mapped) {
    mapped <- lapply(weights, 
      function(w){
        return(estimate_effect(data = data, weights = w, 
                                   matched = TRUE, 
                                   split.model = TRUE,
                                   doubly.robust = FALSE,
                                   estimand = "ATE",
                                   treatment.indicator = "z", 
                                   balance.covariates = balance.covariates,
                                   outcome = "y"))
      }
    )
    mapped.est <- sapply(mapped, function(h) h$estimate)
    mapped.confint <- lapply(mapped, confint)
    mapped.ci <- list(lwr = sapply(mapped.confint, function(m) m$CI[1]),
                      upr = sapply(mapped.confint, function(m) m$CI[2]))
    
  } else {
    mapped <- mapped.est <- rep(NA_real_, noutcome)
    mapped.ci <- list(lwr =  rep(NA_real_, noutcome), upr =  rep(NA_real_, noutcome))
  }
  
  
  
  method.labels <- rep(c("Naive", "Logistic",
                         "CBPS",
                         "SBW", "NNM",
                         "COT"), 2)
  
  
    bias   <- c(linear.est, mapped.est)
    lower.ci <- c(linear.ci$lwr, mapped.ci$lwr)
    upper.ci <- c(linear.ci$upr, mapped.ci$upr)
    in.ci <- lower.ci <=0 & upper.ci >= 0
    return( data.frame(method = method.labels,
                       bias = bias,
                       ci.lwr = lower.ci,
                       ci.upr = upper.ci,
                       in.ci  = in.ci)
    )
}


#### load data and clean ####
arraynum <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
iter <- as.numeric(Sys.getenv("ITER"))

jobid <- Sys.getenv('SLURM_ARRAY_JOB_ID')
data <- "crash"

# set seed
dataGen <- CRASH3$new()
dataGen$gen_data()


start <- proc.time()
output <- sample.eval(data = dataGen, 
                         power = 2,
                      mapped = TRUE)
end <- proc.time()
print(end - start)
date <- gsub(" ", "_", as.name(as.character(Sys.time())))
date <- gsub(":", "=", date)
term <- paste0(c(date, ".rds"), collapse="")
otfl <- paste0(c("crash_est",data,jobid,arraynum,iter, term),collapse="_")
otdr <- file.path("Output", data)
otfn <- file.path(otdr, otfl)
if(!dir.exists(otdr)) {
  drfl <- strsplit(otdr, "/")[[1]]
  for(fn in seq_along(drfl)) {
    curfile <- paste0(drfl[1:fn], collapse="/")
    if(dir.exists(curfile)) next
    dir.create(paste0(drfl[1:fn], collapse="/"))
  }
}

saveRDS(output, 
        file = otfn)

warnings()

q("no")

library(dplyr)
print(output %>% group_by(method, outcome, 
                          estimator, information, target.method,
                          ct) %>% summarise(mse = mean(bias^2),
                                            in.ci = mean(in.ci)), n = 50)
