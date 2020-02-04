set.seed(092834324) #from random.org

#### Load Packages ####
library(causalOT)
library(dplyr)
library(sbw)

#### Sim param ####
n <- 2^9
p <- 6
nsims <- 1000
overlap <- "low"
design <- "B"
std_mean_diff <- c(0.01, 0.05, 0.1, 0.2)
solver <- "cplex" # "gurobi"

#### get simulation functions ####
original <- Hainmueller$new(n = n, p = p, 
                            design = design, overlap = overlap)

#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
# cl <- FALSE
times <- proc.time()
# debugonce(compare_sbw_mine)
output <- compare_sbw_mine(original, nsims = nsims, 
                              standardized.mean.difference = std_mean_diff,
                              solver = solver,
                              parallel = cl)
run.time <- (proc.time() - times)
print(run.time)
parallel::stopCluster(cl)

#### print results ####

output$outcome %>% filter(estimate == "ATE") %>%
  group_by(weighting.method,doubly.robust, matched,
           standardized.mean.difference) %>%
  summarize(bias = mean(values),
            RMSE = sqrt(mean(values^2)))
