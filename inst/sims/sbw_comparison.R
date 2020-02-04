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
outfile <- file.path("Output/sbw_comparison.rds")
if (file.exists(outfile)) {
  output <- readRDS(outfile)
} else {
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
  saveRDS(output, file = outfile)
}
#### print results ####

print(output$outcome %>% filter(estimate == "ATT") %>%
  group_by(doubly.robust, matched,
           standardized.mean.difference,
           weighting.method) %>%
  summarize(bias = mean(values),
            RMSE = sqrt(mean(values^2))),
  n = 40)

print(output$outcome %>% filter(estimate == "ATC") %>%
  group_by(doubly.robust, matched,
           standardized.mean.difference,
           weighting.method) %>%
  summarize(bias = mean(values),
            RMSE = sqrt(mean(values^2))),
n = 40)

print(output$outcome %>% filter(estimate == "ATE") %>%
  group_by(doubly.robust, matched,
           standardized.mean.difference,
           weighting.method) %>%
  summarize(bias = mean(values),
            RMSE = sqrt(mean(values^2))),
n = 40)
