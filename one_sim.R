set.seed(224893390) #from random.org

#### Load Packages ####
library(causalOT)
library(dplyr)

#### Sim param ####
n <- 2^9
p <- 6
nsims <- 100
overlap <- "low"
design <- "B"
distance <- c("Lp", "mahalanobis")
power <- c(1,2)
ground_power <- 2
std_mean_diff <- c(0.001, 0.01, 0.1)

#### get simulation functions ####
original <- Hainmueller$new(n = n, p = p, 
                            design = design, overlap = overlap)

#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
# cl <- FALSE
times <- proc.time()
output <- sim.function(original, nsims, ground_p = ground_power, p = power, 
                       standardized.mean.difference = std_mean_diff,
                       distance = distance, parallel = cl)
parallel::stopCluster(cl)
print(proc.time() - times)
saveRDS(output, file = "low_B_20200126.rds")

#### Calculate summary stat ####
print(output$outcome %>% 
  filter(estimate == "ATE") %>% 
  group_by(weighting.method, doubly.robust, matched, 
           standardized.mean.difference, wasserstein.power, metric) %>% 
  summarize(bias = mean(values),
            variance = var(values),
            mse = mean(values^2)),
  n = 39)


output$`ESS/N` %>% 
  group_by(weighting.method, doubly.robust, matched, 
           standardized.mean.difference, wasserstein.power) %>% 
  summarize(average = mean(values),
            variance = var(values))

# output$`2-Wasserstein` %>% 
#   group_by(weighting.method, doubly.robust, matched, 
#            standardized.mean.difference, wasserstein.power) %>% 
#   summarize(average = mean(values),
#             variance = var(values))
# 
# output[[paste(c(power,"Wasserstein"),collapse="-")]] %>% 
#   group_by(estimate, weighting.method) %>% 
#   summarize(average = mean(values),
#             variance = var(values))

#check number sims correct
output$outcome %>% 
  filter(weighting.method == "None") %>%
  summarize(total_sims = n() == nsims)


