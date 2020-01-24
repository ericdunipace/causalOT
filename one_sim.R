set.seed(224893390) #from random.org

#### Load Packages ####
library(causalOT)
library(dplyr)

#### Data Fun ####
n <- 2^9
p <- 6
nsims <- 100
power <- 1
ground_power <- 1
std_mean_diff <- 0.1
overlap <- "low"
design <- "A"
distance <- "Lp"

#### get simulation functions ####
original <- Hainmueller$new(n = n, p = p, 
                            design = design, overlap = overlap)



#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
registerDoParallel(cl)
output <- sim.function(original, nsims, ground_p = ground_power, p = power, 
                       std_mean_diff = std_mean_diff,
                       distance = distance, parallel = cl)
stopCluster(cl)

#### Calculate summary stat ####
output$outcome %>% 
  group_by(estimate, weighting.method, doubly.robust) %>% 
  summarize(bias = mean(values),
            variance = var(values),
            mse = mean(values^2))


output$`ESS/N` %>% 
  group_by(estimate, weighting.method, Population) %>% 
  summarize(average = mean(values),
            variance = var(values))

output$`2-Wasserstein` %>% 
  group_by(estimate, weighting.method) %>% 
  summarize(average = mean(values),
            variance = var(values))

output[[paste(c(power,"Wasserstein"),collapse="-")]] %>% 
  group_by(estimate, weighting.method) %>% 
  summarize(average = mean(values),
            variance = var(values))

#check number sims correct
output$outcome %>% 
  filter(weighting.method == "None") %>%
  summarize(total_sims = n() == nsims)

# saveRDS(output, file = "low_w1_20200123.rds")
