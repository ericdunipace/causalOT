set.seed(-110942628) #from random.org

#### Load Packages ####
library(causalOT)
library(dplyr)

#### Data Fun ####
n <- 2^9
p <- 6
nsims <- 3
power <- c(1,2)
ground_power <- 2
std_mean_diff <- c(0.001, 0.01, 0.1, 1)
overlap <- c("low","high")
design <- c("A","B")
distance <- c("Lp", "mahalanobis")

#### get simulation functions ####
original <- data_sim_holder(design, overlap)
  
for (o in overlap) {
  for (d in design) {
    original[[o]][[d]] <- Hainmueller$new(n = n, p = p, 
                                          design = d, overlap = o)
  }
}

#### output holder function ####
output <- sim_holder(design, overlap, power, distance, std_mean_diff)

#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
for (o in overlap) {
  for (d in design) {
    for(pp in power) {
      for (dist in distance) {
        for(smd in std_mean_diff) {
          output[[o]][[d]][[pp]][[dist]][[as.character(smd)]] <- 
            sim.function(original[[o]][[d]], nsims, 
             ground_p = ground_power, 
             p = pp, 
             std_mean_diff = smd,
             distance = dist, parallel = cl)
        }
      }
    }
  }
}
parallel::stopCluster(cl)
saveRDS(output, file = "all_20200124.rds")

#### Calculate summary stat ####
to_print <- output[[o]][[d]][[pp]][[dist]][[smd]]
to_print$outcome %>% 
  group_by(estimate, weighting.method, doubly.robust) %>% 
  summarize(bias = mean(values),
            variance = var(values),
            mse = mean(values^2))


to_print$`ESS/N` %>% 
  group_by(estimate, weighting.method, Population) %>% 
  summarize(average = mean(values),
            variance = var(values))

to_print$`2-Wasserstein` %>% 
  group_by(estimate, weighting.method) %>% 
  summarize(average = mean(values),
            variance = var(values))

to_print[[paste(c(pp,"Wasserstein"),collapse="-")]] %>% 
  group_by(estimate, weighting.method) %>% 
  summarize(average = mean(values),
            variance = var(values))

#check number sims correct
to_print$outcome %>% 
  filter(weighting.method == "None") %>%
  summarize(total_sims = n() == nsims)

