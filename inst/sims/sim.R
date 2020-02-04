set.seed(-110942628) #from random.org

#### Load Packages ####
library(causalOT)
library(dplyr)

#### Sim Param ####
n <- 2^9
p <- 6
nsims <- 3
power <- c(1,2)
ground_power <- c(1,2)
std_mean_diff <- c(0.001, 0.01, 0.1, 0.2)
overlap <- c("low","high")
design <- c("A","B")
distance <- c("mahalanobis")

#### get simulation functions ####
original <- seeds <- data_sim_holder(design, overlap)

for (o in overlap) {
  for (d in design) {
    original[[o]][[d]] <- Hainmueller$new(n = n, p = p,
                                          design = d, overlap = o)
    seeds[[o]][[d]] <- sample.int(.Machine$integer.max, size = 1)
  }
}

#### output holder function ####
output <- sim_holder(design, overlap)

#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
times <- proc.time()
for (o in overlap) {
  for (d in design) {
    output[[o]][[d]] <- 
            sim.function(original[[o]][[d]], 
             nsims, 
             ground_p = ground_power, 
             p = power, 
             standardized.mean.difference = std_mean_diff,
             distance = distance, parallel = cl,
             seed = seeds[[o]][[d]])
  }
}
run.time <- (proc.time() - times)
print(run.time)
parallel::stopCluster(cl)
saveRDS(output, file = "all_20200126.rds")

#### Calculate summary stat ####
o <- "low"
smd <- 0.01
cat(o, "overlap, design",d, ", using", pp, "-Wasserstein distance, and standardized mean difference of",smd )
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

if(pp!=2) {
  to_print[[paste(c(pp,"Wasserstein"),collapse="-")]] %>% 
  group_by(estimate, weighting.method) %>% 
  summarize(average = mean(values),
            variance = var(values))
}

#check number sims correct
to_print$outcome %>% 
  filter(weighting.method == "None") %>%
  summarize(total_sims = n() == nsims)

