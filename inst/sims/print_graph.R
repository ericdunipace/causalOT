library(causalOT)
library(dplyr)
library(xtable)

####
n           <- "512"
p           <- "6"
design      <- "B"
overlap     <- c("low", "high")
simsett     <- "hainmueller"
getdate     <- "2020-01-31"

#### Get Data ####
outname_low <- paste(c(simsett, "des", design, "overlap", overlap[1], "n", n, "date", getdate), collapse="_")
filepath_low<- file.path("Output", paste0(outname_low,".rds"))
data_low    <- readRDS(filepath_low)
outname_high<- paste(c(simsett, "des", design, "overlap", overlap[2], "n", n, "date", getdate), collapse="_")
filepath_high<- file.path("Output", paste0(outname_high,".rds"))
data_high   <- readRDS(filepath_high)

table_low   <- setup_data_for_table(data_low)
table_high  <- setup_data_for_table(data_high)

logis <- cbind(as.data.frame(table_high$unmatched$IPW), as.data.frame(table_low$unmatched$IPW[,c("Bias","RMSE")]))
logisM <- cbind(as.data.frame(table_high$matched$IPW), as.data.frame(table_low$matched$IPW[,c("Bias","RMSE")]))
other <- cbind(as.data.frame(table_high$unmatched$other), as.data.frame(table_low$unmatched$other[,c("Bias","RMSE")]))
otherM <- cbind(as.data.frame(table_high$unmatched$other), as.data.frame(table_low$unmatched$other[,c("Bias","RMSE")]))

#### graphs ####


#### print tables ####
print(data_high$outcome %>% 
        filter(estimate == "ATE" & wasserstein.power == 1 & metric == "Lp",
               standardized.mean.difference == 0.1) %>% 
        group_by(weighting.method
                 , doubly.robust
                 , matched 
                 # , standardized.mean.difference 
                 # , wasserstein.power
                 # , metric
        ) %>% 
        summarize(Bias = mean(values),
                  # Variance = var(values),
                  RMSE = sqrt(mean(values^2))),
      n = 71)

str(data$outcome %>% filter(estimate == "ATE" & 
                              weighting.method == "SBW" & 
                              standardized.mean.difference == 0.1 &
                              wasserstein.power == 2))

data_low$outcome %>% filter(estimate == "ATE" & 
                          weighting.method == "SBW" & 
                          standardized.mean.difference == 0.1 &
                          wasserstein.power == 2 & 
                          doubly.robust == TRUE & 
                          matched == FALSE &
                          metric == "Lp") %>% sapply(., unique)

data_low$outcome %>% filter(estimate == "ATE" & 
                          weighting.method == "Logistic" & 
                          standardized.mean.difference == 0.1 &
                          wasserstein.power == 2 & 
                          doubly.robust == TRUE & 
                          matched == FALSE &
                          metric == "Lp") %>% summarize(avg = mean(values), var = sqrt(mean(values^2)))

data_low$outcome %>% filter(estimate == "ATE" & 
                              weighting.method == "SBW" & 
                              standardized.mean.difference == 0.1 &
                              wasserstein.power == 2 & 
                              doubly.robust == TRUE & 
                              matched == FALSE &
                              metric == "Lp") %>% summarize(avg = mean(values), var = sqrt(mean(values^2)))

data_low$outcome %>% filter(estimate == "ATE" & 
                              weighting.method == "Wasserstein" & 
                              standardized.mean.difference == 0.1 &
                              wasserstein.power == 2 & 
                              doubly.robust == TRUE & 
                              matched == TRUE &
                              metric == "Lp") %>% summarize(avg = mean(values), var = sqrt(mean(values^2)))


data_low$outcome %>% 
  filter(estimate == "ATE" & metric == "Lp",
         standardized.mean.difference == 0.1) %>% 
  group_by(weighting.method
           , doubly.robust
           , matched 
           # , standardized.mean.difference 
           # , wasserstein.power
           # , metric
  ) %>% 
  summarize(Bias = mean(values[wasserstein.power == 1]),
            RMSE = sqrt(mean(values[wasserstein.power == 1]^2)),
            bias_2 = mean(values[wasserstein.power == 2]),
            RMSE_2 = sqrt(mean(values[wasserstein.power == 2]^2)))
