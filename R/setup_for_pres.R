setup_data_for_table <- function(data) {
  logistic <- data$outcome %>% 
    filter(estimate == "ATE" & metric == "mahalanobis" &
             weighting.method == "Logistic" & wasserstein.power == 2)
  
  notlog <- data$outcome %>% 
    filter(estimate == "ATE" & metric == "mahalanobis" &
             weighting.method != "Logistic" & wasserstein.power == 2)
  
  #### Gen tibbles ####
  
  ltibble <- logistic %>% 
    filter (matched == FALSE) %>%
    group_by(doubly.robust
             # , matched
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  
  colnames(ltibble)[1:2] <- c("Estimator", "$\\alpha$ truncation level")
  ltibble$Estimator <- ifelse(ltibble$Estimator, "DR Hajek", "Hajek")
  
  nl.tibbleDR <- notlog %>% 
    filter (matched == FALSE & weighting.method != "Wasserstein") %>%
    group_by(doubly.robust
             , weighting.method
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  
  nl.tibbleDRW <- notlog %>% 
    filter (matched == FALSE & weighting.method == "Wasserstein" & standardized.mean.difference== 0.1) %>%
    group_by(doubly.robust
             , weighting.method
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  nl.tibbleDRW$standardized.mean.difference[nl.tibbleDRW$weighting.method == "Wasserstein"] <- NA
  nl.tibbleDR <- rbind(nl.tibbleDR, nl.tibbleDRW)
  nl.tibbleDR <- nl.tibbleDR %>% arrange(desc(-doubly.robust))
  colnames(nl.tibbleDR)[1:3] <- c("Estimator", "Weights", "Std. Mean Difference")
  nl.tibbleDR$Estimator <- ifelse(nl.tibbleDR$Estimator, "DR Hajek", "Hajek")
  
  ltibbleM <- logistic %>% 
    filter (matched == TRUE) %>%
    group_by(doubly.robust
             # , matched
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  
  colnames(ltibbleM)[1:2] <- c("Estimator", "$\\alpha$ truncation level")
  ltibbleM$Estimator <- ifelse(ltibbleM$Estimator, "Bias Adj. Matched", "Matched" )
  
  nl.tibbleDRM <- notlog %>% 
    filter (matched == TRUE & weighting.method != "Wasserstein") %>%
    group_by(doubly.robust
             , weighting.method
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  
  nl.tibbleDRWM <- notlog %>% 
    filter (matched == TRUE & weighting.method == "Wasserstein" & standardized.mean.difference== 0.1) %>%
    group_by(doubly.robust
             , weighting.method
             , standardized.mean.difference
    ) %>%
    summarize(Bias = mean(values),
              # Variance = var(values),
              RMSE = sqrt(mean(values^2)))
  nl.tibbleDRWM$standardized.mean.difference[nl.tibbleDRWM$weighting.method == "Wasserstein"] <- NA
  nl.tibbleDRM <- rbind(nl.tibbleDRM, nl.tibbleDRWM)
  nl.tibbleDRM <- nl.tibbleDRM %>% arrange(desc(-doubly.robust))
  colnames(nl.tibbleDRM)[1:3] <- c("Estimator", "Weights", "Std. Mean Difference")
  nl.tibbleDRM$Estimator <- ifelse(nl.tibbleDRM$Estimator, "Bias Adj. Matched", "Matched")
  
  return(list(unmatched = list(IPW = ltibble, other = nl.tibbleDR),
         matched = list(IPW = ltibbleM, other = nl.tibbleDRM)))
}