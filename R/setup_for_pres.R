setup_data_for_table <- function(high, low, metrics = "mahalanobis", 
                                 estimands = "ATE",
                                 p = 2) {
  
  high_out <- setup_data_work(high, metrics = metrics,
                  estimands = estimands, p = p)
  
  low_out <- setup_data_work(low, metrics = metrics,
                             estimands = estimands, p = p)
  
  
  
  return(list(unmatched = grp_tibble %>% filter(Matched == FALSE),
              matched = grp_tibble %>% filter(Matched == TRUE)))
}

setup_data_work <- function(data, metrics = "mahalanobis", 
                                 estimands = "ATE",
                                 p = 2) {
  
  drop.single <- function(dat, col) {
    xx <- dat[[col]]
    if(length(unique(xx)) == 2) {
      if(any(unique(xx) == "(Missing)")) {
        return(NULL)
      }
    }
    return(xx)
  }
  
  data$outcome$delta <- forcats::fct_explicit_na(data$outcome$delta)
  data$outcome$metric <- forcats::fct_explicit_na(data$outcome$metric)
  data$outcome$method <- factor(data$outcome$method, levels = unique(data$outcome$method))
  data$outcome$estimand <- factor(data$outcome$estimand, levels = unique(data$outcome$estimand))
  
  logistic <- data$outcome %>% 
    filter(estimand %in% estimands & method == "Logistic" )
  
  wass <- data$outcome %>% 
    filter(estimand %in% estimands & metric %in% metrics &
             grepl("Wasserstein", method) & wass_p == p)
  other <- data$outcome %>% 
    filter(estimand %in% estimands & 
             metric %in% c("(Missing)", metrics[1]) &
             (method %in% c("SBW","RKHS","RKHS.dose") ))
  
  combine <- rbind(logistic, other, wass)
  
  #### Gen tibbles ####
  grp_tibble <- combine %>% 
    # filter (match %in% matched) %>%
    group_by(estimand
             , model.augmentation
             , match
             , method
             , delta
             , metric
    ) %>%
    summarize(Bias = mean(estimate)*100,
              # Variance = var(values),
              RMSE = sqrt(mean(estimate^2))*100)
  
  new.lab <- c("Estimand", "Estimator", "Matched", "Method", "$\\alpha$ truncation level", "Metric")
  colnames(grp_tibble)[1:length(new.lab)] <- new.lab
  grp_tibble$Estimator <- ifelse(grp_tibble$Estimator, "DR Hajek", "Hajek")
  
  
  grp_tibble$Metric <- drop.single(grp_tibble, "Metric")
  grp_tibble$"$\\alpha$ truncation level" <- drop.single(grp_tibble, "$\\alpha$ truncation level")
  
  
  return(list(unmatched = grp_tibble %>% filter(Matched == FALSE),
         matched = grp_tibble %>% filter(Matched == TRUE)))
}