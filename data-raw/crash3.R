#### clean and load CRASH-3 data ####

crash <- readxl::read_xlsx("../datasets/crash3/CRASH-3_dataset_anonymised_for_Freebird.xlsx")
crash$sex <- factor(crash$sex)
crash$timeSinceInjury <- difftime(crash$timeSinceInjury, 
                                  time2 = as.POSIXct("1899-12-31 00:00:00", tz = "UTC"),
                                  units = "hours")
crash$systolicBloodPressure <- ifelse(is.na(crash$systolicBloodPressure) & crash$sbpStatus == "Too low to be recorded" & !is.na(crash$sbpStatus),
                                      min(crash$systolicBloodPressure), crash$systolicBloodPressure)

crash$gcs <- as.numeric(sapply(strsplit(crash$gcsEyeOpening, " "), `[`, 1)) +
  as.numeric(sapply(strsplit(crash$gcsMotorResponse, " "), `[`, 1)) +
  as.numeric(sapply(strsplit(crash$gcsVerbalResponse, " "), `[`, 1))

crash$gcsEyeOpening     <- as.numeric(sapply(strsplit(crash$gcsEyeOpening, " "), `[`, 1))
crash$gcsMotorResponse  <- as.numeric(sapply(strsplit(crash$gcsMotorResponse, " "), `[`, 1))
crash$gcsVerbalResponse <- as.numeric(sapply(strsplit(crash$gcsVerbalResponse, " "), `[`, 1))
crash$gcsTiming <- factor(crash$gcsTiming)
crash$gcsTiming <- forcats::fct_recode(crash$gcsTiming,
                                       "after" = "After intubation/sedation",
                                       "before" = "Before intubation/sedation")
crash$pupilReact <- as.numeric(factor(crash$pupilReact, 
                                      ordered = TRUE, 
                                      labels = c("Unable to assess",
                                                 "None react",
                                                 "One reacts",
                                                 "Both react")))

crash$majorExtracranial <- factor(crash$majorExtracranial)
crash$intraCranialBleeding <- factor(crash$intraCranialBleeding)

for (i in match("epidural", colnames(crash)):match("eligible", colnames(crash))) {
  crash[, i] <- as.numeric(crash[, i] == "Yes")
  crash[crash$intraCranialBleeding == "No CT scan available" &
          !is.na(crash$intraCranialBleeding),i] <- 0
}



crash$time <- NA
crash$time <- difftime(crash$timerandtodeath, as.POSIXct("1899-12-31 00:00:00", tz = "UTC"),
                       units = "days")
crash$time[is.na(crash$timerandtodeath) & !is.na(crash$dischargeStatus)] <- 28

crash$death <- as.numeric(!is.na(crash$timerandtodeath))

subset.form <- formula( ~ siteId + sex + age + timeSinceInjury +
                          systolicBloodPressure + 
                          gcsEyeOpening +
                          gcsMotorResponse +
                          gcsVerbalResponse +
                          gcsTiming + 
                          intraCranialBleeding +
                          pupilReact +  majorExtracranial +
                          epidural + subdural + subarachnoid +
                          parenchymal + intraventricular + 
                          time + death + 0
)


csub <- model.frame(subset.form, data = crash)

analysis.form <- update(subset.form, death ~ z + .  - time - death - siteId -
                          intraCranialBleeding + 1)

csub$intraCranialBleeding[csub$intraCranialBleeding == "No CT scan available"] <- NA

sites <- expand.grid(unique(sort(csub$siteId)), 
                           unique(sort(csub$siteId)))
sites  <- sites[apply(sites,1, function(x) !any(duplicated(x))), ]

large.sites <- as.numeric(names(table(crash$siteId))[table(crash$siteId) > 100])

sites <- sites[sites$Var1 %in% large.sites & sites$Var2 %in% large.sites, ]

nrow(sites)

#remove groups without overlap problems
overlap.fun <- function(data, site1, site2) {
  
  subset.form <- formula(data)
  
  analysis.form <- update(subset.form, death ~ z + .  - time - death - siteId -
                            intraCranialBleeding + 1)
  
  
  data$z <- ifelse(data$siteId == site1, 1, 0)
  
  data$z[!data$siteId %in% c(site1, site2)] <- NA_real_
  z <- data[,"z"]
  
  # death.table <- with(data, table(z, death))
  # 
  # death.test <- diffpropci(x1 = death.table[2,2],
  #                          n1 = sum(z == 1),
  #                          x2 = death.table[1,2],
  #                          n2 = sum(z == 0),
  #                          conf.level = 0.95)
  # time.test  <- t.test(time ~ z, data = data)$conf.int
  
  noct <- cbind(y = data[is.na(data$intraCranialBleeding) & !is.na(z),"death"],
                model.matrix(update(analysis.form, . ~ .  - epidural
                                    - subdural - subarachnoid - parenchymal - intraventricular), 
                             data[is.na(data$intraCranialBleeding) & !is.na(z),])[,-1, drop = FALSE])
  noct.time <- data[is.na(data$intraCranialBleeding) & !is.na(z),"time"]
  
  
  ct <- cbind(y = data[!is.na(data$intraCranialBleeding) & !is.na(z),"death"],
              model.matrix(analysis.form, data[!is.na(data$intraCranialBleeding) & !is.na(z),])[,-1])
  ct.time <- data[!is.na(data$intraCranialBleeding) & !is.na(z),"time"]
  
  
  
  ct.conf <- colnames(ct)[-c(1:2)]
  noct.conf <- colnames(noct)[-c(1:2)]
  
  n0.ct <- sum(1 - ct[,"z"])
  n1.ct <- sum(ct[,"z"])
  
  n0.noct <- sum(1 - noct[,"z"])
  n1.noct <- sum(noct[,"z"])
  
  
  if (n0.ct > 2 && n1.ct > 2) {
    ct.est <- c(est1 = mean(ct[ct[,"z"] == 1, "y"], na.rm = TRUE) - mean(ct[ct[,"z"] == 0, "y"], na.rm = TRUE)
                , est2 = mean(ct.time[ct[,"z"] == 1], na.rm = TRUE) - mean(ct.time[ct[,"z"] == 0], na.rm = TRUE)
    )
    ct.ci.death <- qnorm(c(0.025, 0.975)) * sd(ct[ct[,"z"] == 1, "y"], na.rm = TRUE)/sqrt(n1.ct)
    ct.ci.time <- qnorm(c(0.025, 0.975)) * sd(ct.time[ct[,"z"] == 1], na.rm = TRUE)/sqrt(n1.ct)
    
    ct.inci <- c(ct.est[1] >= ct.ci.death[1] && ct.est[1] <= ct.ci.death[2],
                 ct.est[2] >= ct.ci.time[1] && ct.est[2] <= ct.ci.time[2])
    
    mb.ct <- mean(mean_bal(data = ct, balance.covariates = ct.conf,
                           outcome = "y", treatment.indicator = "z"))
    wb.ct <- wasserstein_p( a = rep(1/n0.ct, n0.ct),
                            b = rep(1/n1.ct, n1.ct),
                            cost = cost_fun(x = ct[, -c(1,2)],
                                            z = ct[, "z"],
                                            metric = "sdLp",
                                            estimand = "ATT",
                                            power = 1
                            ),
                            p = 1)
  } else {
    ct.est <- ct.inci <- rep(NA_real_, 2)
    mb.ct <- wb.ct <- NA_real_
  }
  
  if (n0.noct > 2 && n1.noct > 2) {
    noct.est <- c(est1 = mean(noct[noct[,"z"] == 1, "y"], na.rm = TRUE) - mean(noct[noct[,"z"] == 0, "y"], na.rm = TRUE)
                  , est2 = mean(noct.time[noct[,"z"] == 1], na.rm = TRUE) - mean(noct.time[noct[,"z"] == 0], na.rm = TRUE)
    )
    
    noct.ci.death <- qnorm(c(0.025, 0.975)) * sd(noct[noct[,"z"] == 1, "y"], na.rm = TRUE)/sqrt(n1.noct)
    noct.ci.time <- qnorm(c(0.025, 0.975)) * sd(noct.time[noct[,"z"] == 1], na.rm = TRUE)/sqrt(n1.noct)
    
    noct.inci <- c(noct.est[1] >= noct.ci.death[1] && noct.est[1] <= noct.ci.death[2],
                 noct.est[2] >= noct.ci.time[1] && noct.est[2] <= noct.ci.time[2])
    
    
    
    mb.noct <- mean(mean_bal(data = noct, balance.covariates = noct.conf,
                             outcome = "y", treatment.indicator = "z"))
    wb.noct <- wasserstein_p( a = rep(1/n0.noct, n0.noct),
                              b = rep(1/n1.noct, n1.noct),
                              cost = cost_fun(x = noct[, -c(1,2)],
                                              z = noct[, "z"],
                                              metric = "sdLp",
                                              estimand = "ATT",
                                              power = 1),
                              p = 1)
  } else {
    noct.est <- noct.inci <- rep(NA_real_, 2)
    mb.noct <- wb.noct <- NA_real_
  }
  
  return(data.frame(death.ct = ct.est[1],
                    death.noct = noct.est[1],
                    time.ct = ct.est[2],
                    time.noct = noct.est[2],
                    inci.death.ct = ct.inci[1],
                    inci.time.ct = ct.inci[2],
                    inci.death.noct = noct.inci[1],
                    inci.time.noct = noct.inci[2],
                    mean.bal.ct = mb.ct,
                    mean.bal.noct = mb.noct,
                    wass.bal.ct = wb.ct, 
                    wass.bal.noct = wb.noct,
                    site1 = site1,
                    site2 = site2
                   )
                    
  )
}

eval.list <- lapply(1:nrow(sites), function(i) overlap.fun(csub, 
                                              site1 = sites[i, 1],
                                              site2 = sites[i, 2]))
eval <- do.call("rbind", eval.list)


set.seed(824006901) #random.org
sites$arraynum <- rep(1:nrow(sites), length.out = nrow(sites))
sites$iter     <- 1#rep(1:3, each = 1000)[1:nrow(sites)]
sites$seed     <- c(sample.int(.Machine$integer.max, nrow(sites)))

data <- list(crash = csub,
             sitepairs = sites)

saveRDS(object = data, file = "../datasets/crash3/crash3.rds")
