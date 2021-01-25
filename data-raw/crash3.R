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

set.seed(824006901) #random.org
sites$arraynum <- rep(1:nrow(sites), length.out = nrow(sites))
sites$iter     <- 1#rep(1:3, each = 1000)[1:nrow(sites)]
sites$seed     <- c(sample.int(.Machine$integer.max, nrow(sites)))

data <- list(crash = csub,
             sitepairs = sites)

saveRDS(object = data, file = "../datasets/crash3/crash3.rds")
