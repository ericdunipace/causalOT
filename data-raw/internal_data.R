#### setup lalonde data ####
devtools::install_github("jjchern/lalonde")

library(lalonde)
library(dplyr)

lalonde_nsw <- lalonde::nsw_dw
# usethis::use_data(lalonde_dat, overwrite = TRUE, internal = TRUE)

lalonde_full <- lalonde_nsw %>% filter(treat == 1) %>% 
  rbind(lalonde::cps_controls)


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

subset.form <- formula( death ~ siteId + sex + age + timeSinceInjury +
                          systolicBloodPressure + 
                          gcsEyeOpening +
                          gcsMotorResponse +
                          gcsVerbalResponse +
                          gcsTiming + 
                          # intraCranialBleeding + # everyone has intracranial bleeding
                          pupilReact +
                          epidural + subdural + subarachnoid +
                          parenchymal + intraventricular + 
                          1
)


crash3 <-  model.frame(subset.form, data = crash, subset = crash$eligible==1)



#### Save all internal data ####
usethis::use_data(lalonde_full, lalonde_nsw, crash3, overwrite = TRUE, internal = TRUE)
