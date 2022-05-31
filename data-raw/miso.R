library(causalOT)
library(forcats)
library(dplyr)
library(foreign)

#### misoprostol ####
set.seed(972569665) #from random.org
# data from https://doi.org/10.7910/DVN/ETHH4N
miso <- foreign::read.spss(file = "../datasets/misoprostol/final_pphdatabase_merge_topost2.sav",
                           to.data.frame = TRUE)
mcn  <- colnames(miso)
start.col   <- match("caseID", mcn)
end.col     <- match("q4_3", mcn)
outcome.col <- match("addbldloss", mcn)
keep.col    <- c(start.col:end.col, outcome.col)
old.cn      <- factor(mcn[keep.col])
analy.vars  <- c("patID", "sitecode", "tx", "age","age_est",
                 "no_educ", "cur_married", "woman_prof",
                 "husband_prof",
                 "woman_educ", "married","num_livebirth",
                 "woman_occup", "woman_occup_2", 
                 "husband_occup", "husband_occup_2",
                 "hb_test",
                 "gest_age",
                 "prev_pph", "hosp_deliv", 
                 "datetime_delivery",
                 "singleton","epidural",
                 "induced_labor", "augmented_labor", "oxytocin_3rd_stage",
                 "other_uterotonics", "early_cordclamp", "control_cordtraction",
                 "uterine_massage", "datetime_placenta",
                 "bloodlossattx",
                 "datetime_pph_tx", "pph_dx_method",
                 "cum_blood_20m"
)
analysis.form.B <- as.formula(cum_blood_20m ~ tx + age + 
                                no_educ +
                                num_livebirth +
                                cur_married +
                                gest_age + prev_pph + hb_test + 
                                induced_labor +
                                augmented_labor +
                                oxytocin_3rd_stage + early_cordclamp +
                                control_cordtraction + uterine_massage + 
                                placenta + bloodlossattx 
)

analysis.form.A <- update(analysis.form.B, ~ . - induced_labor - augmented_labor -
                            oxytocin_3rd_stage)

# pphloss blood loss at pph dx
# q1_3 time of pph tx
# q2_3 blood loss at pph tx
# q4_3 cum blood loss 20 min after tx
new.cn      <- forcats::fct_recode(old.cn,
                                   "patID" = "ptID_1",
                                   "tx" = "txtarm",
                                   "staff_initials" = "staff_1",
                                   "date" = "date_1",
                                   "csection_planned" = "q1_1",
                                   "prosta_allergy" = "q2_1",
                                   "conditions_inhibit" = "q3_1",
                                   "specific_cond_inhibit" = "q3a_1",
                                   "elligible" = "q4_1",
                                   "age" = "q5_1",
                                   "age_est" = "q5a_1",
                                   "woman_educ" = "q6_1",
                                   "woman_educ_other" = "q6a_1",
                                   "married" = "q7_1",
                                   "num_livebirth" = "q8_1",
                                   "num_livechild" = "q9_1",
                                   "woman_occup" = "q10_1",
                                   "woman_occup_2" = "q10a_1",
                                   "husband_occup" = "q11_1",
                                   "husband_occup_2" = "q11a_1",
                                   "gest_age" = "q13_1",
                                   "prev_pph" = "q14_1",
                                   "hb_test" = "q17_1",
                                   "hosp_deliv" = "q18_1",
                                   "datetime_delivery" = "q2_2",
                                   "singleton" = "q3_2",
                                   "epidural" = "q4_2",
                                   "induced_labor" = "q5_2",
                                   "augmented_labor" = "q5a_2",
                                   "oxytocin_3rd_stage" = "q6_2",
                                   "other_uterotonics" = "q7_2",
                                   "early_cordclamp" = "q8_2",
                                   "control_cordtraction" = "q9_2",
                                   "uterine_massage" = "q10_2",
                                   "datetime_placenta" = "q11_2",
                                   "datetime_pph_tx" = "q1_3",
                                   "pph_dx_method" = "q14_2",
                                   "cum_blood_20m" = "q4_3",
                                   "bloodlossattx" = "q2_3"
)

# miso <- miso[, keep.col]
colnames(miso) <- as.character(new.cn)
miso$no_educ <- ifelse(miso$woman_educ == "no education", 1, 0)
miso$cur_married <- ifelse(miso$married == "married/cohabiting", 1, 0)
miso$woman_prof <- ifelse(miso$woman_occup == "unemployed", 1, 0)
miso$husband_prof <- ifelse(miso$husband_occup == "unemployed" |
                              miso$husband_occup == "not applicable", 1, 0)

mB <- miso[miso$study == "study B, oxy prophylaxis", analy.vars] #prophylaxis

# clean data
mB$age[is.na(mB$age)] <- mB$age_est[is.na(mB$age)]
mB <- mB[!is.na(mB$woman_occup), ]

mB$woman_occup_2 <- factor(mB$woman_occup_2)
mB$datetime_delivery <- as.POSIXct(mB$datetime_delivery, 
                                   origin = "1582-10-14", tz = "GMT")
mB$datetime_placenta <- as.POSIXct(mB$datetime_placenta, 
                                   origin = "1582-10-14", tz = "GMT")
mB$datetime_pph_tx   <- as.POSIXct(mB$datetime_pph_tx, 
                                   origin = "1582-10-14", tz = "GMT")
mB$placenta <- ifelse(mB$datetime_pph_tx >= mB$datetime_placenta, 1, 0)
mB$placenta_time <- difftime(mB$datetime_placenta, mB$datetime_delivery,
                             units = "mins")
mB$tx  <- ifelse(mB$tx == "800mcg miso SL", 1, 0 )

#### All sites
mBf <- model.frame(analysis.form.B, 
                   data = mB)
pph  <- cbind(cum_blood_20m = model.response(mBf), 
             model.matrix(terms(mBf), mBf)[,-1])
pph_sitecode <- mB %>% 
  mutate(sitecode = droplevels(sitecode)) %>% 
  filter(rownames(mB) %in% rownames(pph)) %>% 
  .$sitecode
pph <- cbind(pph[, matrixStats::colSds(mb) != 0], sitecode = as.numeric(pph_sitecode))

usethis::use_data(pph, overwrite = TRUE, internal = FALSE)

#### Egypt vs others in B ####
mBf <- model.frame(analysis.form.B, data = mB)
mb  <- cbind(cum_blood_20m = model.response(mBf), 
             model.matrix(terms(mBf), mBf)[,-1])


mBf %>% 
   mutate(sitecode = mb_sitecode) %>% 
   group_by(sitecode) %>% 
   summarise(stat = t.test(x = cum_blood_20m[tx == 1],
                           y = cum_blood_20m[tx == 0])$statistic,
             mean.y1 = mean(cum_blood_20m[tx == 1]),
             mean.y0 = mean(cum_blood_20m[tx == 0]),
             tx.effect = mean(cum_blood_20m[tx == 1]) - mean(cum_blood_20m[tx == 0]),
             CI.lwr = t.test(x = cum_blood_20m[tx == 1],
                             y = cum_blood_20m[tx == 0])$conf.int[1],
             CI.upr = t.test(x = cum_blood_20m[tx == 1],
                             y = cum_blood_20m[tx == 0])$conf.int[2],
             n0 = sum(tx == 0),
             n1 = sum(tx == 1))

mb_egypt_1 <- mb[mb_sitecode == "Cairo, Egypt" &
                    mb[,"tx"] == 1 , -match(c("tx"), colnames(mb))]
mb_egypt_0 <- mb[mb_sitecode == "Cairo, Egypt" &
                    mb[,"tx"] == 0, -match(c("tx"), colnames(mb))]

mb_other_1 <- mb[mb_sitecode != "Cairo, Egypt" &
                    mb[,"tx"] == 1 , -match(c("tx"), colnames(mb))]
mb_other_0 <- mb[mb_sitecode != "Cairo, Egypt" &
                    mb[,"tx"] == 0, -match(c("tx"), colnames(mb))]

mb_att <- cbind(rbind(cbind(tx = 1, mb_egypt_1), cbind(tx = 0, mb_other_0)))

mb_att <- mb_att[, matrixStats::colSds(mb_att) != 0]

saveRDS(object = mb_att, file = "../datasets/misoprostol/miso_egypt.rds")
pph_egypt <- mb_att
usethis::use_data(pph_egypt, overwrite = TRUE, internal = FALSE)


