set.seed(224893390) #from random.org

#### Load Packages ####
library(ROI)
library(ROI.plugin.ecos)
library(ROI.plugin.qpoases)
library(ROI.plugin.cplex)
library(ggplot2)
library(slam)
library(limbs)
library(causalOT)
library(doParallel)

#### Data Dim ####
n <- 2^9
p <- 6
nsims <- 100
power <- 2
std_mean_diff <- 0.1
overlap <- "low"
design <- "A"

#### get simulation functions ####
original <- Hainmueller$new(n = n, p = p, 
                      design = design, overlap = overlap)

#### dist mat ####
dist <- "Lp"
cost.calc <- switch(dist, "Lp" = causalOT::cost_calc_lp,
                    "mahalanobis" = causalOT::cost_mahalanobis)


#### Simulations ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
registerDoParallel(cl)
# registerDoSEQ()


low_overlap <- foreach(sim = 1:nsims) %dopar% {
  library(ROI)
  library(ROI.plugin.ecos)
  library(ROI.plugin.qpoases)
  library(ROI.plugin.cplex)
  library(slam)
  library(limbs)
  library(causalOT)
  library(doParallel)
  #library(transport)
  
  #### Gen Data ####
  target <- original$clone()
  target$gen_x()
  target$gen_z()
  target$gen_y()
  x <- target$get_x()
  z <- target$get_z()
  y <- target$get_y()
  
  x0 <- target$get_x0()#x[z==0,,drop=FALSE]
  x1 <- target$get_x1()#x[z==1,, drop = FALSE]
  y1 <- y[z==1]
  y0 <- y[z==0]
  
  ns <- target$get_n()
  n1 <- ns["n1"]
  n0 <- ns["n0"]
  
  #### setup holder data ####
  
  weights <- list(ATT = list(), 
                  ATC = list())
  mass <- list(original = list(mass0 = rep(1/n0, n0), 
                               mass1 = rep(1/n1,n1)),
               ATT = list(mass0 = list(sbw = list(), 
                                       Cwass = list(),
                                       wass=list()), 
                     mass1 = rep(1/n1,n1)),
               ATC = list(mass0 = rep(1/n0, n0)), 
                          mass1 = list(sbw = list(), 
                                       Cwass = list(),
                                       wass=list()))
  cost <- cost.calc(x0,x1, power)
  wass <- list("original" = NULL, 
               list(ATT = list(sbw = NULL, wass  = NULL)),
               list(ATC = list(sbw = NULL, wass  = NULL)),
               list(ATE = list(sbw = NULL, wass  = NULL)))
  W2 <- list("pre-match" = NULL, "ATT sbw" = NULL, "ATT constrained Wass" = NULL,
             "ATT Wass" = NULL,
             "ATC sbw" = NULL, "ATC constrained Wass" = NULL,
             "ATC Wass" = NULL,
             "ATE sbw" = NULL, "ATE Wass" = NULL, "ATE constrained Wass" = NULL,
             "feasible ATE sbw" = NULL, "feasible ATE Wass" = NULL, "feasible ATE constrained Wass" = NULL)
  est.list <- list(ATT = list(),
                   ATC = list(),
                   ATE = list(),
                   fATE = list()
  )
  dr.list <- list(H = est.list,
                  DRH = est.list)
  outcome <- list(sbw = dr.list,
                  Cwass = dr.list,
                  wass = dr.list
                  )
  
  #### original wass ####
  wass$original <- wasserstein_p(mass$original$mass0, mass$original$mass1, p = power,
                                          cost = cost)
  W2$`pre-match`<- wasserstein_p(mass$original$mass0, mass$original$mass1, p = 2,
                                          cost = cost)
  
  #### ATT ####
  #SBW
  weights$ATT$sbw <- calc_weight(target, std_mean_diff, method = "SBW", estimate = "ATT")
  wass$ATT$sbw <- wasserstein_p(a = weights$ATT$sbw, b = NULL, p = power, tplan = NULL, cost = cost)
  if(power == 2) {
    W2$`ATT sbw` <- wass$ATT$sbw
  } else {
    W2$`ATT sbw` <- wasserstein_p(a = weights$ATT$sbw, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  # Constrained Wasserstein
  weights$ATT$Cwass <- calc_weight(target, wass$ATT$sbw, method = "Constrained Wasserstein", estimate = "ATT")
  wass$ATT$Cwass <- wasserstein_p(a = weights$ATT$Cwass, b = NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`ATT constrained Wass` <-  wass$ATT$Cwass
  } else {
    W2$`ATT constrained Wass` <- wasserstein_p(a = weights$ATT$Cwass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  #Wasserstein
  weights$ATT$wass <- calc_weight(target, wass$ATT$sbw, method = "Wasserstein", estimate = "ATT")
  wass$ATT$wass <- wasserstein_p(a = weights$ATT$wass, b = NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`ATT Wass` <-  wass$ATT$wass
  } else {
    W2$`ATT Wass` <- wasserstein_p(a = weights$ATT$wass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  #### ATC ###
  #SBW
  weights$ATC$sbw <- calc_weight(target, std_mean_diff, method = "SBW", estimate = "ATC")
  wass$ATC$sbw <- wasserstein_p(weights$ATC$sbw, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`ATC sbw` <- wass$ATC$sbw
  } else {
    W2$`ATC sbw` <- wasserstein_p(a = weights$ATC$sbw, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  # Constrained Wasserstein
  weights$ATC$Cwass <- calc_weight(target, wass$ATC$sbw, method = "Constrained Wasserstein", estimate = "ATC")
  wass$ATC$Cwass <- wasserstein_p(weights$ATC$Cwass, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`ATC constrained Wass` <-  wass$ATC$Cwass
  } else {
    W2$`ATC constrained Wass` <- wasserstein_p(a = weights$ATC$Cwass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  # Wasserstein
  weights$ATC$wass <- calc_weight(target, wass$ATC$sbw, method = "Wasserstein", estimate = "ATC")
  wass$ATC$wass <- wasserstein_p(weights$ATC$wass, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`ATC Wass` <-  wass$ATC$wass
  } else {
    W2$`ATC Wass` <- wasserstein_p(a = weights$ATC$wass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  #### feasible ATE ####
  weights$fATE$sbw <- calc_weight(target, std_mean_diff, method = "SBW", estimate = "feasible")
  wass$fATE$sbw <- wasserstein_p(weights$fATE$sbw, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`feasible ATE sbw` <- wass$fATE$sbw
  } else {
    W2$`feasible ATE sbw` <- wasserstein_p(a = weights$fATE$sbw, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  # Constrained Wasserstein
  weights$fATE$Cwass <- calc_weight(target, wass$fATE$sbw, method = "Constrained Wasserstein", estimate = "feasible")
  wass$fATE$Cwass <- wasserstein_p(weights$fATE$Cwass, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`feasible ATE constrained Wass` <-  wass$ATC$Cwass
  } else {
    W2$`feasible ATE constrained Wass` <- wasserstein_p(a = weights$ATC$Cwass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  # Wasserstein
  weights$fATE$wass <- calc_weight(target, wass$fATE$sbw, method = "Wasserstein", estimate = "feasible")
  wass$fATE$wass <- wasserstein_p(weights$fATE$wass, b=NULL, p = power, tplan = NULL, cost = cost)
  
  if(power == 2) {
    W2$`feasible ATE Wass` <-  wass$ATC$wass
  } else {
    W2$`feasible ATE Wass` <- wasserstein_p(a = weights$ATC$wass, b = NULL, p = 2, tplan = NULL, cost = cost)
  }
  
  #### ATE ####
  weights$ATE$sbw <- convert_ATE(weights$ATC$sbw, weights$ATT$sbw)
  weights$ATE$Cwass <- convert_ATE(weights$ATC$Cwass, weights$ATT$Cwass)
  weights$ATE$wass <- convert_ATE(weights$ATC$wass, weights$ATT$wass)
  
  #### outcome calculation ####
  naive <- outcome_model(data = target, weights = list(w0=rep(1/n0,n0),
                                                       w1=rep(1/n1,n1)),
                         doubly.robust = FALSE, hajek = FALSE,
                         target = "ATE") #mean(y1) - mean(y0)
  # naive_att <- outcome_model(data = target, weights = list(w0=rep(1/n0,n0),
  #                                                          w1=rep(1/n1,n1)),
  #                            doubly.robust = FALSE, hajek = FALSE,
  #                            target = "ATT")#mean(y1) - mean(predict(lm(y0~x0), data.frame(x1)))
  # naive_atc <- outcome_model(data = target, weights = list(w0=rep(1/n0,n0),
  #                                                          w1=rep(1/n1,n1)),
  #                            doubly.robust = FALSE, hajek = FALSE,
  #                            target = "ATT")
  for(model in c("sbw", "Cwass", "wass")){
    for(dr in c(TRUE,FALSE)){
      dr.sel <- ifelse(dr, "DRH","H")
      for(est in c("ATT","ATC","ATE","fATE")){
        est.sel <- ifelse(est == "fATE","feasible",est)
        outcome[[model]][[dr.sel]][[est]] <- 
          outcome_model(data = target, weights = weights[[est]][[model]],
                        doubly.robust = dr, hajek = TRUE,
                        target = est.sel)
      }
        
    }
    
  }
  outcome$SBW$H$ATT <- outcome_model(data = target, weights = weights$ATT$sbw,
                                     doubly.robust = FALSE, hajek = TRUE,
                                     target = "ATT")
  outcome$SBW$H$ATC <- outcome_model(data = target, weights = weights$ATC$sbw,
                                     doubly.robust = FALSE, hajek = TRUE,
                                     target = "ATC")
  outcome$SBW$H$ATE <- outcome_model(data = target, weights = weights$ATE$sbw,
                                     doubly.robust = FALSE, hajek = TRUE,
                                     target = "ATE")
  outcome$SBW$H$fATE <- outcome_model(data = target, weights = weights$fATE$sbw,
                                     doubly.robust = FALSE, hajek = TRUE,
                                     target = "feasible")
  outcome$SBW$DRH$ATT <- outcome_model(data = target, weights = weights$ATT$sbw,
                                     doubly.robust = TRUE, hajek = TRUE,
                                     target = "ATT")
  outcome$SBW$DRH$ATC <- outcome_model(data = target, weights = weights$ATC$sbw,
                                     doubly.robust = TRUE, hajek = TRUE,
                                     target = "ATC")
  outcome$SBW$DRH$ATE <- outcome_model(data = target, weights = weights$ATE$sbw,
                                     doubly.robust = TRUE, hajek = TRUE,
                                     target = "ATE")
  outcome$SBW$DRH$fATE <- outcome_model(data = target, weights = weights$fATE$sbw,
                                       doubly.robust = TRUE, hajek = TRUE,
                                       target = "feasible")
  sbw_drh <- outcome_model_DRH(x,y,z, weights = list(z0 = weights, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  sbw_h <- outcome_model_H(x,y,z, weights = list(z0 = weights, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  
  w2_drh <- outcome_model_DRH(x,y,z, weights = list(z0 = w0_marg_w2, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  w2_h <- outcome_model_H(x,y,z, weights = list(z0 = w0_marg_w2, z1 = rep(1/sum(z), sum(z))), est = "ATE")
  
  
  sbw_drh_ate <- outcome_model_DRH(x,y,z, weights = list(z0 = weights, z1 = weights_ate), est = "ATE")
  sbw_h_ate <- outcome_model_H(x,y,z, weights = list(z0 = weights, z1 = weights_ate), est = "ATE")
  
  w2_drh_ate <- outcome_model_DRH(x,y,z, weights = list(z0 = w0_marg_w2, z1 = w1_marg_w2_ate), est = "ATE")
  w2_h_ate <- outcome_model_H(x,y,z, weights = list(z0 = w0_marg_w2, z1 = w1_marg_w2_ate), est = "ATE")
  
  # print(c("Naive" = naive,
  #         "SBW Hajek" = sbw_h_ate,
  #         "W2 Hajek" =w2_h_ate,
  #         "SBW DR Hajek" = sbw_drh_ate,
  #         "W2 DR Hajek" = w2_drh_ate))
  out <- list(ATT = list("Naive" = naive_att,
                          "SBW Hajek" = sbw_h,
                          "W2 Hajek" =w2_h,
                          "SBW DR Hajek" = sbw_drh,
                          "W2 DR Hajek" = w2_drh),
              ATC = list("Naive" = naive_att,
                         "SBW Hajek" = sbw_h,
                         "W2 Hajek" =w2_h,
                         "SBW DR Hajek" = sbw_drh,
                         "W2 DR Hajek" = w2_drh),
              ATE = list("Naive" = naive,
                              "SBW Hajek" = sbw_h_ate,
                              "W2 Hajek" =w2_h_ate,
                              "SBW DR Hajek" = sbw_drh_ate,
                              "W2 DR Hajek" = w2_drh_ate),
              `feasible ATE` = list("Naive" = naive_att,
                   "SBW Hajek" = sbw_h,
                   "W2 Hajek" =w2_h,
                   "SBW DR Hajek" = sbw_drh,
                   "W2 DR Hajek" = w2_drh),
              W2 = W2,
              ESS = list(control = list(pre = sum(z==0), SBW = 1/sum(weights^2), Wass = 1/sum(w0_marg_w2^2)),
                         treated = list(pre = sum(z==1), SBW = 1/sum(weights_ate^2), Wass = 1/sum(w1_marg_w2_ate^2))
              )
              )
  return(out) 
}
stopCluster(cl)

#### Calculate summary stat ####

ATT <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ATT)))
ATE <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ATE)))
W2 <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$W2)))
ESS <- list()
ESS$control <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ESS$control)))
ESS$treated <- do.call("rbind", lapply(low_overlap, function(l) as.data.frame(l$ESS$treated)))

colMeans(ATT)
colMeans(ATE)
colMeans(W2)
colMeans(ESS)

colVar(ATT)#*(nsims-1)/nsims
colVar(ATE)#*(nsims-1)/nsims
colVar(W2)
colVar(ESS)

colMeans(ATT^2)
colMeans(ATE^2)


# saveRDS(low_overlap, file = "low_w2_2020019.rds")
