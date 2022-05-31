#### Load Packages ####
library(causalOT)
library(forcats)
library(dplyr)
library(ggplot2)
library(xtable)
library(ggsci)
# library(rlang)

#### Load data ####
data("pph", package = "causalOT")
reticulate::use_python("/usr/bin/python3", required = TRUE)

#### setup treatment indicators ####
bal.cov <- colnames(pph)[-c(1:2)]
outcome <- "cum_blood_20m"
tx.ind <- "tx"
continuous.covar <- c("age","gest_age","num_livebirth","hb_test","bloodlossattx")
binary.covar <- colnames(pph)[sapply(colnames(pph), function(x) 2 == length(unique(pph[,x])))][-1]

#### Estimation function ####
miso_estimate <- function(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site) {
  
  weight.by.tx <- function(data, continuous.covar, binary.covar, outcome, tx.ind, estimand = "ATT") {
    bal.cov <- c(continuous.covar, binary.covar)
    power <- 2
    z <- data[,tx.ind]
    n0 <- sum(1 - z)
    n1 <- sum(z)
    
    cost <- cost_fun( x = cbind(scale(data[,continuous.covar], scale = FALSE)
                                %*% causalOT:::inv_sqrt_mat(cov(data[,continuous.covar]))
                                ,  data[,binary.covar]), z = z, metric = "Lp", power = power,
                      estimand = estimand)
    scaled_data <- cbind(scale(data[,continuous.covar], scale = FALSE)
                        %*% causalOT:::inv_sqrt_mat(cov(data[,continuous.covar]))
                        ,  data[,binary.covar], tx = z)
    colnames(scaled_data) <- c(continuous.covar, binary.covar,tx.ind)
    
    weights <- list()
    weights$GLM  = calc_weight(data = data, 
                               constraint = 1e-4,
                               estimand = estimand, 
                               method = "Logistic",
                               formula = "z ~ . + .*. + I(.^2)",
                               treatment.indicator = tx.ind, 
                               balance.covariates = bal.cov,
                               outcome = outcome)
    
    weights$CBPS  = calc_weight(data = data, 
                                estimand = estimand, 
                                method = "CBPS",
                                formula = "z~. + .*. + I(.^2) + 0",
                                treatment.indicator = tx.ind, 
                                balance.covariates = bal.cov,
                                outcome = outcome)
    
    
    weights$SBW  = calc_weight(data = data, 
                               estimand = estimand, 
                               method = "SBW",
                               grid.search = TRUE,
                               solver = "mosek",
                               grid.length = 20,
                               grid = seq(0, 2, length.out = 20),
                               formula = "~. + .*. + I(.^2) + 0",
                               treatment.indicator = tx.ind, 
                               balance.covariates = bal.cov,
                               outcome = outcome)
    
    
    weights$SCM  = calc_weight(data = data,
                               estimand = estimand,
                               method = "SCM",
                               penalty = "none",
                               add.divergence = FALSE,
                               solver = "mosek",
                               treatment.indicator = tx.ind,
                               balance.covariates = bal.cov,
                               outcome = outcome)
    
    weights$NNM  = calc_weight(data = scaled_data,
                               estimand = estimand,
                               method = "NNM",
                               p = ncol(scaled_data)/2+1,
                               treatment.indicator = tx.ind,
                               balance.covariates = bal.cov,
                               outcome = outcome)
    
    
    weights$COT =  calc_weight(data = scaled_data, 
                               estimand = estimand, 
                               method = "Wasserstein",
                               solver = "mosek",
                               grid.search = TRUE,
                               p = power, grid.length = 8,
                               metric = "Lp",
                               add.divergence = TRUE,
                               penalty = "entropy",
                               treatment.indicator = tx.ind, 
                               balance.covariates = bal.cov,
                               wass.method = "sinkhorn_geom",
                               n.boot = 1e2, niter = 1e3,
                               tol = 1e-10, stepsize = 1e-2,
                               verbose = FALSE)
    
    weights$COT$gamma <- weights$NNM$gamma <- weights$SCM$gamma <- NULL
    epsilon <- 1/log(sum(abs(data[,bal.cov]))*nrow(data)^(1/length(bal.cov))) / median(cost)
    for (i in 1:length(weights)) {
      weights[[i]]$gamma <- causalOT:::calc_gamma(weights[[i]], cost = cost, p = 1, 
                                                  niter = 1e4, epsilon = epsilon)
    }
    
    return(weights)
    
  }
  estimate.by.tx <- function(data, weights, continuous.covar, binary.covar, outcome, tx.ind, estimand = "ATT", orig.tx, orig.ci) {
    bal.cov <- c(continuous.covar, binary.covar)
    power <- 2
    z <- data[,tx.ind]
    n0 <- sum(1 - z)
    n1 <- sum(z)
    nsave <- if(estimand == "ATT") {
      n1
    } else {
      n0
    }
    
    cost <- cost_fun( x = cbind(scale(data[,continuous.covar], scale = FALSE)
                                %*% causalOT:::inv_sqrt_mat(cov(data[,continuous.covar]))
                                ,  data[,binary.covar]), z = z, metric = "Lp", power = power,
                      estimand = estimand)
    scaled_data <- cbind(scale(data[,continuous.covar], scale = FALSE)
                         %*% causalOT:::inv_sqrt_mat(cov(data[,continuous.covar]))
                         ,  data[,binary.covar], tx = z)
    colnames(scaled_data) <- c(continuous.covar, binary.covar,tx.ind)
    
    hajek <- lapply(weights,
                            function(w) estimate_effect(data = data, weights = w, 
                                                        matched = FALSE, 
                                                        split.model = TRUE,
                                                        doubly.robust = FALSE,
                                                        estimand = estimand,
                                                        treatment.indicator = tx.ind, 
                                                        balance.covariates = bal.cov,
                                                        outcome = outcome))
    bp <- lapply(weights,
                         function(w) {
                           temp <- w
                           temp$args$power = 1
                           estimate_effect(data = data, weights = temp, 
                                                     matched = TRUE, 
                                                     split.model = TRUE,
                                                     doubly.robust = FALSE,
                                                     estimand = estimand,
                                                     p = 1,
                                                     cost = cost,
                                                     treatment.indicator = tx.ind, 
                                                     balance.covariates = bal.cov,
                                                     outcome = outcome)})
    
    ci.hajek <-  lapply(hajek, confint, verbose = TRUE, method = "asymptotic",
                                    model = "lm", formula = list(treated = "y ~ .",
                                                                 control = "y ~ ."))
    
    ci.bp <- lapply(bp, confint, verbose = TRUE, method = "asymptotic",
                            model = "lm", formula = list(treated = "y ~ .",
                                                         control = "y ~ ."))
    
    init.wass <- sinkhorn(x = scaled_data[z==1,], y = scaled_data[z==0,],
               a = rep(1/n1,n1), b = rep(1/n0,n0),
               power = 2, blur = max(cost)^2, metric = "Lp")$loss
  
    final.wass <- sapply(weights, function(w) {
      sinkhorn(x = scaled_data[z==1,], y = scaled_data[z==0,],
                          a = w$w1, b = w$w0,
               power = 2, blur = max(cost)^2, metric = "Lp")$loss
    })
    
    return(rbind(data.frame(estimator = "hajek",
                            method = names(weights),
                            estimate = sapply(hajek, function(e) e$estimate),
                             sd = sapply(ci.hajek, function(e) e$SD),
                             cover.orig = sapply(ci.hajek, function(e) e$CI[1] < orig.tx && e$CI[2] > orig.tx),
                             est.in.ci = sapply(hajek, function(e) orig.ci[1] < e$estimate && orig.ci[2] > e$estimate),
                             bias = sapply(hajek, function(e) (e$estimate - orig.tx)),
                             # mse = sapply(hajek, function(e) (e$estimate - orig.tx))^2 + sapply(ci.hajek, function(e) e$SD^2),
                            init.wass = init.wass,
                            final.wass = final.wass,
                            n = nsave,
                            estimand = estimand),
                data.frame(estimator = "bp",
                           method = names(weights),
                           estimate = sapply(bp, function(e) e$estimate),
                          sd = sapply(ci.bp, function(e) e$SD),
                          cover.orig = sapply(ci.bp, function(e) e$CI[1] < orig.tx && e$CI[2] > orig.tx),
                          est.in.ci = sapply(bp, function(e) orig.ci[1] < e$estimate && orig.ci[2] > e$estimate),
                          bias = sapply(bp, function(e) (e$estimate - orig.tx)),
                          # mse = sapply(bp, function(e) (e$estimate - orig.tx))^2 + sapply(ci.bp, function(e) e$SD^2),
                          init.wass = init.wass,
                          final.wass = final.wass,
                          n = nsave,
                          estimand = estimand)))
  }
  
  pph_control  <- pph %>% subset(pph[,tx.ind] == 0 & pph[,"sitecode"] != base.site) #control pool
  pph_treated  <- pph %>% subset(pph[,tx.ind] == 1 & pph[,"sitecode"] != base.site) #treated pool
  
  site_treated <- pph %>% subset(pph[,tx.ind] == 1 & pph[,"sitecode"] == base.site)
  site_control <- pph %>% subset(pph[,tx.ind] == 0 & pph[,"sitecode"] == base.site)
  
  original.tx <- mean(site_treated[,outcome] - mean(site_control[,outcome]))
  original.ci <- t.test(x = site_treated[,outcome], y = site_control[,outcome], alternative = "two.sided")$conf.int
  
  treat_data <- rbind(pph_control, site_treated)
  control_data <- rbind(site_control, pph_treated)
  
  weight_control <- weight.by.tx(data = control_data, continuous.covar = continuous.covar, 
                                 binary.covar = binary.covar, outcome = outcome, tx.ind = tx.ind, estimand = "ATC")
  weight_treated <- weight.by.tx(data = treat_data, continuous.covar = continuous.covar, 
                                 binary.covar = binary.covar, outcome = outcome, tx.ind = tx.ind, estimand = "ATT")
  
  
  estimates_control <- estimate.by.tx(data = control_data, weights = weight_control, continuous.covar = continuous.covar, 
                                 binary.covar = binary.covar, outcome = outcome, tx.ind = tx.ind, estimand = "ATC",
                                 orig.tx = original.tx, orig.ci = original.ci)
  estimates_treated <- estimate.by.tx(data = treat_data, weights = weight_treated, continuous.covar = continuous.covar, 
                                 binary.covar = binary.covar, outcome = outcome, tx.ind = tx.ind, estimand = "ATT",
                                 orig.tx = original.tx, orig.ci = original.ci)
  overall.estimates <- rbind(estimates_treated,  estimates_control)
  overall.estimates$sitecode <- base.site
  return(overall.estimates)
  
}

theme_cot <- function(base_size = 11, base_family = "", 
                      base_line_size = base_size/22, 
                      base_rect_size = base_size/22,
                      legend.position = 'bottom',
                      legend.box = "horizontal",
                      legend.justification = "center",
                      legend.margin = ggplot2::margin(0,0,0,0),
                      legend.box.margin = ggplot2::margin(-10,-10,0,-10)) { 
  ggplot2::`%+replace%`(ggplot2::theme_bw(base_size = base_size, base_family = "", 
                                          base_line_size = base_line_size, 
                                          base_rect_size = base_rect_size),
                        ggplot2::theme(
                          plot.title = ggplot2::element_text(hjust = 1),
                          panel.grid.minor = ggplot2::element_blank(),
                          panel.grid.major = ggplot2::element_blank(),
                          strip.background = ggplot2::element_blank(),
                          strip.text.x = ggplot2::element_text(face = "bold"),
                          strip.text.y = ggplot2::element_text(face = "bold"),
                          legend.position = legend.position,
                          legend.box = legend.box,
                          legend.justification = legend.justification,
                          legend.margin = legend.margin,
                          legend.box.margin = legend.box.margin
                          # change stuff here
                        ))
}


#### get estimates for each site ####
set.seed(269570109) #random.org
#site code
# [1] "Cairo, Egypt"   2 "Turkey"        3  "hocmon, Vietnam" 4 "cuchi, Vietnam" 
# [5] "Burkina Faso"
# debugonce(miso_estimate)
egypt <- miso_estimate(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site = 1) #egypt
turkey <- miso_estimate(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site = 2) #turkey
hocmon <- miso_estimate(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site = 3) # hocmon, Vietnam
cuchi <- miso_estimate(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site = 4) # cuchi, Vietnam
burkina <- miso_estimate(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site = 5) # Burkina Faso

#### Get overall estimates ###
overall <- rbind(egypt, turkey, hocmon, cuchi, burkina) %>% 
  mutate(
    # rmse = sqrt(mse),
         method = forcats::fct_relevel(method, "GLM","CBPS","SBW","SCM", "NNM","COT"),
         site = factor(sitecode, labels = c("Egypt","Turkey","Hocmon,\nVietnam","Cuchi,\nVietnam","Burkina Faso"))) 
combined_est <- overall %>% group_by( method, estimator) %>% 
  summarize(
    estimate = weighted.mean(estimate,n),
    bias = weighted.mean(bias,n),
    # rmse = sqrt(weighted.mean(mse, n)),
    cover.orig = mean(cover.orig),
    est.in.ci = mean(est.in.ci),
    n = sum(n)
    )

print(combined_est)

original.ttest <- t.test(x = pph[pph[,tx.ind] == 1,outcome], pph[pph[,tx.ind] == 0,outcome])
original.tau <- original.ttest$estimate[1] - original.ttest$estimate[2]
original.ci <- original.ttest$conf.int
ci.mean.0 <- original.ci - original.tau


overall %>% 
  ggplot(aes(x = method, y = bias, fill = site, color = site, size = n)) +
  geom_point() +
  facet_grid(rows = vars(estimand), cols = vars(estimator)) + theme_cot() +
  theme(panel.grid.major.y = element_line()) + scale_color_jama()

#### plot of estimates ####

pdf(file = "figures/pph.pdf", 
    width = 5.5, height = 2)
print(combined_est %>% 
  filter(estimator == "hajek") %>% 
    mutate(method = forcats::fct_relevel(method,
                                         "COT",
                                         "NNM",
                                         "SCM",
                                         "SBW",
                                         "CBPS",
                                         "GLM"
    )) %>% 
  ggplot(aes(x = estimate, y = method
  )) + 
  geom_vline(xintercept = original.tau, color = "black", alpha = 0.4) +
  geom_vline(xintercept = original.ci[1], color = "black", linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = original.ci[2], color = "black", linetype = 2, alpha = 0.4) +
  geom_point(size = 2)  + 
  scale_x_continuous(minor_breaks = seq(-50,75,25),
                     breaks = seq(-50,75,25)) +
  # scale_color_manual(values = c("#0062ff",rep("black",5))) + 
  theme_bw() +
  ylab("") + theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.major.x = element_blank()
                   ) +
  xlab("difference in blood loss after 20 minutes (mL)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0)
         
  ) 
  )
dev.off()

#### coverage table ####
xtab <- combined_est %>% 
  filter(estimator == "hajek") %>% 
  select(method, cover.orig, est.in.ci) %>% 
  mutate(cover.orig = (100 * cover.orig),
         est.in.ci = (100 * est.in.ci)) %>% 
  rename("\\% C.I. covering original effect" = cover.orig, 
         "\\% of estimates in original C.I." = est.in.ci,
         "Method" = method) %>% 
  xtable(
    caption = "For each method, the table displays the percentage of times that the calculated 95\\% confidence interval (C.I.) covered the true treatment effect and whether the estimated treatment effect was inside the original C.I. from the study. The weighting methods under examination are logistic regression (GLM), Covariate Balancing Propensity Score (CBPS), Stable Balancing Weights (SBW), Synthetic Control Method (SCM), Nearest Neighbor Matching (NNM), and Causal Optimal Transport (COT).",
    label = "tab:pph",
    digits = 0)

align(xtab) <- "ll|p{1.75in}p{1.75in}"
print(xtab, 
      sanitize.colnames.function = function(x){x},
      include.rownames = FALSE,
      file =  "tables/pph.tex")
