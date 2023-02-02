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
# reticulate::use_python("/usr/bin/python3", required = TRUE)

#### setup treatment indicators ####
bal.cov <- colnames(pph)[-c(1:2)]
outcome <- "cum_blood_20m"
tx.ind <- "tx"
continuous.covar <- c("age","gest_age","num_livebirth","hb_test","bloodlossattx")
binary.covar <- colnames(pph)[sapply(colnames(pph), function(x) 2 == length(unique(pph[,x])))][-1]

#### Estimation function ####
miso_estimate <- function(pph, continuous.covar, binary.covar, outcome, tx.ind, base.site) {
  
  weight.by.tx <- function(data, continuous.covar, binary.covar, outcome, tx.ind, estimand = "ATT") {
    # browser()
    inv_sqrt_mat <- function(X, symmetric = FALSE) {
      p <- ncol(X)
      decomp <- eigen(as.matrix(X), symmetric = symmetric)
      return(tcrossprod(decomp$vectors %*% diag(1/sqrt(abs(decomp$values)), p, p), decomp$vectors))
    }
    
    bal.cov <- c(continuous.covar, binary.covar)
    power <- 2
    z <- data[,tx.ind]
    n0 <- sum(1 - z)
    n1 <- sum(z)
    
    center    <- switch(estimand,
                        "ATC" = colMeans(data[z==0,continuous.covar]),
                        "ATT" = colMeans(data[z==1,continuous.covar]))
    
    covar     <- switch(estimand,
                        "ATC" = cov(data[z==0,continuous.covar]),
                        "ATT" = cov(data[z==1,continuous.covar]))
    
    scaled_x  <- cbind(scale(data[,continuous.covar], center = center, scale = FALSE)
                       %*% inv_sqrt_mat(covar)
                       ,  data[,binary.covar])
    x <- data[,bal.cov]
    
    newdat <- data.frame(cbind(z = z, x))
    
    sq_form <- paste0("z ~. + .*. + ", paste0("I(",colnames(x),"^2)", collapse = " + "))
    form <- "z ~."
    
    weights <- list()
    weights$GLM  = calc_weight(x = df2dataHolder(treatment.formula = sq_form,
                                                 data = newdat
                                                 ), 
                               estimand = estimand, method = "Logistic")
    
    weights$CBPS  = calc_weight(x = df2dataHolder(treatment.formula = sq_form,
                                                  data = newdat
                                                  ), 
                                estimand = estimand, method = "CBPS")
    
    
    weights$SBW  = calc_weight(x = df2dataHolder(treatment.formula = sq_form,
                                                 data = newdat
                                                ), 
                              estimand = estimand, method = "SBW", 
                              options = list(verbose = FALSE,
                                             delta = seq(0, 2, length.out = 20)))
    
    
    weights$SCM  = calc_weight(x = scaled_x, z = z,
                               estimand = estimand, method = "SCM",
                               options = list(verbose = FALSE))
    
    weights$NNM  = calc_weight(x = scaled_x, z = z,
                               estimand = estimand, method = "NNM")
    # browser()
    # debugonce(causalOT:::grid_select, signature = signature(object = "ANY", w = "list"))
    # debugonce(causalOT:::cot_solve, signature = signature(object = "gridSearch"))
    weights$COT =  calc_weight(x = scaled_x, z = z,
                               estimand = estimand, method = "COT",
                               options = list(line_search_fn = "strong_wolfe"))
    
    # m1 <- Measure(x = scaled_x[z==1,], adapt = "weights")
    # m0 <- Measure(x = scaled_x[z==0,])
    # m  <- Measure(x = scaled_x)
    # OT1 <- OTProblem(m1, m0)
    # OT1$set_ot_arguments(debias = TRUE)
    # OT1$solve(torch_args = list(line_search_fn = "strong_wolfe"))
    # OT1$choose_hyperparameters()
    # OT1$info()$hyperparam.metrics
    # OT1$selected_lambda
    # OT1$.__enclos_env__$private$final_loss
    
    return(weights)
    
  }
  estimate.by.tx <- function(data, weights, continuous.covar, binary.covar, outcome, tx.ind, estimand = "ATT", orig.tx, orig.ci) {
    # browser()
    inv_sqrt_mat <- function(X, symmetric = FALSE) {
      p <- ncol(X)
      decomp <- eigen(as.matrix(X), symmetric = symmetric)
      return(tcrossprod(decomp$vectors %*% diag(1/sqrt(abs(decomp$values)), p, p), decomp$vectors))
    }
    
    bal.cov <- c(continuous.covar, binary.covar)
    power <- 2
    z <- data[,tx.ind]
    n0 <- sum(1 - z)
    n1 <- sum(z)
    center    <- switch(estimand,
                        "ATC" = colMeans(data[z==0,continuous.covar]),
                        "ATT" = colMeans(data[z==1,continuous.covar]))
    
    covar     <- switch(estimand,
                        "ATC" = cov(data[z==0,continuous.covar]),
                        "ATT" = cov(data[z==1,continuous.covar]))
    
    scaled_x  <- cbind(scale(data[,continuous.covar], center = center, scale = FALSE)
                       %*% inv_sqrt_mat(covar)
                       ,  data[,binary.covar])
    x <- data[,bal.cov]
    y <- data[,outcome]
    nsave <- nrow(x)
    
    max_cost <- sum( (apply(scaled_x,2,min) - 
                        apply(scaled_x,2,max))^2 )
    
    newdat <- data.frame(cbind(z = z, x))
    
    hajek <- lapply(weights, estimate_effect,
                    x = x,
                    y = y,
                    model.function = NULL,
                    estimate.separately = TRUE,
                    aurgment.estimate = FALSE,
                    normalize.weights = TRUE)
    bp <- lapply(weights, estimate_effect,
                 x = scaled_x,
                 y = y,
                 p = 3,
                 model.function = barycentric_projection,
                 estimate.separately = TRUE,
                 aurgment.estimate = FALSE,
                 normalize.weights = TRUE,
                 line_search_fn = "strong_wolfe")
    
    ci.hajek <-  lapply(hajek, confint)
    
    ci.bp <- lapply(bp, confint)
    
    sd.hajek <-  sqrt(sapply(hajek, vcov))
    
    sd.bp <- sqrt(sapply(bp, vcov))
    
    
    init.wass <- ot_distance(x1 = scaled_x[z==1,], x2 = scaled_x[z==0,],
               a = rep(1/n1,n1), b = rep(1/n0,n0),
               p = 2, penalty = max_cost, diameter = max_cost)
  
    final.wass <- sapply(weights, function(w) {
      ot_distance(x1 = scaled_x[z==1,], x2 = scaled_x[z==0,],
                  a = w@w1, b = w@w0,
                  p = 2, penalty = max_cost, diameter = max_cost)
    })
    
    return(rbind(data.frame(estimator = "hajek",
                            method = names(weights),
                            estimate = sapply(hajek, coef),
                             sd = c(sd.hajek),
                             cover.orig = sapply(ci.hajek, function(e) e[1] < orig.tx && e[2] > orig.tx),
                             est.in.ci = sapply(hajek, function(e) orig.ci[1] < coef(e) && orig.ci[2] > coef(e)),
                             bias = sapply(hajek, function(e) (coef(e) - orig.tx)),
                             # mse = sapply(hajek, function(e) (e$estimate - orig.tx))^2 + sapply(ci.hajek, function(e) e$SD^2),
                            init.wass = init.wass,
                            final.wass = final.wass,
                            n = nsave,
                            estimand = estimand),
                data.frame(estimator = "bp",
                           method = names(weights),
                           estimate = sapply(bp, coef),
                          sd = c(sd.bp),
                          cover.orig = sapply(ci.bp, function(e) e[1] < orig.tx && e[2] > orig.tx),
                          est.in.ci = sapply(bp, function(e) orig.ci[1] < coef(e) && orig.ci[2] > coef(e)),
                          bias = sapply(bp, function(e) (coef(e) - orig.tx)),
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

#### distributional balance ####