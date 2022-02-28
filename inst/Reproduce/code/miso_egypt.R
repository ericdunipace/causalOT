#### Load Packages ####
library(causalOT)
library(forcats)
library(dplyr)
library(ggplot2)
library(xtable)
library(ggsci)

#### Load data ####
data("pph", package = "causalOT")
reticulate::use_python("/usr/local/bin/python3", required = TRUE)
theme_cot <- causalOT:::theme_cot

#### setup treatment indicators ####
bal.cov <- colnames(pph)[-c(1:2)]
outcome <- "cum_blood_20m"
tx.ind <- "tx"

#### comparisons ####
orig.ci <- c(-27.1 ,   27.8 )
orig.tx <- 0.369

#### cost matrices ####
power <- 2
z <- pph[,tx.ind]
n0 <- sum(1 - z)
n1 <- sum(z)

continuous.covar <- c("age","gest_age","num_livebirth","hb_test","bloodlossattx")
bin.covar <- colnames(pph)[sapply(colnames(pph), function(x) 2 == length(unique(pph[,x])))][-1]

cost <- cost_fun( x = cbind(scale(pph[,continuous.covar], scale = FALSE)
                            %*% causalOT:::inv_sqrt_mat(cov(pph[,continuous.covar]))
                            ,  pph[,bin.covar]), z = z, metric = "Lp", power = power,
                  estimand = "ATT")
scaled_pph <- cbind(scale(pph[,continuous.covar], scale = FALSE)
                   %*% causalOT:::inv_sqrt_mat(cov(pph[,continuous.covar]))
                   ,  pph[,bin.covar], tx = z)
colnames(scaled_pph) <- c(continuous.covar, bin.covar,"tx")

#### a priori balance ####
std.mean.diff <- mean_bal(data = pph, treatment.indicator = tx.ind, 
                          balance.covariates = bal.cov,
                          outcome = outcome)

sw <- causalOT:::get_sample_weight(sw = NULL, z = z)
wass.diff <- sinkhorn_geom(a = sw$a, b = sw$b, blur = 1e3,
                           x = scaled_pph[z == 0,c(continuous.covar, bin.covar)],
                           y = scaled_pph[z == 1,c(continuous.covar, bin.covar)],
                           p = power)$loss


#### Calculate weights ####
#tune mean balance with sbw
set.seed(652180044)
SBW_tune  = calc_weight(data = pph, 
                           estimand = "ATT", 
                           method = "SBW",
                           grid.search = TRUE,
                           solver = "mosek",
                           grid.length = 20,
                           grid = seq(0, 2, length.out = 20),
                           formula = "~. + .*. + I(.^2) + 0",
                           treatment.indicator = tx.ind, 
                           balance.covariates = bal.cov,
                           outcome = outcome)

weights <- list()
set.seed(604274520) #random.org
weights$COT =  calc_weight(data = scaled_pph, 
                           estimand = "ATT", 
                           method = "Wasserstein",
                           solver = "mosek",
                           grid.search = TRUE,
                           p = power,
                           metric = "Lp",
                           add.divergence = TRUE,
                           penalty = "entropy",
                           treatment.indicator = tx.ind, 
                           balance.covariates = bal.cov,
                           wass.method = "sinkhorn",
                           n.boot = 1e2, 
                           tol = 1e-10, stepsize = 1e-2,
                           verbose = TRUE)

set.seed(-887238551) #random.org
weights$COT_m =  calc_weight(data = scaled_pph, 
                          estimand = "ATT", 
                          method = "Wasserstein",
                          solver = "mosek",
                          grid.search = TRUE,
                          p = power,
                          metric = "Lp",
                          cost = cost,
                          add.divergence = TRUE,
                          penalty = "entropy",
                          treatment.indicator = tx.ind, 
                          balance.covariates = bal.cov,
                          wass.method = "sinkhorn",
                          n.boot = 1e2, 
                          formula = "~. + 0",
                          balance.constraints = SBW_tune$args$constraint,
                          tol = 1e-10, stepsize = 1e-2,
                          verbose = TRUE)


weights$COT$gamma <- NULL
weights$COT_m$gamma <- NULL
for (i in 1:length(weights)) {
  weights[[i]]$gamma <- causalOT:::calc_gamma(weights[[i]], cost = cost, p = 1, niter = 0)
}


#### a posteriori balance ####
std.mean.diff.post <- sapply(weights, mean_bal,
                             data = pph, treatment.indicator = tx.ind, 
                             balance.covariates = bal.cov,
                             outcome = outcome)
wass.diff.post <- sapply(weights, function(w) sinkhorn_geom(a = w$w0, b = w$w1, blur = 1e3, metric = "Lp",
                                                            x = scaled_pph[z == 0,c(continuous.covar, bin.covar)],
                                                            y = scaled_pph[z == 1,c(continuous.covar, bin.covar)],
                                                            p = power)$loss)


#### Plot change in Balance ####
pretty.names <- forcats::fct_recode(factor(names(std.mean.diff)),
                                    "Age" = "age",
                                    "No education" = "no_educ",
                                    "Num. live births" = "num_livebirth",
                                    "Currently married" = "cur_married",
                                    "Gestational age" = "gest_age",
                                    "Previous PPH" = "prev_pphyes",
                                    "Hemoglobin" = "hb_test",
                                    "Labor induced" = "induced_laboryes",
                                    "Labor augmented" = "augmented_laboryes",
                                    "Early cord clamp" = "early_cordclampyes",
                                    "Control cord traction" = "control_cordtractionyes",
                                    "Uterine massage" = "uterine_massageyes",
                                    "Placenta delivered" = "placenta",
                                    "Blood loss at treatment" = "bloodlossattx")
pretty.names <- factor(pretty.names, levels = sort(as.character(pretty.names)))

mean.change <- data.frame(pre = std.mean.diff,
                          post = c(std.mean.diff.post[,c("COT","COT_m")]),
                          method = rep(c("COT","COT_m"), each = nrow(std.mean.diff.post)),
                          variable = rep(pretty.names, 2),
                          Balance = factor(rep("Before", nrow(std.mean.diff.post) * 2), levels = c("Before", "After"))
                          )
mc2 <- mean.change
mc2$Balance <- factor(rep("After", nrow(std.mean.diff.post) * 2), levels = c("Before", "After"))

wass.change <- data.frame(pre = sqrt(wass.diff),
                          post = sqrt(c(wass.diff.post[c("COT","COT_m")])),
                          method = c("COT","COT_m"),
                          Balance =  factor(rep("Before", 2), levels = c("Before", "After"))
                          
)
wchange2 <- wass.change
wchange2$Balance <- factor(rep("After",  2), levels = c("Before", "After"))

wass.plot <-
  wass.change %>% 
  mutate(method = forcats::fct_recode(method, "COT, means" = "COT_m")) %>% 
  ggplot(aes( xmin = pre, xmax = post, 
             x = pre, y = method,
             # yend = variable,
             shape = Balance,
             # color = constraints
             )) +
  geom_linerange(position = position_dodge2(width = 0,
                                                          reverse = TRUE),
                 show.legend = FALSE
                 # ,color = rev(pal_jama()(4))
                 ) +
  geom_point(size = 2, 
             # color = rev(pal_jama()(2)),
             position = position_dodge2(width = 0, reverse = TRUE)) +
  geom_point(data = wchange2 %>% 
               mutate(method = forcats::fct_recode(method, "COT, means" = "COT_m")), 
             aes(
                 x = post,
                 y = method,
                 # yend = variable,
                 
                 shape = Balance), 
             # shape = 2,
             size = 2, 
             # color = rev(pal_jama()(2)),
             position = position_dodge2(width = 0, reverse = TRUE),
             show.legend = FALSE) +
    scale_x_continuous(limits = c(0,2),
                       expand = c(0,0),
                       breaks = seq(0,1.5,0.5),
                       name = "2-Sinkhorn Divergence") +
    scale_y_discrete(name = "Constraints") + 
    scale_shape_manual(
      values = c(16,2),
                name = "Balance",
                         labels = c("After","Before" )) +
    theme_bw() +
    guides(
           shape = guide_legend(order = 1, nrow = 1)) +
    theme(legend.box = "horizontal",
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,0,-10),
          legend.position = "bottom",
          legend.title.align = 0.5,
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())

pph.plot <-
  mean.change %>% 
  mutate(method = forcats::fct_recode(method, "COT, means" = "COT_m")) %>% 
  ggplot(aes(y = variable, xmin = pre, xmax = post, 
             x = pre,
             # yend = variable,
             color = method,
             shape = Balance)) +
  geom_linerange(position = position_dodge2(width = .6,
                                            reverse = TRUE),
                 show.legend = FALSE) +
  geom_point(size = 2, position = position_dodge2(width = .6, reverse = TRUE)) +
  scale_color_jama(name = "Constraints") +
  geom_point(data = mc2 %>% 
               mutate(method = forcats::fct_recode(method, "COT, means" = "COT_m")), aes(y = variable, 
                             x = post,
                             color = method,
                             shape = Balance), 
             # shape = 2,
            
             size = 2, 
             position = position_dodge2(width = .6, reverse = TRUE),
             show.legend = FALSE) +
  scale_x_continuous(limits = c(0,2.2),
                     expand = c(0.01,0),
                     breaks = seq(0,2.2,0.2),
                     name = "Standardized difference in means") +
  ylab("Variable") + 
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(
    values = c(16,2),
    name = "Balance",
    labels = c("After","Before" )) +
  theme_cot() +
  guides(colour = guide_legend(order = 2, nrow = 2), 
         shape = guide_legend(order = 1, nrow = 2)) +
  theme(legend.box = "horizontal",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.position = "bottom",
        legend.title.align = 0.5)


# save figures
pdf(file = "figures/miso_mb.pdf", 
    width = 5.5, height = 4)
print(pph.plot)
dev.off()

pdf(file = "figures/miso_wp.pdf", 
    width = 5.5, height = 1.85)
print(wass.plot)
dev.off()

#### calculate estimates ####

fill <- rep(NA_real_, length(weights))
estimates <- data.frame(`Linear` = fill,
                        `Linear, mapped` = fill,
                        `Augmented, gp` = fill,
                        # `Augmented, mapped` = fill,
                        `Weighted OLS` = fill,
                        check.names = FALSE)
rownames(estimates) <- c(
                         "COT",
                         "COT, means")

lin.est <- lapply(weights,
                 function(w) estimate_effect(data = pph, weights = w, 
                                             matched = FALSE, 
                                             split.model = TRUE,
                                             doubly.robust = FALSE,
                                             estimand = "ATT",
                                             treatment.indicator = tx.ind, 
                                             balance.covariates = bal.cov,
                                             outcome = outcome))

lin.map.est <- lapply(weights,
                      function(w) estimate_effect(data = pph, weights = w, 
                                                  matched = TRUE, 
                                                  split.model = TRUE,
                                                  doubly.robust = FALSE,
                                                  estimand = "ATT",
                                                  p = 1,
                                                  cost = cost,
                                                  treatment.indicator = tx.ind, 
                                                  balance.covariates = bal.cov,
                                                  outcome = outcome))

wols.est <- lapply(weights,
                                   function(w) estimate_effect(data = pph, weights = w, 
                                                               matched = FALSE, 
                                                               split.model = FALSE,
                                                               estimand = "ATT",
                                                               treatment.indicator = tx.ind, 
                                                               balance.covariates = bal.cov,
                                                               outcome = outcome))

dr.gp.est <- lapply(weights, #note this requires rstan
                 function(w) estimate_effect(data = cbind(scaled_pph, cum_blood_20m = pph[,"cum_blood_20m"]), 
                                             weights = w, 
                                             matched = FALSE, 
                                             estimand = "ATT",
                                             model = "gp",
                                             treatment.indicator = tx.ind, 
                                             balance.covariates = bal.cov,
                                             outcome = outcome))


estimates$Linear <- sapply(lin.est,
                  function(w) w$estimate)

estimates$`Linear, mapped` <- sapply(lin.map.est,
                    function(w) w$estimate)


estimates$`Augmented, gp` <- sapply(dr.gp.est,
                              function(w) w$estimate)


estimates$`Weighted OLS` <- sapply(wols.est,
                                   function(w) w$estimate)


print(estimates)


#### calculate CI ####

  set.seed(137567173)
  ci.linear <- lapply(lin.est, confint, verbose = TRUE, method = "asymptotic",
                      model = "gp", formula = list(treated = "y ~ .",
                                                   control = "y ~ ."))
  
  set.seed(862020444)
  ci.linear.map <- lapply(lin.map.est, confint, verbose = TRUE, method = "asymptotic",
                          model = "gp", formula = list(treated = "y ~ .",
                                                       control = "y ~ ."))
  
  set.seed(926512351)
  ci.wols <- lapply(wols.est, confint, method = "asymptotic")
  
  
  set.seed(851708157)
  ci.dr.gp <- lapply(dr.gp.est, confint, verbose = TRUE, method = "asymptotic",
                  model = "gp")
  
  
  ci.obj <- list(linear = ci.linear,
                 linear.map = ci.linear.map,
                 dr.gp = ci.dr.gp,
                 wols = ci.wols)
  
  

#### Combine estimates and CI ####
est.mat <- as.data.frame(matrix(NA_character_, nrow = length(lin.est),
                  ncol = 2*length(ci.obj)))
rownames(est.mat) <- rownames(estimates)
colnames(est.mat) <- NULL
colnames(est.mat)[seq(1,2*length(ci.obj),2)] <- colnames(estimates)

est.mat[,1] <- estimates$Linear
est.mat[,2] <- sapply(ci.obj$linear, function(ci) paste0("(", paste(format(round(ci$CI, digits = 1), digits = 4,nsmall = 1), collapse = ", "), ")"))

est.mat[,3] <- estimates$`Linear, mapped`
est.mat[,4] <- sapply(ci.obj$linear.map, function(ci) paste0("(", paste(format(round(ci$CI, digits = 1), digits = 4,nsmall = 1), collapse = ", "), ")"))


est.mat[,5] <- estimates$`Weighted OLS`
est.mat[,6] <- sapply(ci.obj$wols, function(ci) paste0("(", paste(format(round(ci$CI, digits = 1), digits = 4,nsmall = 1), collapse = ", "), ")"))

est.mat[,7] <- estimates$`Augmented, gp`
est.mat[,8] <- sapply(ci.obj$dr.gp, function(ci) paste0("(", paste(format(round(ci$CI, digits = 1), digits = 4,nsmall = 1), collapse = ", "), ")"))


for (i in seq(1,8,2)) {
  est.mat[[i]] <- round(as.numeric(est.mat[[i]]), digits = 1)
}
colnames(est.mat) <- 1:ncol(est.mat)
est.tab <- est.mat %>% mutate(method = c(
                                         "COT",
                                         "COT, means")) %>% 
  relocate(c("method"))

colnames(est.tab) <- c("method","Hajek","",
                       "Barycentric Projection", "", 
                       "Weighted OLS", "", 
                       "Augmented, GP", "")

x.est <- xtable(est.tab,
                digits = rep(1, ncol(est.tab) + 1),
                align = "ll|lrlrlrlr|")

cn.est <- colnames(est.tab)
addtorow.est <- list()
addtorow.est$pos <- list(0)
addtorow.est$command <- paste0(paste0('& & & \\multicolumn{2}{c}{', cn.est[!is.na(cn.est)], '}', collapse=''), '\\\\')

#add to figures folder
pap.tab <- est.tab[
                   ,c(1,2:3, 8:9, 6:7, 4:5)]
colnames(pap.tab)[seq(3,9,2)] <- ""
colnames(pap.tab)[c(1, 4)] <- c("", "Augmented")
addtorow.pap <- list(pos = list(-1), 
                     command = paste0("\\hline", paste0('& \\multicolumn{2}{c}{', c("Hajek", "Augmented", "Weighted OLS", "Barycentric projection"), '}', collapse=''), '\\\\'))
colnames(pap.tab) <- c("Method", rep(c("Est.", "C.I."), 4))
  print(xtable(pap.tab,
               digits = rep(1, ncol(pap.tab) + 1),
               align = "ll|rlrlrlrl",
               caption = "Estimates and confidence intervals for optimal transport weighting methods applied to a modification of the data in \\cite{Blum2010}. We have constructed a control group for the misoprostol receiving patients at Egyptian site using the controls from the four other sites. The original treatment effect at the Egypt site was $0.369$ mL with a 95\\% C.I. of $(-27.1,  27.8)$. The augmented method uses a gaussian process as the outcome model.",
               label = "tab:miso"), add.to.row=addtorow.pap, include.colnames=TRUE, include.rownames = FALSE,
        file = "tables/miso.tex",
        hline.after = c(0, 2),
        table.placement = "tbh",
        size="\\fontsize{9pt}{10pt}\\selectfont")

#### est and ci plot ####
est.df <- data.frame(estimate = c(estimates$Linear,
                                  estimates$`Augmented, gp`,
                                  estimates$`Weighted OLS`,
                                  estimates$`Linear, mapped`),
                     ci.lwr = c(sapply(ci.obj$linear, `[[`,"CI")[1,],
                                sapply(ci.obj$dr.gp, `[[`,"CI")[1,],
                                sapply(ci.obj$wols, `[[`,"CI")[1,],
                                sapply(ci.obj$linear.map, `[[`,"CI")[1,]
                     ),
                     ci.upr = c(sapply(ci.obj$linear, `[[`,"CI")[2,],
                                sapply(ci.obj$dr.gp, `[[`,"CI")[2,],
                                sapply(ci.obj$wols, `[[`,"CI")[2,],
                                sapply(ci.obj$linear.map, `[[`,"CI")[2,]),
                     "weight method" = c("COT", "COT, means"),
                     "estimator" = rep(c("Hajek","Augmented"
                                         , "WOLS"
                                         ,"Barycentric Projection"), each = 2),
                     check.names = FALSE
)



# put in paper file
pdf(file = "figures/miso_att.pdf", 
    width = 5.5, height = 2)
print(est.df %>% 
        mutate(estimator = forcats::fct_relevel(factor(estimator), "Hajek", "Augmented",
                                                "WOLS",
                                                "Barycentric Projection")) %>% 
        mutate(estimator = forcats::fct_recode(estimator, "Barycentric\nProjection" = "Barycentric Projection")) %>% 
        ggplot(aes(x = estimate, y = `estimator`,
                   color = `weight method`
        )) + 
        geom_vline(xintercept = orig.tx, color = "black", alpha = 0.4) +
        geom_vline(xintercept = orig.ci[1], color = "black", linetype = 2, alpha = 0.4) +
        geom_vline(xintercept = orig.ci[2], color = "black", linetype = 2, alpha = 0.4) +
        geom_point(position = position_dodge2(width = 0.6, reverse = TRUE),
                   size = 2)  + 
        geom_linerange(aes(xmin = ci.lwr,
                           xmax = ci.upr),
                       position = position_dodge2(width = 0.6,
                                                  reverse = TRUE),
                       size = 0.75, show.legend = FALSE) +
        scale_color_jama() + 
        theme_bw() +
        ylab("") + theme(panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank()) +
        xlab("difference in blood loss after 20 minutes (mL)") +
        theme(plot.title = element_text(hjust = 0.5),
                            legend.position = "right",
                            legend.box = "vertical",
                            legend.box.margin = margin(0,0,0,0)
        ) )
dev.off()
