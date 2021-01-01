#### Load package ####
library(causalOT)

#### Set Seed ####
set.seed(245153142) #from random.org

#### Set Up Data ####
#Load LaLonde R6 generator
lalonde <- LaLonde$new(design = "Full")

#Just gets the fixed data
lalonde$gen_data()


#### Calculate Weights ####
weights <- list(naive = list(w0 = rep(1, lalonde$get_n()[1]),
                             w1 = rep(1, lalonde$get_n()[2])),
                ipw = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "Logistic", formula = "z~."),
                sbw = calc_weight(lalonde, constraint = 0,
                                  estimand = "ATT", formula = "~.+0",
                                  method = "SBW", solver = "mosek"),
                # in principle can use grid search but it just selects
                # constraint = 0, so we can save time by just using
                # that constraint. To check, uncomment the lines below
                # and comment the sbw lines above
                # sbw = calc_weight(lalonde, constraint = NULL, 
                #                   estimand = "ATT", formula = "~.+0",
                #                   method = "SBW", grid.search = TRUE,
                #                   grid.length = 20, solver = "mosek",
                #                   n.boot = 1000),
                nnm  = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "NNM", p = 1, metric = "mahalanobis"),
                ot  = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "Constrained Wasserstein", grid.search = TRUE,
                                  p = 1, metric = "mahalanobis", solver = "mosek",
                                  wass.method = "greenkhorn", wass.iter = 1e4,
                                  n.boot = 1000),
                mot = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "Wasserstein", 
                                  grid.search = TRUE, add.joint = TRUE,
                                  p = 1, metric = "mahalanobis", 
                                  solver = "mosek",
                                  wass.method = "greenkhorn", wass.iter = 1e4,
                                  n.boot = 1000),
                mot.means = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "Wasserstein", grid.search = TRUE, add.joint = TRUE,
                                  formula = "~.+0",
                                  balance.constraints = 0.1, 
                                  p = 1, metric = "mahalanobis", 
                                  solver = "mosek",
                                  wass.method = "greenkhorn", wass.iter = 1e4,
                                  n.boot = 1000
                                  )
)

#### Calculate ATT ####
# raw weight estimators
effects <- lapply(weights, estimate_effect, data = lalonde,
                  hajek = TRUE, doubly.robust = FALSE,
                     target = "ATT", split.model = TRUE,
                     formula = NULL)

# if not using mahlanobis with power = 2, can run this and expect
# different results
match.effects <- lapply(weights, estimate_effect, data = lalonde,
                        hajek = TRUE, doubly.robust = FALSE,
                        matched = TRUE,
                        target = "ATT", split.model = TRUE,
                        formula = NULL)
names(match.effects) <- paste("matched",names(match.effects),sep = ".")
estimates <- cbind(sapply(effects, function(e) e$estimate),
                   sapply(match.effects, function(e) e$estimate))
colnames(estimates) <- c("Linear","Matched/Mapped")

# doubly robust lm/weight
effects.dr <- lapply(weights, estimate_effect, data = lalonde,
                     hajek = TRUE, doubly.robust = TRUE,
                     target = "ATT", model = "lm", split.model = TRUE,
                     formula = NULL)

match.effects.dr <- lapply(weights, estimate_effect, data = lalonde,
                        hajek = TRUE, doubly.robust = TRUE,
                        matched = TRUE,
                        target = "ATT", split.model = TRUE,
                        formula = NULL)
names(match.effects.dr) <- paste("matched",names(match.effects.dr),sep = ".")
estimates.dr <- cbind(sapply(effects.dr, function(e) e$estimate),
                   sapply(match.effects.dr, function(e) e$estimate))
colnames(estimates.dr) <- c("Linear","Matched/Mapped")


# weighted lm
effects.lm <- lapply(weights, estimate_effect, data = lalonde,
                  hajek = TRUE, doubly.robust = FALSE,
                  target = "ATT", model = "lm", split.model = FALSE,
                  formula = NULL)

estimates.lm <- sapply(effects.lm, function(e) e$estimate)
names(estimates.lm) <- names(effects.lm)
estimates.lm <- as.matrix(estimates.lm)
colnames(estimates.lm) <- "Weighted OLS"

#print all estimates
print(estimates)
print(estimates.dr)
print(as.matrix(estimates.lm))


#### Double check implementation ####
# library(sbw)
# dat <- lalonde_full %>% select(!data_id)
# 
# # debugonce(sbw)
# # sbw_compare <- sbw(dat = dat, ind = "treat", out = "re78",
# #                    bal = list(bal_cov = c('age','education','black','hispanic','married','nodegree','re74','re75'),
# #                               bal_tol = 0, bal_alg = FALSE),
# #                    sol = list(sol_nam = "mosek"))
# sbw_wt_comparerds <- readRDS("ll_weights.rds")
# sbw_wt_compare <- list(w0 = sbw_wt_comparerds[dat$treat==0],
#                        w1 = sbw_wt_comparerds[dat$treat==1])
# sum(dat$re78[dat$treat == 0] * sbw_wt_compare$w0)
# mean(dat$re78[dat$treat == 1])
# mean(dat$re78[dat$treat == 1]) - sum(dat$re78[dat$treat == 0] * sbw_wt_compare$w0)
# 
# # debugonce(sbw)
# # sbw_compare <- sbw(dat = dat, ind = "treat", out = "re78",
# #                    bal = list(bal_cov = c('age','education','black','hispanic','married','nodegree','re74','re75'),
# #                               bal_alg = TRUE),
# #                    sol = list(sol_nam = "mosek"))
# sbw_wt_comparerds2 <- readRDS("ll_weights_tune.rds")
# sbw_wt_compare2 <- list(w0 = sbw_wt_comparerds2[dat$treat==0],
#                        w1 = sbw_wt_comparerds2[dat$treat==1])
# sum(dat$re78[dat$treat == 0] * sbw_wt_compare2$w0)
# mean(dat$re78[dat$treat == 1])
# mean(dat$re78[dat$treat == 1]) - sum(dat$re78[dat$treat == 0] * sbw_wt_compare2$w0)
# 
