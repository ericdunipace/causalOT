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
weights <- list(ipw = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                   method = "Logistic"),
                sbw = calc_weight(lalonde, constraint = NULL, 
                                  estimand = "ATT",
                                   method = "SBW", grid.search = TRUE,
                                   grid.length = 100, solver = "mosek"),
                nnm  = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "NNM", p = 2, metric = "mahalanobis"),
                ot  = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                   method = "Constrained Wasserstein", grid.search = TRUE,
                   p = 2, metric = "mahalanobis", solver = "mosek",
                   wass.method = "greenkhorn", wass.iter = 1000),
                mot = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                   method = "Wasserstein", 
                                  grid.search = TRUE, add.joint = TRUE,
                                   solver = "mosek",
                                   wass.method = "greenkhorn", wass.iter = 1000),
                mot.means = calc_weight(lalonde, constraint = NULL, estimand = "ATT",
                                  method = "Wasserstein", grid.search = TRUE, add.joint = TRUE,
                                  formula = "~.+0",
                                  balance.constraints = 0.1, 
                                  solver = "mosek",
                                  wass.method = "greenkhorn", wass.iter = 1000
                                  )
)

#### Calculate ATT ####

