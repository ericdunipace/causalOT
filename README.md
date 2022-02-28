# causalOT: Optimal transport methods for causal inference

This R package implements the methods described in [*Optimal transport methods for causal inference*](http://arxiv.org/abs/2109.01991).

## Installation
This package can be installed in a few ways.

### 1. devtools
Using the `remotes` in `R`, one can install the package with 
```
remotes::install_github("ericdunipace/causalOT")
```

### 2. downlaod and install
After downloading the git package using `git clone` or by downloading the .zip file from the button above (Code -> Download Zip), you can install the package with
```
devtools::install("causalOT")
```

## Additional required software
If using the Sinkhorn divergence based optimal transport weights, you will need to intall python 3 and the `geomloss` package, which depens on `numpy`, `scipy`, `torch`, and `pykeops`.

## Useage
The functions in the package are built to construct weights to make distributions more same and estimate causal effects. The primary method we recommend is by using optimal transport weights which balance distributions by design.

### Estimating weights
The weights can be estimated as below with Causal Optimal Transport and using the `calc_weight` function in the package. We select optimal hyperparameters through our bootstrap-based algorithm.
```
library(causalOT)
set.seed(1111)

hainmueller <- Hainmueller$new(n = 128)
hainmueller$gen_data()

weights <- calc_weight(data = hainmueller, method = "Wasserstein",
                       add.divergence = TRUE, grid.search = TRUE,
                       verbose = TRUE)

```
These weights will balance distributions, making estimates of treatment effects the same. We can then estimate effects with 
```
tau_hat <- estimate_effect(data = hainmueller, weights = weights,
                           hajek = TRUE, doubly.robust = FALSE,
                           estimand = "ATE")
```
This creates an object of class `causalEffect` which can be fed into the native `R` function `confint` to calculate asymptotic confidence intervals.
```
ci_tau <- confint(object = tau_hat, level = 0.95, 
                  method = "asymptotic",
                  model = "lm",
                  formula = list(control = "y ~.", treated = "y ~ ."))
```
This then gives the following estimate and C.I.
```
print(tau_hat$estimate)
#  0.2574432
print(ci_tau$CI)
# -0.08814654  0.60303297
```
