
supported_methods <- function() {
  c("Wasserstein", likelihood_methods(), grid_search_methods())
}

likelihood_methods <- function() {
  c("Logistic","Probit","CBPS")
}

grid_search_methods <- function() {
  c(balancing_function_methods(), balancing_distributions_methods())
}

balancing_function_methods <- function() {
  c("SBW","EntropyBW")
}

balancing_distributions_methods <- function() {
  c("SCM", cot_methods())
}

cot_methods <- function() {
  c("COT","EnergyBW", "NNM")
}