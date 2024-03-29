% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treatment_effect.R
\name{coef.causalEffect}
\alias{coef.causalEffect}
\title{Extract treatment effect estimate}
\usage{
\method{coef}{causalEffect}(object, ...)
}
\arguments{
\item{object}{An object of class \linkS4class{causalEffect}}

\item{...}{Not used}
}
\value{
A number corresponding to the estimated treatment effect
}
\description{
Extract treatment effect estimate
}
\examples{
# set-up data
set.seed(1234)
data <- Hainmueller$new()
data$gen_data()

# calculate quantities
weight <- calc_weight(data, method = "Logistic", estimand = "ATE")
tx_eff <- estimate_effect(causalWeights = weight)

all.equal(coef(tx_eff), c(estimate = tx_eff@estimate))
}
