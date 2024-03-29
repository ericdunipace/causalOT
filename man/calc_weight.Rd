% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_weight.R
\name{calc_weight}
\alias{calc_weight}
\title{Estimate causal weights}
\usage{
calc_weight(
  x,
  z,
  estimand = c("ATC", "ATT", "ATE"),
  method = supported_methods(),
  options = NULL,
  weights = NULL,
  ...
)
}
\arguments{
\item{x}{A numeric matrix of covariates. You can also pass an object of class \link{dataHolder} or \link{DataSim}, which will make argument \code{z} not necessary,}

\item{z}{A binary treatment indicator.}

\item{estimand}{The estimand of interest. One of "ATT","ATC", or "ATE".}

\item{method}{The method to estimate the causal weights. Must be one of the methods returned by \code{\link[=supported_methods]{supported_methods()}}.}

\item{options}{The options for the solver. Specific options depend on the solver you will be using and you can use the solver specific options functions as detailed below..}

\item{weights}{The sample weights. Should be \code{NULL} or have a weight for each observation in the data. Normalized to sum to one.}

\item{...}{Not used at this time.}
}
\value{
An object of class \link[=causalWeights-class]{causalWeights}
}
\description{
Estimate causal weights
}
\details{
We detail some of the particulars of the function arguments below.
\subsection{Causal Optimal Transport (COT)}{

This is the.main method of the package. This method relies on various solvers depending on the particular options chosen. Please see \link[=cotOptions]{cotOptions()} for more details.
}

\subsection{Energy Balancing Weights (EnergyBW)}{

This is equivalent to COT with an infinite penalty parameter, \code{options(lambda = Inf)}. Uses the same solver and options as COT, \link[=cotOptions]{cotOptions()}.
}

\subsection{Nearest Neighbor Matching with replacement (NNM)}{

This is equivalent to COT with a penalty parameter = 0, \code{options(lambda = 0)}. Uses the same solver and options as COT, \link[=cotOptions]{cotOptions()}.
}

\subsection{Synthetic Control Method (SCM)}{

The SCM method is equivalent to an OT problem from a different angle. See \code{\link[=scmOptions]{scmOptions()}}.
}

\subsection{Entropy Balancing Weights (EntropyBW)}{

This method balances chosen functions of the covariates specified in the data argument, \code{x}. See \code{\link[=entBWOptions]{entBWOptions()}} for more details. Hainmueller (2012).
}

\subsection{Stable Balancing Weights (SBW)}{

Entropy Balancing Weights with a different penalty parameter, proposed by Zuizarreta (2012). See \code{\link[=sbwOptions]{sbwOptions()}} for more details
}

\subsection{Covariate Balancing Propensity Score (CBPS)}{

The CBPS method of Imai and Ratkovic. Options argument is passed to the function \link[CBPS:CBPS]{CBPS()}.
}

\subsection{Logistic Regression or Probit Regression}{

The main methods historically for implementing inverse probability weights. Options are passed directly to the \code{glm} function from \code{R}.
}
}
\examples{
set.seed(23483)
n <- 2^5
p <- 6
#### get data ####
data <- Hainmueller$new(n = n, p = p)
data$gen_data()
x <- data$get_x()
z <- data$get_z()

if (torch::torch_is_installed()) {
# estimate weights
weights <- calc_weight(x = x,
                                 z = z, 
                                 estimand = "ATE",
                                 method = "COT",
                                 options = list(lambda = 0))
#we can also use the dataSim object directly
weightsDS <- calc_weight(x = data,
                                 z = NULL,
                                 estimand = "ATE",
                                 method = "COT",
                                 options = list(lambda = 0))
all.equal(weights@w0, weightsDS@w0)
all.equal(weights@w1, weightsDS@w1)
}
}
\seealso{
\code{\link[=estimate_effect]{estimate_effect()}}
}
