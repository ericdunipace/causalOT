% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cotOOP.R
\name{OTProblem}
\alias{OTProblem}
\title{Object Oriented OT Problem}
\usage{
OTProblem(measure_1, measure_2, ...)
}
\arguments{
\item{measure_1}{An object of class \link{Measure}}

\item{measure_2}{An object of class \link{Measure}}

\item{...}{Not used at this time}
}
\value{
An R6 object of class "OTProblem"
}
\description{
Object Oriented OT Problem
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{device}}{the \code{\link[torch:torch_device]{torch::torch_device()}} of the data.}
\item{\code{dtype}}{the \link[torch:torch_dtype]{torch::torch_dtype} of the data.}
\item{\code{selected_delta}}{the delta value selected after \code{choose_hyperparameters}}
\item{\code{selected_lambda}}{the lambda value selected after \code{choose_hyperparameters}}
}
\if{html}{\out{</div>}}
}

\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{loss}}{prints the current value of the objective. Only availble after the \href{#method-OTProblem-solve}{\code{OTProblem$solve()}} method has been run}
\item{\code{penalty}}{Returns a list of the lambda and delta penalities that will be iterated through. To set these values, use the \href{#method-OTProblem-setup_arguments}{\code{OTProblem$setup_arguments()}} function.}
}
\if{html}{\out{</div>}}
}

\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-OTProblem-add}{\code{OTProblem$add()}}
\item \href{#method-OTProblem-subtract}{\code{OTProblem$subtract()}}
\item \href{#method-OTProblem-multiply}{\code{OTProblem$multiply()}}
\item \href{#method-OTProblem-divide}{\code{OTProblem$divide()}}
\item \href{#method-OTProblem-setup_arguments}{\code{OTProblem$setup_arguments()}}
\item \href{#method-OTProblem-solve}{\code{OTProblem$solve()}}
\item \href{#method-OTProblem-choose_hyperparameters}{\code{OTProblem$choose_hyperparameters()}}
\item \href{#method-OTProblem-info}{\code{OTProblem$info()}}
\item \href{#method-OTProblem-clone}{\code{OTProblem$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-add"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-add}{}}}
\subsection{Method \code{add()}}{
adds \code{o2} to the OTProblem
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$add(o2)}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{o2}}{A number or object of class OTProblem}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-subtract"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-subtract}{}}}
\subsection{Method \code{subtract()}}{
subtracts \code{o2} from OTProblem
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$subtract(o2)}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{o2}}{A number or object of class OTProblem}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-multiply"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-multiply}{}}}
\subsection{Method \code{multiply()}}{
multiplies OTProblem by \code{o2}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$multiply(o2)}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{o2}}{A number or an object of class OTProblem}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-divide"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-divide}{}}}
\subsection{Method \code{divide()}}{
divides OTProblem by \code{o2}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$divide(o2)}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{o2}}{A number or object of class OTProblem}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-setup_arguments"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-setup_arguments}{}}}
\subsection{Method \code{setup_arguments()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$setup_arguments(
      lambda,
      delta,
      grid.length = 7L,
      cost.function = NULL,
      p = 2,
      cost.online = "auto",
      debias = TRUE,
      diameter = NULL,
      ot_niter = 1000L,
      ot_tol = 0.001
    )}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lambda}}{The penalty parameters to try for the OT problems. If not provided, function will select some}
\item{\code{delta}}{The constraint paramters to try for the balance function problems, if any}
\item{\code{grid.length}}{The number of hyperparameters to try if not provided}
\item{\code{cost.function}}{The cost function for the data. Can be any function that takes arguments \code{x1}, \code{x2}, \code{p}. Defaults to the Euclidean distance}
\item{\code{p}}{The power to raise the cost matrix by. Default is 2}
\item{\code{cost.online}}{Should online costs be used? Default is "auto" but "tensorized" stores the cost matrix in memory while "online" will calculate it on the fly.}
\item{\code{debias}}{Should debiased OT problems be used? Defaults to TRUE}
\item{\code{diameter}}{Diameter of the cost function.}
\item{\code{ot_niter}}{Number of iterations to run the OT problems}
\item{\code{ot_tol}}{The tolerance for convergence of the OT problems}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
NULL
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ ot$setup_arguments(lambda = c(1000,10))
    }
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-solve"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-solve}{}}}
\subsection{Method \code{solve()}}{
Solve the OTProblem at each parameter value. Must run setup_arguments first.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$solve(
      niter = 1000L,
      tol = 1e-05,
      optimizer = c("torch", "frank-wolfe"),
      torch_optim = torch::optim_lbfgs,
      torch_scheduler = torch::lr_reduce_on_plateau,
      torch_args = NULL,
      osqp_args = NULL,
      quick.balance.function = TRUE
    )}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{niter}}{The nubmer of iterations to run solver at each combination of hyperparameter values}
\item{\code{tol}}{The tolerance for convergence}
\item{\code{optimizer}}{The optimizer to use. One of "torch" or "frank-wolfe"}
\item{\code{torch_optim}}{The \code{torch_optimizer} to use. Default is \link[torch:optim_lbfgs]{torch::optim_lbfgs}}
\item{\code{torch_scheduler}}{The \link[torch:lr_scheduler]{torch::lr_scheduler} to use. Default is \link[torch:lr_reduce_on_plateau]{torch::lr_reduce_on_plateau}}
\item{\code{torch_args}}{Arguments passed to the torch optimizer and scheduler}
\item{\code{osqp_args}}{Arguments passed to \code{\link[osqp:osqpSettings]{osqp::osqpSettings()}} if appropriate}
\item{\code{quick.balance.function}}{Should \code{\link[osqp:osqp]{osqp::osqp()}} be used to select balance function constraints (delta) or not. Default true.}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ ot$solve(niter = 1, torch_optim = torch::optim_rmsprop)
    }
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-choose_hyperparameters"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-choose_hyperparameters}{}}}
\subsection{Method \code{choose_hyperparameters()}}{
Selects the hyperparameter values through a bootstrap algorithm
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$choose_hyperparameters(
      n_boot_lambda = 100L,
      n_boot_delta = 1000L,
      lambda_bootstrap = Inf
    )}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_boot_lambda}}{The number of bootstrap iterations to run when selecting lambda}
\item{\code{n_boot_delta}}{The number of bootstrap iterations to run when selecting delta}
\item{\code{lambda_bootstrap}}{The penalty parameter to use when selecting lambda. Higher numbers run faster.}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ ot$choose_hyperparameters(n_boot_lambda = 10, 
                                             n_boot_delta = 10, 
                                             lambda_bootstrap = Inf)
    }
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-info"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-info}{}}}
\subsection{Method \code{info()}}{
Provides diagnostics after solve and choose_hyperparameter methods have been run.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$info()}\if{html}{\out{</div>}}
}
\subsection{Returns}{
a list with slots
\itemize{
\item \code{loss} the final loss values
\item \code{iterations} The number of iterations run for each combination of parameters
\item \code{balance.function.differences} The final differences in the balance functions
\item \code{hyperparam.metrics} A list of the bootstrap evalustion for delta and lambda values}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ ot$info()
    }
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OTProblem-clone"></a>}}
\if{latex}{\out{\hypertarget{method-OTProblem-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OTProblem$clone(deep = FALSE)}\if{html}{\out{</div>}}
}
\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

\examples{
## ------------------------------------------------
## Method `OTProblem(measure_1, measure_2)`
## ------------------------------------------------

if (torch::torch_is_installed()) {
  # setup measures
  x <- matrix(1, 100, 10)
  m1 <- Measure(x = x)
  
  y <- matrix(2, 100, 10)
  m2 <- Measure(x = y, adapt = "weights")
  
  z <- matrix(3,102, 10)
  m3 <- Measure(x = z)
  
  # setup OT problems
  ot1 <- OTProblem(m1, m2)
  ot2 <- OTProblem(m3, m2)
  ot <- 0.5 * ot1 + 0.5 * ot2
  print(ot)

## ------------------------------------------------
## Method `OTProblem$setup_arguments`
## ------------------------------------------------

  ot$setup_arguments(lambda = 1000)

## ------------------------------------------------
## Method `OTProblem$solve`
## ------------------------------------------------

  ot$solve(niter = 1, torch_optim = torch::optim_rmsprop)

## ------------------------------------------------
## Method `OTProblem$choose_hyperparameters`
## ------------------------------------------------

  ot$choose_hyperparameters(n_boot_lambda = 1,
                            n_boot_delta = 1, 
                            lambda_bootstrap = Inf)

## ------------------------------------------------
## Method `OTProblem$info`
## ------------------------------------------------

ot$info()
}
}
