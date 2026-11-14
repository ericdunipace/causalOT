# causalOT 1.1

## Breaking changes
* `vcov()` for `causalEffect` objects now computes the variance from the full
  plug-in influence function. Standard errors and confidence intervals will differ
  from earlier versions. `vcov(object, type = "components")` gives the
  previous calculation (which ignores cross terms)
* `cotOptions(opt.direction)` is deprecated. Use `cotOptions(solver)` instead.
* `cotOptions()` now defaults to `solver = "fixed_point"` and `niter = 500`
  for all solvers. 

## New features
* Updated convergence criterion to prevent needless long runs
* `cotOptions(polish = "auto")` refines COT weights with a monotone Newton
  polish using an exact Hessian. It is on by default and skips fits with
  more than 2,000 source points. 
* New solvers: `cotOptions(solver = "mirror")` runs a mirror 
  descent solver for debiased COT weights and `"fixed_point"` runs fast fixed
  point iterations also for debiased COT weights. 
* Lambda selection is now done via an RBF kernel. You can have the software
  choose thiswith `lambda = "auto"` before the estimation of weights or after 
  with weight estimation with option `lambda.selector = "auto"` in `cotOptions()`.
* The lambda grid, selection scores, and selected value are stored in
  `causalWeights@info$metric`.
* `Measure()` gains a `unit.id` argument to identify observations.
* `cotOptions()` now accepts `ot_niter` and `ot_tol` to control the inner
  Sinkhorn iterations.

## Minor improvements and fixes
* Bootstrap scores normalize resampled weights in `eval_lambda`, `eval_delta`,
  and `sbw_oop_bs_`; lambda scoring preserves the resampled candidate weights.
* `cotOptions()` now warns about arguments in `...` that no solver uses.
  Previously, misspelled arguments were silently ignored.

## Miscellaneous
* Updated maintainer e-mail address

# causalOT 1.0.4

## Minor improvements and fixes
* Fixed test that gave error on some linux builds

# causalOT 1.0.3

## Minor improvements and fixes
* Updated `R` required version to 4.1.0
* Fixed issue with html tags in the manual causing a CRAN note

# causalOT 1.0.2

## Minor improvements and fixes
* `rkeops` was removed to CRAN. Constructed a work-around to prevent having to remove all optional code that uses `rkeops` by moving `rkeops` build to an additional repository.

# causalOT 1.0.1

## Minor improvements and fixes
* Bug squashed under the hood arising from an update in an imported package

# causalOT 1.0

## Breaking changes
* Simplified the arguments to the main `calc_weight()` function
* Removed direct Python dependency
* Changed estimation function `estimate_effect()` arguments

## New features
* Added a new object-oriented set of functions for advanced users. Described in the new vignette, "Object Oriented COT Objects."
* Added `torch` based data structures allowing GPU support
* Added support for `rkeops` methods to allow for GPU use without the large memory footprint (CRAN version only for Unix machines)
* Added quadratic programming based on the `osqp` package where relevant
* Added internal optimal transport optimization removing reliance on `approxOT` package
* Added `vcov`, `confint`, and `coef` methods for the output of `estimate_effect()`
* Added a barycentric projection outcome estimation method via the function `barycentric_projection()` and a corresponding predict method.
* Added a `plot` method for the `causalWeights` and `causalEffects` objects

## Minor improvements and fixes
* Added a `NEWS.md` file to track changes to the package.
