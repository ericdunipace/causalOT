# causalOT 1.1

## Breaking changes
* `vcov()` for `causalEffect` objects now computes the variance from the full
  plug-in influence function, including the cross terms that the previous
  calculation left out. The sandwich estimator is no longer used for
  outcome-model fits. Standard errors and confidence intervals will differ
  from earlier versions. `vcov(object, type = "components")` gives the
  previous three-component calculation.
* The default `opt.direction` for `cotOptions()` is now `"primal"`, using the
  LBFGS solver. The previous default was `opt.direction = "dual"`, which
  remains available. The dual solver default is still
  `dual.solver = "fixed_point_joint"`.
<!-- TODO: revisit the default solver once solver = "mirror" supports
     multiple lambdas and balance constraints. -->
* The default dual solver for `"COT"` is now `dual.solver = "fixed_point_joint"`.

## New features
* `cotOptions(opt.direction = "primal", solver = "mirror")` runs a mirror
  descent solver for debiased COT weights. It requires a single finite
  `lambda` and does not yet support balance constraints.
* New dual solvers `dual.solver = "fixed_point_joint"` and `"joint"`. Their
  tuning parameters (`fixed_point_eta`, `n_potential_steps`, etc.) can be
  passed through `cotOptions(...)`.
* `cotOptions()` gains `lambda.selector` and `lambda.selector.control` for
  choosing how lambda is selected from the fitted grid. The default,
  `"bootstrap_legacy"`, is the same bootstrap used in earlier versions.
  Leave-one-out and mass-bootstrap selectors are available as experimental
  alternatives.
* The lambda grid, selection scores, and selected value are stored in
  `causalWeights@info$metric`.
* `Measure()` gains a `unit.id` argument.
* `cotOptions()` now accepts `ot_niter` and `ot_tol` to control the inner
  Sinkhorn iterations.

## Minor improvements and fixes
* `cotOptions()` now warns about arguments in `...` that no solver uses.
  Previously, misspelled arguments were silently ignored.

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
