<!-- badges: start -->
[![R-CMD-check](https://github.com/ericdunipace/causalOT/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/ericdunipace/causalOT/actions/workflows/check-standard.yaml)
<!-- badges: end -->

# causalOT: Optimal transport methods for causal inference

This R package implements the methods described in [*Optimal transport methods for causal inference*](https://arxiv.org/abs/2109.01991).

## Installation
This package can be installed in a few ways.

### 1. devtools
Using the `remotes` package in `R`, one can install the package with 
```
remotes::install_github("ericdunipace/causalOT")
```

### 2. downlaod and install
After downloading the git package using `git clone` or by downloading the .zip file from the button above (Code -> Download Zip) and unzipping, you can install the package with
```
devtools::install("path/to/causalOT")
```

## Additional required software
If using Causal Optimal Transport based on the Sinkhorn divergence, you will need to install Python 3.x and the `geomloss` package, which depends on `numpy`, `scipy`, `torch`, and `pykeops`.

## Useage
The functions in the package are built to construct weights to make distributions more same and estimate causal effects. The primary method we recommend is by using optimal transport weights which balance distributions by design. For more information about using this package, see the vignette "Using causalOT".

## Reproducing the paper
In the folder `inst/Reproduce` you can find code and an RMarkdown file to reproduce the figures present in the [paper](https://arxiv.org/abs/2109.01991).

## Package author
[Eric Dunipace](https://ericdunipace.github.io)

## License
This package is licensed under GPL 3.0.

