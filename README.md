
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spbsampling <img src="man/figures/logo.png" align="right" />

[![CRAN
version](http://www.r-pkg.org/badges/version/Spbsampling)](https://cran.r-project.org/package=Spbsampling)
[![downloads](https://cranlogs.r-pkg.org/badges/Spbsampling)](https://cran.r-project.org/package=Spbsampling)

An [R](https://www.r-project.org) package for *spatially balanced
sampling*.

A sample is *spatially balanced* when is spread on the auxiliary space.
The **Spbsampling** package provides functions to draw this kind of
samples. It contains fast implementations (C++ via **Rcpp** and
**RcppArmadillo**) of the included sampling methods, and related
functions to deal with distance matrix standardization and spatial
balance index.

For details regarding the implemented sampling designs, look at the
references section.

Authors: Francesco Pantalone, Roberto Benedetti, Federica Piersimoni.

Maintainer: Francesco Pantalone.

## Installation

You can install the released version of Spbsampling from
[CRAN](https://CRAN.R-project.org)

``` r
install.packages("Spbsampling")
```

or the development version from GitHub

``` r
# using devtools
# install.packages("devtools")
devtools::install_github("FrancescoPantalone/Spbsampling")
# or using remotes
# install.packages("remotes")
remotes::install_github("FrancescoPantalone/Spbsampling")
```

## References

Benedetti R, Piersimoni F (2017). “A spatially balanced design with
probability function proportional to the within sample distance.”
*Biometrical Journal*, 59(5), 1067–1084.
<https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201600194>

Benedetti R, Piersimoni F (2017). “Fast Selection of Spatially Balanced
Samples.” *arXiv*. <https://arxiv.org/abs/1710.09116>

## License

GPL-3
