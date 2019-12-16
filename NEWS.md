# Spbsampling 1.3.1

* In the functions `pwd()` and `swd()` replaced `arma::sort_index()` with `arma::stable_sort_index()`.
* In the functions `pwd()`, `swd()` and `hpwd()` fixed the sample size check.
* In the functions `pwd()`, `swd()` and `hpwd()` changed the name of the parameter `bexp` in `beta`.
* In the functions `stprod()` and `stsum()` changed the name of the parameter `vec`in `con`.

# Spbsampling 1.3.0

## New features

* In the functions `pwd()` and `hpwd()` a new parameter called `bexp` has been introduced, which regulates the amount of spread of the samples. Default value set to `10`.
* In the function `swd()` the parameter `bexp` is now set to double and with default equal to `10`.

## Minor improvements

* Minor improvements of guides:
  * some updates to the examples;
  * removed some non-ASCII characters;
  * inserted links to reference websites.


# Spbsampling 1.2.0

* Changed output of functions `pwd()`, `swd()` and `hpwd()` from a matrix of dimension `(nsamp * nrepl) * 2` to a matrix of dimension `nrepl * nsamp`.
* Removed **Rdpack** from `Imports`.
* Minor improvements of guides.
* Changed name of file `heurprod.cpp` to `hpwd.cpp` 

# Spbsampling 1.1.0

## New features

* Added a new function `sbi()` for compute the spatial index of a sample, implemented in C++ through the Armadillo library, using **Rcpp** and **RcppArmadillo**.
* In the functions `pwd()` and `swd()` the parameter `nrepl` is now set to default at `1`.
* In the functions `pwd()` and `swd()` the parameter `niter` is now set to default at `10`.
* In the functions `stprod()` and `stsum()` the parameter `differ` is now set to default at `1e-15`.
* In the functions `stprod()` and `stsum()` the parameter `niter` is now set to default at `1000`.
* In the functions `stprod()` and `stsum()` there is no print for each convergence step anymore. 

## Improvements under the hood

* The functions `pwd()` and `swd()` now take advantages of the Armadillo library.
* The function `stprod()` is now implemented in C++ through the Armadillo library, using **Rcpp** and **RcppArmadillo**.
* The function `stsum()` is now implemented in C++ through the Armadillo library, using **Rcpp** and **RcppArmadillo**.
* The function `hpwd()` is now implemented in C++ through the Armadillo library, using **Rcpp** and **RcppArmadillo**.
* In all the functions there are now checks about the correct inputs.

## Minor improvements

* Correction of some typos along the guides.

## Miscellaneous

* Added a `NEWS.md` file to track changes to the package.
* Updated `DESCRIPTION` file.

# Spbsampling 1.0.0

* First release on CRAN.
