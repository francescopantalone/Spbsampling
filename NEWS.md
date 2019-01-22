# Spbsampling 1.1.0

## New features

* Added a new function `sbi()` for compute the spatial index of a sample, implemented in C++ through the Armadillo library, using **Rcpp** and **RcppArmadillo**.
* In the function `pwd()` and `swd()` the parameter `nrepl` is now set to default at `1`.
* In the function `pwd()` and `swd()` the parameter `niter` is now set to default at `10`.
* In the function `stprod()` and `stsum()` the parameter `differ` is now set to default at `1e-15`.
* In the function `stprod()` and `stsum()` the parameter `niter` is now set to default at `1000`.
* In the function `stprod()` and `stsum()` there is no print for each convergence step anymore. 

## Improvement under the hood

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
