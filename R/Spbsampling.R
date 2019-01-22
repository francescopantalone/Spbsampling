#' Spatially balanced sampling designs.
#'
#' Selection of spatially balanced samples.
#' In particular, the implemented sampling designs allow to select probability samples spread over the population of interest, in any dimension and using any distance function (e.g. Euclidean distance, Manhattan distance).
#' The implementation has been done in \code{C++} through the use of \code{Rcpp} and \code{RcppArmadillo}.
#' @author
#' Francesco Pantalone, Roberto Benedetti, Federica Piersimoni
#'
#' Maintainer: Francesco Pantalone \email{pantalone.fra@@gmail.com}
#' @references
#' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
#'
#' \insertRef{fast_selection}{Spbsampling}
#' @docType package
#' @name Spbsampling
#'
#' @useDynLib Spbsampling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

