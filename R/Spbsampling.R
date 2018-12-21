#' Spatially balanced sampling designs.
#'
#' The \code{Spbsampling} provides functions to draw spatially balanced samples. It contains fast implementations (\code{C++} via \code{Rcpp}) of the included sampling methods.
#' In particular, the algorithms to draw spatially balanced samples are \code{\link{pwd}} and \code{\link{swd}}. These methods use a probability distribution, proportional to the within sample distance.
#' Moreover, there is a function \code{\link{hpwd}}, which is a heuristic method to achieve approximated \code{\link{pwd}} samples.
#' Finally, there are two functions, \code{\link{stprod}} and \code{\link{stsum}}, useful to standardize distance matrixes in order to achieve fixed sample size using, respectively, the functions \code{\link{pwd}} and \code{\link{swd}}.
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

