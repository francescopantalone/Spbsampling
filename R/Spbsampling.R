#' Spatially balanced sampling designs
#'
#' Selection of spatially balanced samples. In particular, the implemented
#' sampling designs allow to select probability samples well spread over the
#' population of interest, in any dimension and using any distance function
#' (e.g. Euclidean distance, Manhattan distance). The implementation has been
#' done in \code{C++} through the use of \code{Rcpp} and \code{RcppArmadillo}.
#' @author
#' Francesco Pantalone, Roberto Benedetti, Federica Piersimoni
#'
#' Maintainer: Francesco Pantalone \email{pantalone.fra@@gmail.com}
#' @references
#' Pantalone F, Benedetti R, Piersimoni F (2022). An R Package for
#' Spatially Balanced Sampling. \emph{Journal of Statistical Software, Code Snippets},
#' \strong{103}(2), 1-22.
#' <doi:10.18637/jss.v103.c02>
#'
#' Benedetti R, Piersimoni F (2017). A spatially balanced design with
#' probability function proportional to the within sample distance.
#' \emph{Biometrical Journal}, \strong{59}(5), 1067-1084.
#' \doi{10.1002/bimj.201600194}
#'
#' Benedetti R, Piersimoni F (2017). Fast Selection of Spatially Balanced Samples. \emph{arXiv}.
#' \url{https://arxiv.org/abs/1710.09116}
#' @docType package
#' @name Spbsampling
#'
#' @useDynLib Spbsampling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

