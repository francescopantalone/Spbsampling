#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Standardize a symmetric matrix (distances) to fixed row (column) totals
//'
//' Standardizes a distance matrix to fixed rows and columns
//' products. The function iteratively constrains the rows sums of the matrix
//' to know totals, and in order to keep the symmetry of the matrix, at each
//' iteration performs an average with its transpose. When the known totals are
//' all equal to a constant (e.g. 1), this method provides a simple and
//' accurate way to scale a distance matrix to a doubly stochastic matrix.
//'
//' The standardized matrix will not be affected by problems arising from units
//' with different inclusion probabilities caused by undesired features of the
//' spatial distribution of the population, as edge effects and/or isolated
//' points.
//'
//' @param  mat A distance matrix size NxN.
//' @param  vec A vector of row (column) constraints.
//' @param differ A scalar with the maximum accepted difference with the constraint (default = 1e-15).
//' @param niter An integer with the maximum number of iterations (default = 1000).
//' @return Return a standardized distance matrix of size NxN.
//' @references
//' Benedetti R, Piersimoni F (2017). “A spatially balanced design with
//' probability function proportional to the within sample distance.”
//' \emph{Biometrical Journal}, \strong{59}(5), 1067–1084.
//' @examples
//' dis <- as.matrix(dist(cbind(simul2$x, simul2$y))) # distance matrix
//' con <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stsum(dis, con) # standardized matrix
//' @export
// [[Rcpp::export]]

arma::mat stsum(arma::mat mat, arma::vec vec, double differ = 1e-15, int niter = 1000)
{
  if(mat.is_square() == FALSE)
  {
    throw Rcpp::exception("The distance matrix has to be N x N.");
  }
  if(mat.is_symmetric() == FALSE)
  {
    Rcpp::warning("The distance matrix is not symmetric.");
  }
  if(vec.n_elem != mat.n_rows)
  {
    throw Rcpp::exception("The dimension of vec has to be equal to the number of row/column of the distance matrix.");
  }
  if(all(vec) == FALSE)
  {
    throw Rcpp::exception("Vec cannot contain a zero.");
  }
  if(differ < 0)
  {
    throw Rcpp::exception("differ has to be greater or equal than 0.");
  }
  if(niter <= 0)
  {
    throw Rcpp::exception("niter has to be greater than 0.");
  }
  int nsize = mat.n_rows;
  arma::vec rowsums;
  int v = 1;
  double dif = 1;
  while(dif > differ && v < niter)
  {
    rowsums = sum(mat, 1);
    for(int i = 0; i < nsize; i++)
    {
      mat.row(i) /= rowsums(i);
      mat.row(i) *= vec(i);
    }
    mat = (mat + mat.t()) / 2;
    mat.diag().fill(0);
    dif = abs((sum(mat, 1) - vec) / vec).max();
    v++;
  }
  return mat;
}
