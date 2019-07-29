#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Standardize a symmetric matrix (distances) to fixed row (column) products
//'
//' Standardizes a distance matrix to fixed rows and columns
//' products. The function iteratively constrains a logarithmic transformed
//' matrix to know products, and in order to keep the symmetry of the matrix,
//' at each iteration performs an average with its transpose. When the known
//' products are all equal to a constant (e.g. 0), this method provides a
//' simple and accurate way to scale a distance matrix to a doubly stochastic
//' matrix.
//'
//' The standardized matrix will not be affected by problems arising from units
//' with different inclusion probabilities caused by undesired features of the
//' spatial distribution of the population, as edge effects and/or isolated
//' points.
//'
//' @param mat A distance matrix size NxN.
//' @param vec A vector of row (column) constraints.
//' @param differ A scalar with the maximum accepted difference with the constraint (default = 1e-15).
//' @param niter An integer with the maximum number of iterations (default = 1000).
//' @return Returns a standardized distance matrix of size NxN.
//' @references
//' Benedetti R, Piersimoni F (2017). A spatially balanced design with
//' probability function proportional to the within sample distance.
//' \emph{Biometrical Journal}, \strong{59}(5), 1067-1084.
//' \url{https://doi.org/10.1002/bimj.201600194}
//' @examples
//' \dontshow{
//' d <- matrix(runif(200), 100, 2)
//' dis <- as.matrix(dist(d))
//' con <- rep(0, nrow(dis))
//' stand_dist <- stprod(mat = dis, vec = con)
//' }
//' \donttest{
//' dis <- as.matrix(dist(cbind(simul1$x, simul1$y))) # distance matrix
//' con <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(mat = dis, vec = con) # standardized matrix
//' }
//' @export
// [[Rcpp::export]]

arma::mat stprod(arma::mat mat, arma::vec vec, double differ = 1e-15, int niter = 1000)
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
  mat = log(mat);
  mat.diag().fill(0);
  int v = 1;
  double dif = 1;
  while(dif > differ && v < niter)
  {
    rowsums=sum(mat, 1);
    for(int i = 0; i < nsize; i++)
    {
      mat.row(i) -= (rowsums(i) / (nsize - 1));
      mat.row(i) += vec(i) / (nsize - 1);
    }
    mat = (mat + mat.t()) / 2;
    mat.diag().fill(0);
    dif = abs(sum(mat, 1) - vec).max();
    v++;
  }
  mat = exp(mat);
  mat.diag().fill(0);
  return mat;
}
