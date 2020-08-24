#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Product Within Distance (Spatially Balanced Sampling Design)
//'
//' Selects spatially balanced samples through the use of the
//' Product Within Distance design (PWD). To have constant inclusion
//' probabilities \eqn{\pi_{i}=n/N}, where \eqn{n} is sample size and
//' \eqn{N} is population size, the distance matrix has to be standardized with
//' function \code{\link{stprod}}.
//'
//' @param dis A distance matrix NxN that specifies how far all the pairs
//' of units in the population are.
//' @param n Sample size.
//' @param beta Parameter \eqn{\beta} for the algorithm. The higher
//' \eqn{\beta} is, the more the sample is going to be spread (default = 10).
//' @param nrepl Number of samples to draw (default = 1).
//' @param niter Maximum number of iterations for the algorithm. More iterations are
//' better but require more time. Usually 10 is very efficient (default = 10).
//' @return Returns a list with the following components:
//' \itemize{
//' \item \code{s}, a matrix \code{nrepl} x \code{n}, which contains the
//' \code{nrepl} selected samples, each of them stored in a row. In particular,
//' the i-th row contains all labels of units selected in the i-th sample.
//' \item \code{iterations}, number of iterations run by the algorithm.
//' }
//' @references
//' Benedetti R, Piersimoni F (2017). A spatially balanced design with
//' probability function proportional to the within sample distance.
//' \emph{Biometrical Journal}, \strong{59}(5), 1067-1084.
//' \url{https://doi.org/10.1002/bimj.201600194}
//' @examples
//' # Example 1
//' # Draw 1 sample of dimension 15 without constant inclusion probabilities
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' s <- pwd(dis = dis, n = 15)$s  # drawn sample
//' \donttest{
//' # Example 2
//' # Draw 1 sample of dimension 15 with constant inclusion probabilities
//' # equal to n/N, with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' con <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(mat = dis, con = con) # standardized matrix
//' s <- pwd(dis = stand_dist$mat, n = 15)$s  # drawn sample
//'
//' # Example 3
//' # Draw 2 samples of dimension 15 with constant inclusion probabilities
//' # equal to n/N, with N = population size, and an increased level of spread, beta = 20
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' con <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(mat = dis, con = con) # standardized matrix
//' s <- pwd(dis = stand_dist$mat, n = 15, beta = 20, nrepl = 2)$s  # drawn samples
//' }
//' @export
// [[Rcpp::export]]

Rcpp::List pwd (arma::mat dis, int n, double beta = 10, int nrepl = 1, int niter = 10)
{
  int npo = dis.n_rows;
  if(dis.is_square() == FALSE)
  {
    throw Rcpp::exception("The distance matrix has to be N x N.");
  }
  if(dis.is_symmetric() == FALSE)
  {
    Rcpp::warning("The distance matrix is not symmetric.");
  }
  if(n >= npo)
  {
    throw Rcpp::exception("Sample size equal or greater than population size.");
  }
  if(n <= 0)
  {
    throw Rcpp::exception("Sample size negative or 0.");
  }
  if(nrepl <= 0)
  {
    throw Rcpp::exception("nrepl has to be greater than 0.");
  }
  if(niter <= 0)
  {
    throw Rcpp::exception("niter has to be greater than 0");
  }
  arma::mat selez(nrepl, n);
  arma::vec ord(npo);
  arma::vec sor(npo);
  arma::vec codord(npo);
  arma::vec gcod(npo * 3);
  arma::vec rand(npo);
  arma::uvec urand(npo);
  arma::vec q(nrepl);
  arma::vec iter_vec(nrepl);
  double ch,totc,totb;
  ch = 0.0;
  totc = 0.0;
  totb = 0.0;
  int w, k, z, iter;
  dis = arma::pow(dis, beta);
  dis.diag().fill(1);
  dis = log(dis);
  arma::vec cod = arma::linspace(1, npo, npo);
  bool acc = true;
  for(int cc = 1; cc <= nrepl; cc++)
  {
    iter = 0;
    acc = true;
    rand = Rcpp::runif(npo);
    urand = arma::stable_sort_index(rand);
    codord = cod(urand);
    while(iter < niter && acc == true)
    {
      z = 1;
      gcod = Rcpp::runif(npo * 3);
      acc = false;
      while(z <= npo)
      {
        w = trunc(gcod(z - 1) * n) + 1;
        k = trunc(gcod(npo + z - 1) * (npo - n)) + 1 + n;
        totc = 0;
        for (int i = 0; i < n; i++)
        {
          totc = totc + dis(codord(w - 1) - 1, codord(i) - 1);
        }
        totb = 0;
        for (int i = 0; i < n; i++)
        {
          totb = totb + dis(codord(k - 1) - 1, codord(i) - 1);
        }
        totb = totb - dis(codord(k - 1) - 1, codord(w - 1) - 1);
        if (log(gcod(2 * npo + z - 1)) < (totb - totc))
        {
          ch = codord(w - 1);
          codord(w - 1) = codord(k - 1);
          codord(k - 1) = ch;
          acc = true;
        }
        z++;
      }
      iter++;
    }
    selez.row(cc - 1) = codord.subvec(0, n - 1).t();
    iter_vec(cc - 1) = iter;
  }
  return Rcpp::List::create(Rcpp::Named("s") = selez, Rcpp::Named("iterations") = iter_vec);
}
