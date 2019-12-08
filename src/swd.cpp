#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Sum Within Distance (Spatially Balanced Sampling Design)
//'
//' Selects spatially balanced samples through the use of the
//' Sum Within Distance design (SWD). To have a constant inclusion
//' probabilities \eqn{\pi_{i}=nsamp/N}, where \eqn{nsamp} is sample size and
//' \eqn{N} is population size, the distance matrix has to be standardized with
//' function \code{\link{stsum}}.
//'
//' @param dis A distance matrix NxN that specifies how far all the pairs
//' of units in the population are.
//' @param nsamp Sample size.
//' @param beta Parameter \eqn{\beta} for the algorithm. The higher
//' \eqn{\beta} is, the more the sample is going to be spread.
//' @param nrepl Number of samples to draw (default = 1).
//' @param niter Number of iterations for the algorithm. More iterations are
//' better but require more time. Usually 10 is very efficient (default = 10).
//' @return Returns a matrix \code{nrepl} x \code{nsamp}, which contains the
//' \code{nrepl} selected samples, each of them stored in a row. In particular,
//' the i-th row contains all labels of units selected in the i-th sample.
//' @references
//' Benedetti R, Piersimoni F (2017). A spatially balanced design with
//' probability function proportional to the within sample distance.
//' \emph{Biometrical Journal}, \strong{59}(5), 1067-1084.
//' \url{https://doi.org/10.1002/bimj.201600194}
//' @examples
//' # Example 1
//' # Draw 1 sample of dimension 15 without constant inclusion probabilities
//' dis <- as.matrix(dist(cbind(income_emilia$x_coord, income_emilia$y_coord))) # distance matrix
//' s <- swd(dis = dis, nsamp = 15)  # drawn sample
//' \donttest{
//' # Example 2
//' # Draw 1 sample of dimension 15 with constant inclusion probabilities
//' # equal to nsamp/N, with N = population size
//' dis <- as.matrix(dist(cbind(income_emilia$x_coord,income_emilia$y_coord))) # distance matrix
//' con <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stsum(mat = dis, con = vec) # standardized matrix
//' s <- swd(dis = stand_dist, nsamp = 15)  # drawn sample
//'
//' # Example 3
//' # Draw 2 samples of dimension 15 with constant inclusion probabilities
//' # equal to nsamp/N, with N = population size and an increased level of spread, i.e. beta = 20
//' dis <- as.matrix(dist(cbind(income_emilia$x_coord,income_emilia$y_coord))) # distance matrix
//' con <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stsum(mat = dis, con = vec) # standardized matrix
//' s <- swd(dis = stand_dist, nsamp = 15, beta = 20, nrepl = 2)  # drawn samples
//' }
//' @export
// [[Rcpp::export]]

arma::mat swd (arma::mat dis, int nsamp, double beta = 10, int nrepl = 1, int niter = 10)
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
  if(nsamp >= npo)
  {
    throw Rcpp::exception("Sample size equal or greater than population size.");
  }
  if(nsamp <= 0)
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
  arma::mat selez(nrepl, nsamp);
  arma::vec ord(npo);
  arma::vec codord(npo);
  arma::vec gcod(npo * 3);
  arma::vec rand(npo);
  arma::uvec urand(npo);
  double ch, totc, totb;
  ch = 0.0;
  totc = 0.0;
  totb = 0.0;
  int w, k, z, iter;
  arma::vec cod = arma::linspace(1, npo, npo);
  for(int cc = 1;cc <= nrepl; cc++)
  {
    iter = 0;
    rand = Rcpp::runif(npo);
    urand = arma::stable_sort_index(rand);
    codord = cod(urand);
    while(iter < niter)
    {
      z = 1;
      gcod = Rcpp::runif(npo * 3);
      while(z <= npo)
      {
        w = trunc(gcod(z - 1) * nsamp) + 1;
        k = trunc(gcod(npo + z - 1) * (npo - nsamp)) + 1 + nsamp;
        totc = 0;
        for (int i = 0; i < nsamp; i++)
        {
          totc = totc + dis(codord(w - 1) - 1, codord(i) - 1);
        }
        totb = 0;
        for (int i = 0; i < nsamp; i++)
        {
          totb = totb + dis(codord(k - 1) - 1, codord(i) - 1);
        }
        totb = totb - dis(codord(k - 1) - 1, codord(w - 1) - 1);
        if (gcod(2 * npo + z - 1) < pow((totb / totc), beta))
        {
          ch = codord(w - 1);
          codord(w - 1) = codord(k - 1);
          codord(k - 1) = ch;
        }
        z++;
      }
      iter++;
    }
    selez.row(cc - 1) = codord.subvec(0, nsamp - 1).t();
  }
  return selez;
}
