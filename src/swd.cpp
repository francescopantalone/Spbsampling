#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Sum Within Distance (Spatially Balanced Sampling).
//'
//' This is an implemention of a spatially balanced design, with a probability function proportional to the within sample distance, using the sum of distance as an index of the within sample distance (Sum Within Distance, \code{swd} in short).
//' To have a constant inclusion probabilities \eqn{\pi_{i}=nsamp/N}, where \eqn{nsamp} is sample size and \eqn{N} is population size, standardize the distance matrix with function \code{\link{stsum}}.
//'
//' @param dis A distance matrix NxN that specifies how far are all the pairs of units in the population.
//' @param nsamp Sample size.
//' @param bexp Parameter \eqn{\beta} for the algorithm. The higher \eqn{\beta} is, the more the sample is going to be spread.
//' @param nrepl Number of samples to draw (default = 1).
//' @param niter Number of iterations for the algorithm. More iterations are better but require more time. Usually 10 is very efficient (default = 10).
//' @return Return a matrix 2 x \code{nrepl} with \code{nrepl} samples drawn. In particular, the element \eqn{a_{ij}}{a_ij} is the j-th unit of the population drawn in the i-th sample.
//' @references
//' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
//' @examples
//' # Example 1
//' # Draw 20 samples of dimension 15 without constant probabilities and beta = 1
//' dis <- as.matrix(dist(cbind(income_emilia$x_coord, income_emilia$y_coord))) # distance matrix
//' nsamp <- 15  # sample size
//' nrepl <- 20  # number of samples to draw
//' niter <- 10  # number of iterations in the algorithm
//' bexp <- 10   # parameter beta
//' samples <- swd(dis, niter, nsamp, nrepl, bexp)  # drawn samples
//' \donttest{
//' # Example 2
//' # Draw 20 samples of dimension 15 with constant probabilities equal to nsamp/N and beta = 10
//' # with N = population size
//' dis <- as.matrix(dist(cbind(income_emilia$x_coord,income_emilia$y_coord))) # distance matrix
//' nsamp <- 15  # sample size
//' nrepl <- 20  # numbers of samples to drawn
//' niter <- 10  # numbers of iterations in the algorithm
//' bexp <- 10  # parameter beta
//' vec <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stsum(dis, vec, 1e-15, 1000) # standardized matrix
//' samples <- swd(stand_dist, niter, nsamp, nrepl, bexp)  # drawn samples
//' }
//' @export
// [[Rcpp::export]]

arma::mat swd (arma::mat dis, int nsamp, int bexp, int nrepl = 1, int niter = 10)
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
  if(nsamp > npo)
  {
    throw Rcpp::exception("Sample size greater than population size.");
  }
  if(nsamp <= 0)
  {
    throw Rcpp::exception("Sample size negative or 0.");
  }
  if(nrepl < 0)
  {
    throw Rcpp::exception("bexp has to be different from 0.");
  }
  if(nrepl <= 0)
  {
    throw Rcpp::exception("nrepl has to be greater than 0.");
  }
  if(niter <= 0)
  {
    throw Rcpp::exception("niter has to be greater than 0");
  }
  arma::mat selez(nsamp * nrepl, 2);
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
    urand = arma::sort_index(rand);
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
        if (gcod(2 * npo + z - 1) < pow((totb / totc), bexp))
        {
          ch = codord(w - 1);
          codord(w - 1) = codord(k - 1);
          codord(k - 1) = ch;
        }
        z++;
      }
      iter++;
    }
    for (int i=((cc - 1) * nsamp); i < (cc * nsamp); i++)
    {
      selez(i, 0) = cc;
    }
    int j = 0;
    for (int i = ((cc - 1) * nsamp); i < (cc * nsamp); i++)
    {
      selez(i, 1) = codord(j);
      j++;
    }
  }
  return selez;
}
