#include <Rcpp.h>
using namespace Rcpp;

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

NumericMatrix swd (NumericMatrix dis, int nsamp, int bexp, int nrepl = 1, int niter = 10)
{
  NumericMatrix selez(nsamp * nrepl, 2);
  int npo=dis.nrow();
  NumericVector rcopy(npo);
  NumericVector ord(npo);
  NumericVector sor(npo);
  NumericVector codord(npo);
  NumericVector gcod(npo * 3);
  NumericVector rand(npo);
  double ch, totc, totb;
  if(nsamp > npo)
  {
    throw Rcpp::exception("Sample size greater than population size.");
  }
  if(nrepl < 0)
  {
    throw Rcpp::exception("nrepl has to be greater than 0.");
  }
  if(niter < 0)
  {
    throw Rcpp::exception("niter has to be greater than 0");
  }
  ch = 0.0;
  totc = 0.0;
  totb = 0.0;
  int w, k, z, iter;
  NumericVector cod(npo);
  for(int i = 0; i < npo; i++)
  {
    cod(i) = i + 1;
  }
  for(int cc = 1;cc <= nrepl; cc++)
  {
    iter = 0;
    rand = runif(npo);
    rcopy = clone(rand);
    sor = rcopy.sort();
    for(int i = 0; i < npo; i++)
    {
      for(int j = 0; j < npo; j++)
      {
        if(sor(i) == rand(j))
        {
          ord(i) = j + 1;
          rand(j)=-1;
          j = npo;
        }
      }
    }
    for(int i = 0; i < npo; i++)
    {
      codord(i) = cod(ord(i) - 1);
    }
    while(iter < niter)
    {
      z = 1;
      gcod = runif(npo * 3);
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
