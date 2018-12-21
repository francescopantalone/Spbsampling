#include <Rcpp.h>
using namespace Rcpp;

//' Product Within Distance (Spatially Balanced Sampling).
//'
//' This is an implemention of a spatially balanced design, with a probability function proportional to the within sample distance, using the product of distance as an index of the within sample distance (Product Within Distance, \code{pwd} in short).
//' To have constant inclusion probabilities \eqn{\pi_{i}=nsamp/N}, where \eqn{nsamp} is sample size and \eqn{N} is population size, the distance matrix has to be standardized with function \code{\link{stprod}}.
//' Note that there is a parameter \eqn{\beta} which regulates the spread of the sample: The higher \eqn{\beta} is, the more the sample is going to be spread.
//' This parameter is regulated by the exponent of the distance matrix (D^1 -> \eqn{\beta = 1}, D^10 -> \eqn{\beta = 10}).
//'
//' @param dis A distance matrix NxN that specifies how far are all the pairs of units in the population.
//' @param nsamp Sample size.
//' @param nrepl Number of samples to draw (default = 1).
//' @param niter Number of iterations for the algorithm. More iterations are better but require more time. Usually 10 is very efficient (default = 10).
//' @return Return a matrix 2x\code{nrepl} with \code{nrepl} samples drawn. In particular, the element \eqn{a_{ij}}{a_ij} is the j-th unit of the population drawn in the i-th sample.
//' @references
//' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
//' @examples
//' # Example 1
//' # Draw 20 samples of dimension 15 without constant probabilities and with beta = 1
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' nsamp <- 15  # sample size
//' nrepl <- 20  # number of samples to draw
//' niter <- 10  # number of iterations in the algorithm
//' samples <- pwd(dis, niter, nsamp, nrepl)  # drawn samples
//' \donttest{
//' # Example 2
//' # Draw 20 samples of dimension 15 with constant probabilities equal to nsamp/N
//' # with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' nsamp <- 15  # sample size
//' nrepl <- 20  # number of samples to draw
//' niter <- 10  # number of iterations in the algorithm
//' vec <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(dis, vec ,1e-15 ,1000) # standardized matrix
//' samples <- pwd(stand_dist, niter, nsamp, nrepl)  # drawn samples
//'
//' # Example 3
//' # Draw 20 samples of dimension 15 with constant probabilities equal to nsamp/N and beta = 10
//' # with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' dis <- dis^10 # setting beta = 10
//' nsamp <- 15  # sample size
//' nrepl <- 20  # number of samples to draw
//' niter <- 10  # number of iterations in the algorithm
//' vec <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(dis, vec, 1e-15, 1000) # standardized matrix
//' samples <- pwd(stand_dist, niter, nsamp, nrepl)  # drawn samples
//' }
//' @export
// [[Rcpp::export]]

NumericMatrix pwd (NumericMatrix dis, int nsamp, int nrepl = 1, int niter = 10)
{
  NumericMatrix selez(nsamp * nrepl, 2);
  int npo=dis.nrow();
  NumericMatrix discopy(npo, npo);
  NumericVector rcopy(npo);
  NumericVector ord(npo);
  NumericVector sor(npo);
  NumericVector codord(npo);
  NumericVector gcod(npo * 3);
  NumericVector rand(npo);
  double ch,totc,totb;
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
  discopy = clone(dis);
  for(int i = 0; i < npo; i++)
  {
    for(int j = 0; j < npo; j++)
    {
      if(i == j)
      {
        discopy(i, j) = 0;
      }
      else
      {
        discopy(i, j) = log(dis(i, j));
      }
    }
  }
  NumericVector cod(npo);
  for(int i = 0; i < npo; i++)
  {
    cod(i) = i + 1;
  }
  for(int cc = 1; cc <= nrepl; cc++)
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
          totc = totc + discopy(codord(w - 1) - 1, codord(i) - 1);
        }
        totb = 0;
        for (int i = 0; i < nsamp; i++)
        {
          totb = totb + discopy(codord(k - 1) - 1, codord(i) - 1);
        }
        totb = totb - discopy(codord(k - 1) - 1, codord(w - 1) - 1);
        if (log(gcod(2 * npo + z - 1)) < (totb - totc))
        {
          ch = codord(w - 1);
          codord(w - 1) = codord(k - 1);
          codord(k - 1) = ch;
        }
        z++;
      }
      iter++;
    }
    for (int i = ((cc - 1) * nsamp); i < (cc * nsamp); i++)
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
