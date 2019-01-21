#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Heuristic Product Within Distance (Spatially Balanced Sampling).
//'
//' The function \code{hpwd} provides a fast selection of spatially balanced samples: it is a heuristic and a fast implemention of the algorithm \code{\link{pwd}}.
//' It generates samples approximately with the same probabilities of the \code{\link{pwd}} but with a significantly smaller number of steps.
//' In fact, this algorithm randomly selects a sample of size \eqn{n} exactly with \eqn{n} steps, updating at each step the selection probability of not-selected units, depending on their distance from the units, that were already selected in the previous steps.
//' To have constant inclusion probabilities \eqn{\pi_{i}=nsamp/N}, where \eqn{nsamp} is sample size and \eqn{N} is population size, the distance matrix has to be standardized with function \code{\link{stprod}}.
//' Note that there is a parameter \eqn{\beta} which regulates the spread of the sample: The higher \eqn{\beta} is, the more the sample is going to be spread.
//' This parameter is regulated by the exponent of the distance matrix (D^1 -> \eqn{\beta = 1}, D^10 -> \eqn{\beta = 10}).
//'
//' @param dis A distance matrix NxN that specifies how far are all the pairs of units in the population.
//' @param nsamp Sample size.
//' @param nrepl Number of samples to draw (default = 1).
//' @return Return a matrix 2x\code{nrepl} with \code{nrepl} samples drawn. In particular, the element \eqn{a_{ij}}{a_ij} is the j-th unit of the population drawn in the i-th sample.
//' @references
//' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
//'
//' \insertRef{fast_selection}{Spbsampling}
//' @examples
//' # Example 1
//' # Draw 50 samples of dimension 10 without constant probabilities and beta = 1
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' samples <- hpwd(dis, 10, 50) # drawn samples
//' \donttest{
//' # Example 2
//' # Draw 100 samples of dimension 15 with constant probabilities equal to nsamp/N and beta = 1
//' # with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' vec <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(dis, vec, 1e-15, 1000) # standardized matrix
//' samples <- hpwd(stand_dist, 15, 100) # drawn samples
//'
//' # Example 3
//' # Draw 100 samples of dimension 15 with constant probabilities equal to nsamp/N and beta = 10
//' # with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' dis <- dis^10 # setting beta = 10
//' vec <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(dis, vec, 1e-15, 1000) # standardized matrix
//' samples <- hpwd(stand_dist, 15, 100) # drawn samples
//' }
//' @export
// [[Rcpp::export]]

arma::mat hpwd(arma::mat dis, int nsamp, int nrepl=1)
{
  int npop = dis.n_rows;
  if(dis.is_square() == FALSE)
  {
    throw Rcpp::exception("The distance matrix has to be N x N.");
  }
  if(dis.is_symmetric() == FALSE)
  {
    Rcpp::warning("The distance matrix is not symmetric.");
  }
  if(nsamp > npop)
  {
    throw Rcpp::exception("Sample size greater than population size.");
  }
  if(nsamp <= 0)
  {
    throw Rcpp::exception("Sample size negative or 0.");
  }
  if(nrepl <= 0)
  {
    throw Rcpp::exception("nrepl has to be greater than 0.");
  }
  arma::mat selez(nsamp * nrepl, 2);
  arma::vec r(1);
  arma::vec c(npop);
  arma::vec psel(npop);
  double drawn;
  for(int cc = 1; cc <= nrepl; cc++)
  {
    psel.fill(1.0 / npop);
    for(int j = 1; j <= nsamp; j++)
    {
      selez((cc - 1) * nsamp + j - 1, 0) = cc;
      r = Rcpp::runif(1);
      c.fill(0);
      drawn = 0;
      c(0) = psel(0);
      if(c(0) > r(0))
      {
        drawn = 0;
      }
      else
      {
        for(int i = 1; i < npop; i++)
        {
          c(i) = c(i - 1) + psel(i);
          if(c(i) > r(0))
          {
            drawn = i;
            break;
          }
        }
      }
      selez((cc - 1) * nsamp + j - 1, 1) = drawn + 1;
      for(int i = 0; i < npop; i++)
      {
        psel(i) = psel(i) * dis(selez((cc - 1) * nsamp + j - 1, 1) - 1, i);
      }
      psel = psel / sum(psel);
    }
  }
  return selez;
}
