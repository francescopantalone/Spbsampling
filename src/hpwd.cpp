#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Heuristic Product Within Distance (Spatially Balanced Sampling Design)
//'
//' Selects spatially balanced samples through the use of
//' Heuristic Product Within Distance design (HPWD). To have constant inclusion
//' probabilities \eqn{\pi_{i}=n/N}, where \eqn{n} is sample size
//' and \eqn{N} is population size, the distance matrix has to be standardized
//' with function \code{\link{stprod}}.
//'
//' The HPWD design generates samples approximately with the same
//' probabilities of the \code{\link{pwd}} but with a significantly smaller
//' number of steps. In fact, this algorithm randomly selects a sample of size
//' \eqn{n} exactly with \eqn{n} steps, updating at each step the selection
//' probability of not-selected units, depending on their distance from the
//' units that were already selected in the previous steps.
//'
//' @param dis A distance matrix NxN that specifies how far all the pairs
//' of units in the population are.
//' @param n Sample size.
//' @param beta Parameter \eqn{\beta} for the algorithm. The higher
//' \eqn{\beta} is, the more the sample is going to be spread (default = 10).
//' @param nrepl Number of samples to draw (default = 1).
//' @return Returns a matrix \code{nrepl} x \code{n}, which contains the
//' \code{nrepl} selected samples, each of them stored in a row. In particular,
//' the i-th row contains all labels of units selected in the i-th sample.
//' @references
//' Benedetti R, Piersimoni F (2017). A spatially balanced design with
//' probability function proportional to the within sample distance.
//' \emph{Biometrical Journal}, \strong{59}(5), 1067-1084.
//' \url{https://doi.org/10.1002/bimj.201600194}
//'
//' Benedetti R, Piersimoni F (2017). Fast Selection of Spatially Balanced Samples. \emph{arXiv}.
//' \url{https://arxiv.org/abs/1710.09116}
//' @examples
//' # Example 1
//' # Draw 1 sample of dimension 10 without constant inclusion probabilities
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' s <- hpwd(dis = dis, n = 10) # drawn sample
//' \donttest{
//' # Example 2
//' # Draw 1 sample of dimension 15 with constant inclusion probabilities
//' # equal to n/N, with N = population size
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' con <- rep(1, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(mat = dis, con = con) # standardized matrix
//' s <- hpwd(dis = stand_dist$mat, n = 15) # drawn sample
//'
//' # Example 3
//' # Draw 2 samples of dimension 15 with constant inclusion probabilities
//' # equal to n/N, with N = population size, and an increased level of spread, beta = 20
//' dis <- as.matrix(dist(cbind(lucas_abruzzo$x, lucas_abruzzo$y))) # distance matrix
//' con <- rep(0, nrow(dis)) # vector of constraints
//' stand_dist <- stprod(mat = dis, con = con) # standardized matrix
//' s <- hpwd(dis = stand_dist$mat, n = 15, beta = 20, nrepl = 2) # drawn samples
//' }
//' @export
// [[Rcpp::export]]

arma::mat hpwd(arma::mat dis, int n, double beta = 10, int nrepl=1)
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
  if(n >= npop)
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
  arma::mat selez(nrepl, n);
  arma::vec r(1);
  arma::vec c(npop);
  arma::vec psel(npop);
  double drawn;
  dis = arma::pow(dis, beta);
  dis.diag().fill(0);
  for(int cc = 1; cc <= nrepl; cc++)
  {
    psel.fill(1.0 / npop);
    for(int j = 1; j <= n; j++)
    {
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
      selez(cc -1, j - 1) = drawn + 1;
      for(int i = 0; i < npop; i++)
      {
        psel(i) = psel(i) * dis(drawn, i);
      }
      psel = psel / sum(psel);
    }
  }
  return selez;
}
