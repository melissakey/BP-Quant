#include <Rcpp.h> 
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
DataFrame compute_posterior(NumericVector pi_probs, NumericMatrix x_mat, int nu, double tol = 1e-8, int max_output = 65536) {
  double log_prior, tmp_post_prob, x;
  IntegerVector post_prob_keys(max_output);
  NumericVector post_prob(max_output, 0.0);
  // NumericVector::iterator pp =; 
  int k = 0;
  int n_combs = pow(2,nu);
  
  for(int n = 0; n < n_combs; n++) {
    tmp_post_prob = 0;
    for(int i = nu-1; i >= 0; i--) {
      x = fmod(n / pow(2, i), 2);
      log_prior = x * log(pi_probs[i]) +  log(1 - pi_probs[i])* (1 - x);
      tmp_post_prob = tmp_post_prob + log_prior + log(x_mat(i, x + 1));
    }
    if(exp(tmp_post_prob) > tol) {
      post_prob[k] = exp(tmp_post_prob);
      post_prob_keys[k] = n;
      k++;
    }
  }
  
  DataFrame df = DataFrame::create(
    Named("index") = head(post_prob_keys, k),
    Named("post_prob") = head(post_prob, k)
  );
  
  return df;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 nu <- 16
pi_probs <- c(.5, rep(.5/(nu-1), nu-1))
tmp <- runif(nu)
x_mat <- cbind(tmp, 1 - tmp)
tmp <- compute_posterior(pi_probs, x_mat, nu, tol = 1e-40)
*/
