#include <Rcpp.h>
using namespace Rcpp;

// g_cpp function: computes the transformed value of x
// [[Rcpp::export]]
NumericVector g_cpp(NumericVector x, double n0 = 1, double a = 1, double b = 1, double eps = 0) {
  int n = x.size();
  NumericVector result(n);
  result = ifelse(x<n0, pow(x,a) + eps, pow(std::max(0.0,n0),a) + b*(x-std::max(0.0,n0)) + eps);
  return result;
}

// lambda_obj_cpp: computes the objective function for a given l
// [[Rcpp::export]]
double lambda_obj_cpp(double l, double n0, double a, double b, double eps, int mx = 1000) {
  IntegerVector n_int = seq_len(mx);
  NumericVector n (n_int);
  NumericVector lgs = log(1 + l / g_cpp(n, n0, a, b, eps));
  NumericVector cumsum_lgs = cumsum(lgs);
  double result = sum(exp(-cumsum_lgs));
  return result;
}

// lambda_obj_adj_cpp: computes the adjusted objective function (1 - lambda_obj_cpp)
// [[Rcpp::export]]
double lambda_obj_adj_cpp(double l, double n0, double a, double b, double eps, int mx = 1000) {
  return 1 - lambda_obj_cpp(l, n0, a, b, eps, mx);
}

// find_lambda_root: helper function for root finding
double find_lambda_root(double n0, double a, double b, double eps, double l) {
  return lambda_obj_adj_cpp(l, n0, a, b, eps);
}

// find_lambda_cpp: implements root finding using the bisection method
// [[Rcpp::export]]
double find_lambda_cpp(double n0, double a, double b, double eps) {
  double lower = 0;  // lower bound for lambda
  double upper = 10; // upper bound for lambda
  double tol = 1e-7; // tolerance level for convergence
  int max_iter = 100; // maximum iterations
  double mid;
  
  // Bisection method for root-finding
  for (int iter = 0; iter < max_iter; ++iter) {
    mid = (lower + upper) / 2.0;
    double f_mid = find_lambda_root(n0, a, b, eps, mid);
    
    if (std::abs(f_mid) < tol) {
      break; // root found
    } else if (f_mid > 0) {
      upper = mid;
    } else {
      lower = mid;
    }
  }
  
  return mid;  // return the estimated root (lambda)
}

// dgpa_cpp: computes the density function for a given x
// [[Rcpp::export]]
double dgpa_cpp(int x, double n0, double a, double b, double eps, double lambda, bool logd = true) {
  double result;
  if (x == 1) {
    result = -std::log(1 + g_cpp(NumericVector::create(x), n0, a, b, eps)[0] / lambda);
  } else {
    double term1 = -std::log(1 + g_cpp(NumericVector::create(x), n0, a, b, eps)[0] / lambda);
    
    // Create sequence manually
    NumericVector seq_x(x - 1);
    for (int i = 0; i < x - 1; ++i) {
      seq_x[i] = i + 1;
    }
    
    double term2 = sum(log(1 + lambda / g_cpp(seq_x, n0, a, b, eps)));
    result = term1 - term2;
  }
  
  if (!logd) {
    result = exp(result);
  }
  return result;
}

// llh_cpp: computes the log-likelihood
// [[Rcpp::export]]
double llh_cpp(NumericMatrix data, double lambda, double a, double b, double n0, double eps = 0) {
  int n = data.nrow();
  double llh = 0.0;
  
  for (int i = 0; i < n; ++i) {
    int x = data(i, 0);
    int count = data(i, 1);
    llh += count * dgpa_cpp(x, n0, a, b, eps, lambda, true);
  }
  
  return llh;
}

// joint_prior_cpp: computes the joint prior probability
// [[Rcpp::export]]
double joint_prior_cpp(double a, double b, double n0) {
  double log_prior_a = dnorm(NumericVector::create(a), 0, 10, true)[0];
  double log_prior_b = dgamma(NumericVector::create(b), 1.0, 0.01, true)[0];
  double log_prior_n0 = dgamma(NumericVector::create(n0), 1.0, 0.001, true)[0];
  
  return log_prior_n0 + log_prior_a + log_prior_b;
}

// joint_post_cpp: computes the joint posterior probability
// [[Rcpp::export]]
double joint_post_cpp(NumericMatrix data, double lambda, double a, double b, double n0) {
  double llh = llh_cpp(data, lambda, a, b, n0);
  double prior = joint_prior_cpp(a, b, n0);
  return llh + prior;
}


// Helper function declarations
// g_cpp, lambda_obj_cpp, find_lambda_cpp, dgpa_cpp, llh_cpp, joint_prior_cpp, joint_post_cpp

// gpa_mcmc_cpp: MCMC sampling function with adaptive proposal variance
  
// [[Rcpp::export]]
List ppa_mcmc(NumericMatrix data, double init_n0, double init_a, double init_b, double init_lambda, 
                    int n_iter = 10000, double a_sd = 0.1, double b_sd = 0.1, double n0_sd = 1, 
                    int adapt_interval = 1000, double target_accept_rate = 0.234, 
                    double adapt_factor = 1.05) {
  
  // Initialize traces and variables
  NumericVector lambda_trace(n_iter);
  NumericVector a_trace(n_iter);
  NumericVector b_trace(n_iter);
  NumericVector n0_trace(n_iter);
  NumericVector post_trace(n_iter);
  
  double last_a = init_a;
  double last_b = init_b;
  double last_n0 = init_n0;
  double last_lambda = init_lambda;
  double last_post = joint_post_cpp(data, last_lambda, last_a, last_b, last_n0);
  
  int accepted = 0;  // Track accepted proposals
  
  // Run MCMC
  for (int i = 0; i < n_iter; ++i) {
    // Propose new parameters
    double new_a = rnorm(1,last_a, a_sd)[0];
    double new_b = rnorm(1,last_b, b_sd)[0];
    double new_n0 = rnorm(1,last_n0, n0_sd)[0];
    
    if (new_b < 0){
      a_trace[i] = last_a;
      b_trace[i] = last_b;
      n0_trace[i] = last_n0;
      lambda_trace[i] = last_lambda;
      post_trace[i] = last_post;
    } // Skip iteration if b is negative
    
    // Find new lambda using the bisection method
    double new_lambda = find_lambda_cpp(new_n0, new_a, new_b, 0);
    
    // Compute new posterior
    double new_post = joint_post_cpp(data, new_lambda, new_a, new_b, new_n0);
    
    // Accept/reject step
    double logA = new_post - last_post;
    if (std::log(R::runif(0, 1)) < logA) {
      // Accept the new state
      last_a = new_a;
      last_b = new_b;
      last_n0 = new_n0;
      last_lambda = new_lambda;
      last_post = new_post;
      accepted++;
    }
    
    // Store current state
    a_trace[i] = last_a;
    b_trace[i] = last_b;
    n0_trace[i] = last_n0;
    lambda_trace[i] = last_lambda;
    post_trace[i] = last_post;
    
    // Output progress every 1000 iterations
    if ((i + 1) % adapt_interval == 0) {
      Rcpp::Rcout << "Iteration " << i + 1 << "/" << n_iter << ": "
                  << "a = " << last_a << ", "
                  << "b = " << last_b << ", "
                  << "n0 = " << last_n0 << ", "
                  << "lambda = " << last_lambda << ", "
                  << "posterior = " << last_post << std::endl;
      
      // Adaptive proposal variance based on acceptance rate
      double accept_rate = static_cast<double>(accepted) / adapt_interval;
      Rcpp::Rcout << "Acceptance rate = " << accept_rate << std::endl;
      
      // Adjust proposal variances if acceptance rate is too far from target
      if (accept_rate < target_accept_rate) {
        a_sd /= adapt_factor;  // Reduce variance to increase acceptance
        b_sd /= adapt_factor;
        n0_sd /= adapt_factor;
      } else {
        a_sd *= adapt_factor;  // Increase variance to decrease acceptance
        b_sd *= adapt_factor;
        n0_sd *= adapt_factor;
      }
      
      // Reset acceptance count after each adaptation interval
      accepted = 0;
    }
  }
  
  return List::create(Named("a_trace") = a_trace,
                      Named("b_trace") = b_trace,
                      Named("n0_trace") = n0_trace,
                      Named("lambda_trace") = lambda_trace,
                      Named("post_trace") = post_trace,
                      Named("a_sd_final") = a_sd,
                      Named("b_sd_final") = b_sd,
                      Named("n0_sd_final") = n0_sd);
}


