#include <Rcpp.h>
#include <random>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <chrono>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


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
  if(n0>data(n-1,0)){
    return -INFINITY;
  }
  return llh;
}

// joint_prior_cpp: computes the joint prior probability
// [[Rcpp::export]]
double joint_prior_cpp(double a, double b, double n0) {
  double log_prior_a = dnorm(NumericVector::create(a), 0, 10.0, true)[0];
  double log_prior_b = dgamma(NumericVector::create(b), 1.0, 100.00, true)[0];
  double log_prior_n0 = dgamma(NumericVector::create(n0), 1.0, 1000.00, true)[0];
  
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
                    int adapt_interval = 1000, int burn_in=10000) {
  
  // Initialize traces and variables
  NumericVector lambda_trace(n_iter + burn_in);
  NumericVector a_trace(n_iter + burn_in);
  NumericVector b_trace(n_iter + burn_in);
  NumericVector n0_trace(n_iter + burn_in);
  NumericVector post_trace(n_iter + burn_in);
  NumericVector a_sd_trace(n_iter + burn_in);
  NumericVector b_sd_trace(n_iter + burn_in);
  NumericVector n0_sd_trace(n_iter + burn_in);
  
  double last_a = init_a;
  double last_b = init_b;
  double last_n0 = init_n0;
  double last_lambda = init_lambda;
  double last_post = joint_post_cpp(data, last_lambda, last_a, last_b, last_n0);
  double up_count =0;
  
  int accepted = 0;  // Track accepted proposals
  
  // Run MCMC
  for (int i = 0; i < (n_iter + burn_in); ++i) {
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
    a_sd_trace[i] = a_sd;
    b_sd_trace[i] = b_sd;
    n0_sd_trace[i] = n0_sd;
    
    
    // Output progress every 1000 iterations
    if ((i + 1) % adapt_interval == 0) {
      up_count++;
      Rcpp::Rcout << "Iteration " << i + 1 << "/" << n_iter+burn_in << ": "<<std::endl
                  << "a = " << last_a << ", "
                  << "b = " << last_b << ", "
                  << "n0 = " << last_n0 << ", "
                  << "lambda = " << last_lambda << ", "
                  << "posterior = " << last_post << std::endl
                  << "a_sd = "<< a_sd << ","
                  << "b_sd = "<< b_sd << ","
                  << "n0_sd = "<< n0_sd << std::endl;
      
      // Adaptive proposal variance based on acceptance rate
      double accept_rate = static_cast<double>(accepted) / adapt_interval;
      Rcpp::Rcout << "Acceptance rate = " << accept_rate << std::endl;
    if (i>(burn_in + adapt_interval)){
      a_sd = 2.38*sqrt(var(Rcpp::unique(a_trace[Rcpp::Range(0,up_count*adapt_interval)]))/3) + 0.1/sqrt(3);
      b_sd = 2.38*sqrt(var(Rcpp::unique(b_trace[Rcpp::Range(0,up_count*adapt_interval)]))/3) + 0.1/sqrt(3);
      n0_sd = 2.38*sqrt(var(Rcpp::unique(n0_trace[Rcpp::Range(0,up_count*adapt_interval)]))/3) + 0.1/sqrt(3);
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
                      Named("n0_sd_final") = n0_sd,
                      Named("a_sd") = a_sd_trace,
                      Named("b_sd") = b_sd_trace,
                      Named("n0_sd") = n0_sd_trace,
                      Named("data") = data);
                      
}

std::mt19937 rng(std::random_device{}());
// [[Rcpp::export]]
double tgamma_upper(double a, double z){
  return(
    boost::math::tgamma(a) - boost::math::tgamma_lower(a,z)
  );
}


// ggp_tail: computes tail of generalised gamma pareto process
// [[Rcpp::export]]
double ggp_tail(double x, double sigma, double tau,double beta){
  double mult = pow(beta,sigma*tau) / (sigma*tau*boost::math::tgamma(1-sigma));
  double p1 = pow(beta*x,-sigma*tau)*boost::math::tgamma_lower(sigma*(tau-1),beta*x);
  double p2 = -(tgamma_upper(1-sigma,beta*x) - pow(beta*x,-sigma)*exp(-beta*x))/sigma;
  
  return(
    mult*(p1+p2)
  );
}

// [[Rcpp::export]]
double ggp_tail_invgrad(double x, double sigma, double tau,double beta){
  double top = pow(beta, -sigma*tau)*(tgamma_upper(1-sigma,beta*x) - pow(beta*x,-sigma)*exp(-beta*x));
  double bottom = sigma*boost::math::tgamma_lower(sigma*(tau-1),beta*x);
  return(
    (1-top/bottom)/(sigma*tau)
  );
}
// [[Rcpp::export]]
double ggp_rho(double x, double sigma, double tau, double beta){
  return(
    pow(x,-1-sigma*tau)* boost::math::tgamma_lower(sigma*(tau-1),beta*x)/boost::math::tgamma(1-sigma)
  );
}

// [[Rcpp::export]]
double stable_rho(double x, double sigma){
  return(
    sigma*pow(x,-1-sigma)/boost::math::tgamma(1-sigma)
  );
}

// [[Rcpp::export]]
double stable_tail(double x, double sigma){
  return(
    pow(x, -sigma)/boost::math::tgamma(1-sigma)
  );
}



// // [[Rcpp::export]]
// double ggp_mu_integrand(double v, double w, double sigma, double tau, double beta){
//   return(
//     -expm1(-2*w*v)*ggp_rho(v, sigma, tau, beta)
//   );
// }
// 
// // [[Rcpp::export]]
// double ggp_mu(double w, double sigma, double tau, double beta) {
//   boost::math::quadrature::exp_sinh<double> integrator;
//   double error;
//   double L1;
//   double result = integrator.integrate(
//     [&](double v) {
//       return ggp_mu_integrand(v, w, sigma, tau, beta);
//    },
//     0.0,
//     &error,
//     &L1
//   );
//   return result;
// }


// [[Rcpp::export]]
double dggp_rho(double x, double sigma, double tau, double beta, double p){
  return(
    (p*pow(x*beta, -sigma*tau)*boost::math::tgamma_lower(sigma*(tau-1) + 1, beta*x)+
      (1-p)*pow(beta, -sigma)*pow(x, -1-sigma)*exp(-beta*x))/boost::math::tgamma(1-sigma)
  );
}

// [[Rcpp::export]]
double dggp_tail(double x, double sigma, double tau, double beta, double p){
  double rho0 = p/(sigma*tau)*
    (boost::math::tgamma(1-sigma)+  pow(beta*x, -sigma *tau) * boost::math::tgamma_lower(sigma*(tau-1)+1, beta*x)+tgamma_upper(1-sigma, beta*x));
  double rho1 = (1-p)/(sigma)*(-tgamma_upper(1-sigma, beta*x) + pow(beta*x, -sigma)*exp(-beta*x));
  return(
    (rho1+rho0)/boost::math::tgamma(1-sigma)
  );
}


// [[Rcpp::export]]
double gg_tail(double x, double sigma, double beta){
  
  double mult = pow(beta,sigma)/(boost::math::tgamma(1-sigma)*sigma);
  double p1 = pow(beta*x,-sigma)*exp(-beta*x);
  double p2 = tgamma_upper(1-sigma,beta*x);
  return(
    mult*(p1-p2)
  );
}
// [[Rcpp::export]]
double gg_rho(double w, double sigma, double beta){
  return(pow(w, -1-sigma)*exp(-beta*w) / boost::math::tgamma(1-sigma));
}



// // [[Rcpp::export]]
// double rho0(double w, double sigma, double c, double eta){
//   return(
//     eta*pow(w,-1-sigma)*exp(-c*w) / ( pow(c,sigma) * gamma(1-sigma))
//   );
// }
// 
// // [[Rcpp::export]]
// double rho1(double w, double sigma, double tau, double c , double eta){
//   return(
//     eta*pow(w,-1-tau) * boost::math::tgamma_lower(tau - sigma + 1, c*w) / ( pow(c,tau) * gamma(1-sigma))
//   );
// }

// [[Rcpp::export]]
double rho0hat(double w, double sigma, double c, double eta){
  w = std::max(w,pow(10,-10));
  return(
    (eta/(sigma*boost::math::tgamma(1-sigma)))*
      (pow(c,-sigma)*pow(w, -sigma)*exp(-c*w) - (boost::math::tgamma(1-sigma) - boost::math::tgamma_lower(1-sigma, c*w)))
  );
}


// [[Rcpp::export]]
double rho1hat(double w, double sigma, double tau, double c , double eta){
  w = std::max(w,pow(10,-10));
  double mult = eta/(tau*boost::math::tgamma(1-sigma));
  double t1 = boost::math::tgamma(1-sigma) - boost::math::tgamma_lower(1-sigma, c*w);
  double t2 = pow(c,-tau)*pow(w,-tau)*boost::math::tgamma_lower(tau-sigma+1,c*w);
  return(mult*(t1+t2));
}

// // [[Rcpp::export]]
// NumericMatrix sample_edges(NumericVector ws,double scaling, bool quiet=true){
//   NumericMatrix edges(pow(ws.size(),2),2); 
//   int counter=0;
//   for(int i=0;i<ws.size();i++){
//     if(!quiet){
//       std::cout<<i<<std::endl;
//     }
//     for(int j=0;j<i;j++){
//       double W = 1-exp(-2*scaling*ws[i]*ws[j]);
//       if(W>Rcpp::runif(1)[0]){
//           edges(counter,0) = i+1;
//           edges(counter,1) = j+1;
//           counter=counter+1;
//       }
//     }
//   }
//   
//   return(edges(Range(0,counter-1),_));
//   
// }

// [[Rcpp::export]]
NumericMatrix sample_edges(NumericVector ws,
                           double scaling,
                           bool quiet = true) {
  const int n = ws.size();
  std::vector<int> from, to;
  
  for (int i = 0; i < n; ++i) {
    if (!quiet) Rcpp::Rcout << i << '\n';
    
    // One random-number generation call per row.
    NumericVector u = Rcpp::runif(i);
    const double wi = ws[i];
    
    for (int j = 0; j < i; ++j) {
      const double W =
        -std::expm1(-2.0 * scaling * wi * ws[j]);
        
        if (u[j] < W) {
          from.push_back(i + 1);
          to.push_back(j + 1);
        }
    }
  }
  
  const R_xlen_t m = from.size();
  NumericMatrix edges(m, 2);
  
  for (R_xlen_t k = 0; k < m; ++k) {
    edges(k, 0) = from[k];
    edges(k, 1) = to[k];
  }
  
  return edges;
}
