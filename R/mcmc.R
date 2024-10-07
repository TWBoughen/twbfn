#' MCMC Sampling for PPA Model
#'
#' This function attempts to estimate the parameters of the preferential attachment function
#' from the sequence of degrees of a network given that it was generated using the model
#' from Rudas.
#' 
#' The model is somewhat adaptive in that the proposal variance is increased for all parameters
#' if the acceptance rate is higher than the target acceptance rate, and decreases when it is below.
#'
#' @param data A numeric matrix containing the data, with the first column being the unique degrees and the second being the
#' counts of those degrees.
#' @param init_n0 Initial value for the parameter \code{n0}.
#' @param init_a Initial value for the parameter \code{a}.
#' @param init_b Initial value for the parameter \code{b}.
#' @param init_lambda Initial value for the parameter \code{lambda}.
#' @param n_iter Number of iterations (default is 10000).
#' @param a_sd Initial standard deviation for \code{a} proposals (default is 0.1).
#' @param b_sd Initial standard deviation for \code{b} proposals (default is 0.1).
#' @param n0_sd Initial standard deviation for \code{n0} proposals (default is 1).
#' @param adapt_interval Interval for adapting the proposal variance (default is 1000).
#' @param target_accept_rate Target acceptance rate for the adaptive MCMC (default is 0.234).
#' @param adapt_factor Adaptation factor for the proposal variance (default is 1.05).
#'
#' @return A list containing the results.
#' 
#' @export
ppa_mcmc <- function(data, init_n0, init_a, init_b, init_lambda, 
                     n_iter = 10000, a_sd = 0.1, b_sd = 0.1, n0_sd = 1, 
                     adapt_interval = 1000, target_accept_rate = 0.234, 
                     adapt_factor = 1.05, burn_in=10000) {
  
  # Call the C++ function using .Call()
  result <- .Call('_twbfn_ppa_mcmc', data, init_n0, init_a, init_b, init_lambda, 
                  n_iter, a_sd, b_sd, n0_sd, adapt_interval, target_accept_rate, adapt_factor, burn_in)
  
  return(result)
}


#' @export
ppa_mle <- function(dat){
  return(optim(c(10,1,1),fn=ppa_llh,dat=dat))
}
#' @export
ppa_llh <- function(pars, dat){
  lambda = find_lambda_cpp(pars[1], pars[2], pars[3], 0)
  return(-llh_cpp(dat, lambda, pars[2], pars[3], pars[1]))
}





