#' MCMC Sampling for GPA Model using C++
#'
#' This function wraps the C++ implementation of the MCMC sampler for the GPA model.
#' It performs MCMC sampling using the provided input data and parameters.
#'
#' @param data A numeric matrix containing the data.
#' @param init_n0 Initial value for the parameter \code{n0}.
#' @param init_a Initial value for the parameter \code{a}.
#' @param init_b Initial value for the parameter \code{b}.
#' @param init_lambda Initial value for the parameter \code{lambda}.
#' @param n_iter Number of MCMC iterations (default is 10000).
#' @param a_sd Standard deviation for \code{a} proposals (default is 0.1).
#' @param b_sd Standard deviation for \code{b} proposals (default is 0.1).
#' @param n0_sd Standard deviation for \code{n0} proposals (default is 1).
#' @param adapt_interval Interval for adapting the proposal variance (default is 1000).
#' @param target_accept_rate Target acceptance rate for the adaptive MCMC (default is 0.234).
#' @param adapt_factor Adaptation factor for the proposal variance (default is 1.05).
#'
#' @return A list containing the results of the MCMC simulation.
#' 
#' @export
ppa_mcmc <- function(data, init_n0, init_a, init_b, init_lambda, 
                     n_iter = 10000, a_sd = 0.1, b_sd = 0.1, n0_sd = 1, 
                     adapt_interval = 1000, target_accept_rate = 0.234, 
                     adapt_factor = 1.05) {
  
  # Call the C++ function using .Call()
  result <- .Call('_twbfn_ppa_mcmc', data, init_n0, init_a, init_b, init_lambda, 
                  n_iter, a_sd, b_sd, n0_sd, adapt_interval, target_accept_rate, adapt_factor)
  
  return(result)
}
