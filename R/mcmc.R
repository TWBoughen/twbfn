#' MCMC Sampling for PPA Model
#'
#' This function attempts to estimate the parameters of the preferential attachment function
#' from the sequence of degrees of a network given that it was generated using the model
#' from Rudas.
#' 
#' @param data A numeric matrix containing the data, with the first column being the unique degrees and the second being the
#' counts of those degrees.
#' @param init_n0 Initial value for the parameter \code{n0}.
#' @param init_a Initial value for the parameter \code{a}.
#' @param init_b Initial value for the parameter \code{b}.
#' @param init_lambda Initial value for the parameter \code{lambda}.
#' @param n_iter Number of iterations (default is 10000).s
#' @param a_sd Initial standard deviation for \code{a} proposals (default is 0.1).
#' @param b_sd Initial standard deviation for \code{b} proposals (default is 0.1).
#' @param n0_sd Initial standard deviation for \code{n0} proposals (default is 1).
#' @param adapt_interval Interval for adapting the proposal variance (default is 1000).
#' @param burn_in Period where MCMC is not adaptive. 
#'
#' @return A list containing the results.
#' 
#' @family model fitting functions
#' 
#' @export
ppa_mcmc <- function(data, init_n0, init_a, init_b, init_lambda, 
                     n_iter = 10000, a_sd = 0.1, b_sd = 0.1, n0_sd = 1, 
                     adapt_interval = 1000, burn_in=10000) {
  
  # Call the C++ function using .Call()
  result <- .Call('_twbfn_ppa_mcmc', data, init_n0, init_a, init_b, init_lambda, 
                  n_iter, a_sd, b_sd, n0_sd, adapt_interval, burn_in)
  
  return(result)
}






