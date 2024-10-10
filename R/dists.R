#' PMF of PPA model
#' 
#' Returns log of PMF
#' 
#' @param x Numeric or numeric vector of values for the PMF to be evaluated at
#' @param a Value of a in PA function
#' @param b Value of b in PA function
#' @param n0 Value of n0 in PA function
#' @param eps Value of eps in PA function
#' @param lambda Value of lambda usually found using find_lambda
#' 
#' @family distribution functions
#' 
#' @export
dppa = Vectorize(function(x, n0,a,b,eps,lambda,log=T){
  g = Vectorize(function(x,n0=1,a=1,b=1,eps=0){
    #nom
    # const = n0^(1-a)
    if(x<n0){
      return(x^a + eps)
    }else{
      return(max(0,n0)^a + b*(x-(max(0,n0))) + eps)
    }
  }, vectorize.args = 'x')
  if(x==1){
    return(
      -log(1+g(x,n0,a,b,eps)/lambda)
    )
  }
  return(
    -log(1+g(x,n0,a,b,eps)/lambda) - sum(log(1+lambda/g(1:(x-1),n0,a,b,eps)))
  )
}, vectorize.args = 'x')

#' Survival of PPA model
#' 
#' Returns log of survival
#' 
#' @param x Numeric or numeric vector of values for the PMF to be evaluated at
#' @param a Value of a in PA function
#' @param b Value of b in PA function
#' @param n0 Value of n0 in PA function
#' @param eps Value of eps in PA function
#' @param lambda Value of lambda usually found using find_lambda
#' 
#' @family distribution functions
#' 
#' @export
pppa = Vectorize(function(x, n0,a,b,eps,lambda,log=T){
  g = Vectorize(function(x,n0=1,a=1,b=1,eps=0){
    #nom
    # const = n0^(1-a)
    if(x<n0){
      return(x^a + eps)
    }else{
      return(max(0,n0)^a + b*(x-(max(0,n0))) + eps)
    }
  }, vectorize.args = 'x')
  return(-sum(log(1+lambda/g(1:x,n0,a,b,eps))))
}, vectorize.args = 'x')

#' MLE method for PPA model
#' 
#' This function will obtain the MLE of the parameters in the PPA model for a given data set.
#' 
#' @param dat A matrix with two columns of degrees and their counts
#' 
#' @family model fitting functions
#'
#' @export
ppa_mle <- function(dat){
  return(optim(c(10,1,1),fn=ppa_llh,dat=dat))
}
#'LL for ppa model
#'
#' @export
ppa_llh <- function(pars, dat){
  lambda = find_lambda_cpp(pars[1], pars[2], pars[3], 0)
  return(-llh_cpp(dat, lambda, pars[2], pars[3], pars[1]))
}

