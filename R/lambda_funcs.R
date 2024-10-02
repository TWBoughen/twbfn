lambda_obj = Vectorize(function(l, n0,a,b,eps, mx=1e3){
  g = Vectorize(function(x,n0=1,a=1,b=1,eps=0){
    #nom
    # const = n0^(1-a)
    if(x<n0){
      return(x^a + eps)
    }else{
      return(max(0,n0)^a + b*(x-(max(0,n0))) + eps)
    }
  }, vectorize.args = 'x')
  n = 1:mx
  lgs  = log(1 + l/g(n,n0,a,b,eps))
  return(sum(exp(-cumsum(lgs))))
}, vectorize.args = 'l')
lambda_obj_adj = function(l,n0,a,b,eps,mx=1e3){
  return(1-lambda_obj(l,n0,a,b,eps,mx=mx))
}


#' Finding the value of \code{lambda}
#' 
#' Function to be used to find the value of the constant \code{lambda} in the Rudas formula for a set of parameters.
#' 
#' @export
find_lambda = function(n0,a,b,eps){
  return(uniroot(lambda_obj_adj, interval=c(0,10),n0=n0,a=a,b=b,eps=eps,extendInt = 'upX'))
}