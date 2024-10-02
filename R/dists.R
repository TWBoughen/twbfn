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