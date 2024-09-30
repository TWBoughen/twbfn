#' @export
dppa = Vectorize(function(x, n0,a,b,eps,lambda,log=T){
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
  return(-sum(log(1+lambda/g(1:x,n0,a,b,eps))))
}, vectorize.args = 'x')