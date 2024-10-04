\name{twbfn}
\alias{twbfn}
\alias{twbfn}
\docType{package}
\title{
  Collection of functions for Thomas Boughen's PhD
}
\description{
  This R package has been made in order to collect the functions that I have been using for sampling from certain network models and investigating their properties.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions. 
  
  The piecewise-preferential attachement (PPA) model which is referred to throughout this package
  is a general preferential attachment model as described by Rudas et al.  where the preferential attachment function takes the form:
  
  \deqn{g(n) = \begin{cases}n^a,&n<n_0\\ n_0^a + b(n-n_0),&n\ge n_0 \end{cases}}
}
\author{
Thomas W Boughen,<t.w.boughen1@newcastle.ac.uk>.

Maintainer: Thomas W Boughen <t.w.boughen1@newcastle.ac.uk>
}
\references{
 Rudas, A., Tóth, B. and Valkó, B. (2007), https://doi.org/10.1002/rsa.20137
}
\keyword{ package }
\seealso{
  \bold{Theory}
  
  {\code{\link{dppa}} }{PMF of PPA model}
  
  {\code{\link{pppa}} }{Survival of PPA model}
  
  \bold{Simulation}
  
  \code{\link{sample_ppa} }{Sample from PPA model}
  
  \bold{MCMC}
  
  \code{\link{ppa_mcmc} }{Parameter estimation for PPA model}
  
  \bold{Misc}
  
  \code{\link{deg_count}}{ Get degree counts} 
  
  \code{\link{deg_dist}}{ Get empirical degree distribution}
  
  \code{\link{deg_surv}}{ Get empirical degree survival}
  
  \code{\link{find_lambda} }
}