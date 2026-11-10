
#' @useDynLib twbfn, .registration = TRUE
#' @importFrom Rcpp evalCpp


counts_to_surv = function(df){
  return(
    twbfn::deg_surv(rep(df[,1], df[,2]))
  )
}

dggp_tail_vec = Vectorize(dggp_tail,vectorize.args = 'x')
dggp_rho_vec = Vectorize(dggp_rho, vectorize.args = 'x')

stable_tail_vec = Vectorize(stable_tail,vectorize.args = 'x')
stable_rho_vec = Vectorize(stable_rho, vectorize.args = 'x')

ggp_tail_vec = Vectorize(ggp_tail,vectorize.args = 'x')
ggp_tail_invgrad_vec = Vectorize(ggp_tail_invgrad,vectorize.args = 'x')
ggp_rho_vec = Vectorize(ggp_rho, vectorize.args = 'x')
gg_tail_vec = Vectorize(gg_tail,vectorize.args = 'x')
gg_rho_vec = Vectorize(gg_rho, vectorize.args = 'w')

#' @export
ggp_sample = function(N, sigma, tau, beta){
  
  N_mean = N
  wmin = 1e-5
  wscale = 1
  
  x = 10^seq(-6,6,l=10000)
  y = ggp_tail_vec(x,sigma,tau,beta)
  ggp_tail_inv = approxfun(y,x)
  
  nu = N_mean/ggp_tail_vec(wmin, sigma, tau, beta)
  N = rpois(1,N_mean)
  u = runif(N,0,ggp_tail_vec(wmin, sigma, tau, beta))
  w = ggp_tail_inv(u)
  w = w[!is.na(w)]
  return(w)
  
  
}

stable_tail_inv = function(x, sigma){
  return(
    gamma(1-sigma)^(-1/sigma) * x^(-1/sigma)
  )
}

stable_sample = function(N, sigma){
  
  N_mean = N
  wmin = 1e-5
  wscale = 1
  
  x = 10^seq(-6,6,l=10000)
  y = stable_tail_vec(x,sigma)
  
  nu = N_mean/stable_tail_vec(wmin, sigma)
  N = rpois(1,N_mean)
  u = runif(N,0,stable_tail_vec(wmin, sigma))
  w = stable_tail_inv(u,sigma)
  w = w[!is.na(w)]
  return(w)
}

#' @export
dggp_sample = function(N, sigma, tau, beta,p){
  
  N_mean = N
  wmin = 1e-5
  wscale = 1
  
  x = 10^seq(-6,6,l=10000)
  y = dggp_tail_vec(x,sigma,tau,beta,p)
  dggp_tail_inv = approxfun(y,x)
  dggp_tail_vec_min = dggp_tail_vec(wmin, sigma, tau, beta,p)
  
  nu = N_mean/dggp_tail_vec_min
  N = rpois(1,nu*dggp_tail_vec_min)
  u = runif(N,0,dggp_tail_vec_min)
  w = dggp_tail_inv(u)
  w = w[!is.na(w)]
  return(w)
  
  
}

#' @export
gg_sample = function(N, sigma, beta){
  
  
  N_mean = N
  wmin = 1e-5
  wscale = 1
  
  x = 10^seq(-6,6,l=10000)
  y = gg_tail_vec(x,sigma,beta)
  gg_tail_inv = approxfun(y,x)
  
  nu = N_mean/gg_tail_vec(wmin, sigma, beta)
  N = rpois(1,nu*gg_tail_vec(wmin, sigma, beta))
  u = runif(N,0,gg_tail_vec(wmin, sigma, beta))
  w = gg_tail_inv(u)
  w = w[!is.na(w)]
  return(w)
}
gamma_upper = Vectorize(function(a, z){
  if(a<0){
    return(
      (gamma_upper(a+1, z)-z^a *exp(-z))/a
    )
  }else{
    return(tgamma_upper(a, z))
  }
}, vectorize.args = 'z')
integrand_gg = function(v,w, sigma, beta){
  return(
    -expm1(-2*w*v)*gg_rho(v, sigma, beta)
  )
}
gg_mu = Vectorize(function(w, sigma, beta){
  return(
    pracma::quadinf(integrand_gg, 0, Inf, w=w, sigma=sigma,  beta=beta)$Q
  )
}, vectorize.args = 'w')
integrand_ggp = function(v,w, sigma, tau,beta){
  return(
    -expm1(-2*w*v)*ggp_rho(v, sigma,tau, beta)
  )
}
ggp_mu = Vectorize(function(w, sigma,tau, beta){
  return(
    pracma::quadinf(integrand_ggp, 0, Inf, w=w, sigma=sigma,tau=tau,  beta=beta)$Q
  )
}, vectorize.args = 'w')

#' @export
simggp = function(N, sigma, tau, beta){
  w_ggp = ggp_sample(N,sigma,tau,beta)
  edges_ggp = sample_edges(w_ggp,scaling=1)
  edges_ggp[,1] = as.character(edges_ggp[,1])
  edges_ggp[,2] = as.character(edges_ggp[,2])
  
  G = igraph::graph_from_edgelist(as.matrix(edges_ggp), directed = F)
  return(G)
}
simstable = function(N, sigma){
  w_ggp = stable_sample(N,sigma)
  edges_ggp = sample_edges(w_ggp,scaling=1)
  edges_ggp[,1] = as.character(edges_ggp[,1])
  edges_ggp[,2] = as.character(edges_ggp[,2])
  
  G = igraph::graph_from_edgelist(as.matrix(edges_ggp), directed = F)
  return(G)
}

#' @export
simgg = function(N, sigma, beta){
  w_gg = gg_sample(N,sigma,beta)
  edges_gg = sample_edges(w_gg,scaling=1)
  edges_gg[,1] = as.character(edges_gg[,1])
  edges_gg[,2] = as.character(edges_gg[,2])
  
  G = igraph::graph_from_edgelist(as.matrix(edges_gg), directed = F)
  return(G)
}
exp_integrand_ggp = function(w, j, nu, sigma, tau, beta){
  return(
    pgamma(nu*ggp_mu(w, sigma, tau, beta), shape= j+1, scale=1)*ggp_rho_vec(w, sigma, tau, beta)
  )
}
est_ggp = Vectorize(function(j,p, sigma, tau, beta,nu ){
  return(
    (log(pracma::quadinf(exp_integrand_ggp, 0, Inf, j=j, nu=p*nu, sigma=sigma, tau=tau, beta =beta)$Q)-
       log(pracma::quadinf(exp_integrand_ggp, 0, Inf, j=j, nu=nu, sigma=sigma, tau=tau, beta =beta)$Q))/log(p)
  )
}, vectorize.args = 'j')
integrand_ggp = function(v,w, sigma, tau,beta){
  return(
    -expm1(-2*w*v)*ggp_rho(v, sigma,tau, beta)
  )
}
exp_integrand_ggp = function(w, j, nu, sigma, tau, beta){
  return(
    pgamma(nu*ggp_mu(w, sigma, tau, beta), shape= j+1, scale=1)*ggp_rho_vec(w, sigma, tau, beta)
  )
}
est_ggp = Vectorize(function(j,p, sigma, tau, beta,nu ){
  return(
    (log(pracma::quadinf(exp_integrand_ggp, 0, Inf, j=j, nu=p*nu, sigma=sigma, tau=tau, beta =beta)$Q)-
       log(pracma::quadinf(exp_integrand_ggp, 0, Inf, j=j, nu=nu, sigma=sigma, tau=tau, beta =beta)$Q))/log(p)
  )
}, vectorize.args = 'j')
exp_integrand_gg = function(w, j, nu, sigma, beta){
  return(
    pgamma(nu*gg_mu(w, sigma, beta), shape= j+1, scale=1)*gg_rho_vec(w, sigma, beta)
  )
}
est_gg = Vectorize(function(j,p, sigma, beta,nu ){
  return(
    (log(pracma::quadinf(exp_integrand_gg, 0, Inf, j=j, nu=p*nu, sigma=sigma, beta =beta)$Q)-
       log(pracma::quadinf(exp_integrand_gg, 0, Inf, j=j, nu=nu, sigma=sigma, beta =beta)$Q))/log(p)
  )
}, vectorize.args = 'j')


#' @export
estimator = Vectorize(function(degcounts, j, p=0.5){
  temp = degcounts[degcounts[,1]>j,]
  probs = p*pbeta(p,j+1,temp[,1]-j)
  new = rbinom(nrow(temp), temp[,2],probs)
  est = (log(sum(new)) - log(sum(temp[,2])))/log(p) - 1
  return(est)
}, vectorize.args = c('j','p'))

#' @export
estimator_mean = Vectorize(function(degcounts, j, p=0.5){
  temp = degcounts[degcounts[,1]>j,]
  probs = p*pbeta(p,j+1,temp[,1]-j)
  est = (log(sum(temp[,2]*probs)) - log(sum(temp[,2])))/log(p) - 1
  return(est)
}, vectorize.args = c('j','p'))


#' @export
estimator_var = Vectorize(function(degcounts, j, p){
  temp = degcounts[degcounts[,1]>j,]
  probs = p*pbeta(p,j+1,temp[,1]-j)
  v = sum(temp[,2]*probs*(1-probs))
  mn = sum(temp[,2]*probs)
  return(
    v/(mn^2 * log(p)^2)
  )
}, vectorize.args = 'j')
rep_estimator = function(degcounts, p,reps){
  df = as.data.frame(rep(degcounts[,1], reps))
  names(df) = 'degree'
  df$est  = estimator(degcounts, df$degree,p)
  df = df[!is.nan(df$est),]
  df = df[!is.infinite(df$est),]
  return(df)
}
estimator_compress = Vectorize(function(degcounts, j, p=0.5){
  temp = degcounts[degcounts[,1]>j,]
  p = pmax(p, j/max(temp[,1]))
  probs = p*pbeta(p,j+1,temp[,1]-j)
  new = rbinom(nrow(temp), temp[,2],probs)
  est = (log(sum(new)) - log(sum(temp[,2])))/log(p) - 1
  return(est)
}, vectorize.args = c('j','p'))
rep_estimator_compress = function(degcounts, p,reps){
  df = as.data.frame(rep(degcounts[,1], reps))
  names(df) = 'degree'
  df$est  = estimator_compress(degcounts, df$degree,p)
  df = df[!is.nan(df$est),]
  df = df[!is.infinite(df$est),]
  return(df)
}
rep_estimator_long = function(degcounts, p,reps){
  df = as.data.frame(rep(degcounts[,1], reps))
  names(df) = 'degree'
  est  = estimator(degcounts, rep(1:(max(df$degree)-1), reps),p)
  df = data.frame(degree = rep(1:(max(df$degree)-1), reps), est = est)
  df = df[!is.nan(df$est),]
  df = df[!is.infinite(df$est),]
  return(df)
}
rep_estimator_mean = function(degcounts, p,reps){
  df = as.data.frame(rep(degcounts[,1], reps))
  names(df) = 'degree'
  df$est  = estimator(degcounts, df$degree,p)
  df = df[!is.nan(df$est),]
  df = df[!is.infinite(df$est),]
  mdf = aggregate(data=df, est~degree, FUN=mean)
  return(mdf)
}
#' @export
hill = function(x,ks=NULL){
  X = sort(x, decreasing=FALSE)
  n = length(X)
  if(is.null(ks)){
    ks = 1:max(X)
  }
  hills = (sapply(ks, function(k){
    js = 1:k
    sum(log(X[n-js+1]) - log(X[n-k]))/(k) 
  }))
  
  return(
    cbind(X[n-ks+1], 1/hills)
  )
}

#' @export
naulet_est = function(degcounts, prob){
  df = degcounts
  N = sum(df[,2])
  ests = sapply(prob, FUN = function(p){
    Np = p*sum(df[,2] * (1-(1-p)^df[,1]))
    (log(Np) - log(N))/log(p) - 1
  })
  return(data.frame(p = prob, est = ests))
}

estimator_le = Vectorize(function(degcounts, j, p=0.5){
  temp = degcounts
  probs = p*(pbeta(1-p,pmax(temp[,1]-j,0),j+1) - (1-p)^temp[,1])
  probs[temp[,1]<=j] = 1
  new = rbinom(nrow(temp), temp[,2],probs)
  est = -(log(sum(new)) - log(sum(temp[temp[,1]<=j,2])))/log(p)
  return(est)
}, vectorize.args = c('j','p'))

rep_estimator_le_long = function(degcounts, p,reps){
  df = as.data.frame(rep(degcounts[,1], reps))
  names(df) = 'degree'
  est  = estimator_le(degcounts, rep(1:(max(df$degree)-1), reps),p)
  df = data.frame(degree = rep(1:(max(df$degree)-1), reps), est = est)
  df = df[!is.nan(df$est),]
  df = df[!is.infinite(df$est),]
  return(df)
}

#' @export
true_estimator = function(G, p){
  selected = rbinom(vcount(G), 1, p)
  degs = degree(G)
  sg = subgraph(G, V(G)[selected==1])
  pdegs= degree(sg)
  js= unique(degs)
  ests = sapply(js, FUN = function(j){
    return(
      (log(sum(pdegs>j)) - log(sum(degs>j)))/log(p) - 1
    )
  })
  return(data.frame(degree=js, est=ests))
}

#' @export
rep_true_estimator = function(G, p, reps){
  df = data.frame(degree=NA, est=NA)
  for(i in 1:reps){
    df = rbind(df, true_estimator(G, p))
  } 
  return(df)
}

#' @export
estimator_plot = function(degcounts,ps=NULL){
  if(is.null(ps)){
    ps = seq(0.5,0.9,l=9)
  }
  df= data.frame(x=numeric(), mean=numeric(), sd = numeric(), p=numeric())
  for(p in ps){
    message(p)
    x = degcounts[degcounts[,1]<=(p*max(degcounts[,1])),1]
    mn = estimator_mean(degcounts, x, p)
    std = sqrt(estimator_var(degcounts, x, p))
    temp = data.frame(x=x, mean=mn, sd = std)
    temp$p=p
    df  =rbind(df, temp)
  }
  df = df[!is.na(df$mean),]
  df = df[is.finite(df$mean),]
  ggplot(df) + 
    geom_ribbon(aes(x=x, ymin=mean-sd*1.96, ymax = mean+sd*1.96, group=p, fill=as.character(p)), alpha=0.2)+
    geom_ribbon(aes(x=x, ymin=mean-sd*1.96, ymax = mean+sd*1.96, group=p, fill=as.character(p)), alpha=0.2)+
    geom_line(aes(x=x, y=mean,group=p, color=as.character(p))) + 
    coord_cartesian(ylim=c(1e-2,min(10, max(df$mean*1.5))), xlim=c(1,max(degcounts[,1]))) +
    facet_wrap(vars(p))+
    theme_light()+
    theme(legend.position = 'none', aspect.ratio = 1)
}
