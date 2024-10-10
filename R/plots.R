pppa_par = function(pars,vals){
  #n0, a, b, eps, lambda
  return(exp(pppa(vals, pars[1], pars[2], pars[3], 0, pars[4])))
}
dppa_par = function(pars,vals){
  #n0, a, b, eps, lambda
  return(exp(dppa(vals, pars[1], pars[2], pars[3], 0, pars[4])))
}

#'
#'@export
mcmc_surv_plot = function(res, burn_in, alpha=0.95){
  df = data.frame(a = res[[1]], b=res[[2]], n0=res[[3]], lambda = res[[4]])
  df = df[-(1:burn_in+1e3),]
  x = res$data[,1]
  df = df[!duplicated(df),]
  df = df[,c(3,1,2,4)]
  df=as.matrix(df)
  cl = makeCluster(detectCores()-1)
  clusterEvalQ(cl, library('twbfn'))
  result = parApply(cl = cl, X = df, MARGIN=1, FUN=pppa_par, vals=x)
  stopCluster(cl)
  S_quants = apply(result, 1, quantile, probs = c((1-alpha)/2, 1 - (1-alpha)/2))
  S_mean = apply(result, 1, mean)
  full_df = cbind(x, S_mean, t(S_quants))
  colnames(full_df) = c('x', 'S_mean', paste0('S_', (1-alpha)/2), paste0('S_',1 - (1-alpha)/2))
  plt = ggplot(data=full_df, aes(x=x, y=S_mean)) + geom_line(colour='red') +
    geom_ribbon(aes(ymin = full_df[,paste0('S_', (1-alpha)/2)], ymax = full_df[,paste0('S_',1 - (1-alpha)/2)]),
                alpha=0.4)+
    geom_point(data=NULL, aes(x=res$data[,1], y=1-cumsum(res$data[,2])/sum(res$data[,2])))+
    scale_x_log10(expand=c(0,0)) + scale_y_log10(expand=c(0,0), limits= c(NA,1)) + 
    ylab('1-F(n)') + xlab('degree (n)')
  return(plt)
}
#'
#'@export
mcmc_pmf_plot = function(res, burn_in, alpha=0.95){
  df = data.frame(a = res[[1]], b=res[[2]], n0=res[[3]], lambda = res[[4]])
  df = df[-(1:burn_in+1e3),]
  x = res$data[,1]
  df = df[!duplicated(df),]
  df = df[,c(3,1,2,4)]
  df=as.matrix(df)
  cl = makeCluster(detectCores()-1)
  clusterEvalQ(cl, library('twbfn'))
  result = parApply(cl = cl, X = df, MARGIN=1, FUN=dppa_par, vals=x)
  stopCluster(cl)
  f_quants = apply(result, 1, quantile, probs = c((1-alpha)/2, 1 - (1-alpha)/2))
  f_mean = apply(result, 1, mean)
  full_df = cbind(x, f_mean, t(f_quants))
  colnames(full_df) = c('x', 'f_mean', paste0('f_', (1-alpha)/2), paste0('f_',1 - (1-alpha)/2))
  plt = ggplot(data=full_df, aes(x=x, y=f_mean)) + geom_line(colour='red') +
    geom_ribbon(aes(ymin = full_df[,paste0('f_', (1-alpha)/2)], ymax = full_df[,paste0('f_',1 - (1-alpha)/2)]),
                alpha=0.4)+
    geom_point(data=NULL, aes(x=res$data[,1], y=res$data[,2]/sum(res$data[,2])))+
    scale_x_log10(expand=c(0,0)) + scale_y_log10(expand=c(0,0), limits= c(NA,1)) + 
    ylab('f(n)') + xlab('degree (n)')
  return(plt)
}
