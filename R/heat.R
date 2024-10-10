#'Make heatmap for value of xi
#'
#'@export
xi_heatmap = function(n0=15,focus_xi=0.5,a_lims = c(0,2), b_lims=c(0,2), res=200){
  a_vec = seq(a_lims[1], a_lims[2], l=res)
  b_vec = seq(b_lims[1], b_lims[2], l=res)
  par_grid = as.matrix(expand.grid(a_vec, b_vec, n0))
  
  
  cl = makeCluster(detectCores()-1)
  clusterEvalQ(cl, library('twbfn'))
  result = parApply(cl=cl, X=par_grid, MARGIN=1, FUN = find_lambda_par)
  stopCluster(cl)
  
  df = as.data.frame(par_grid)
  names(df) = c('a', 'b', 'n0')
  df$lambda = result
  df$xi=df$b/df$lambda
  
  plt = ggplot(data=df, aes(x=a, y=b, fill=xi)) + geom_raster() +
    scale_fill_gradient2(midpoint = focus_xi)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))
  return(plt)
}