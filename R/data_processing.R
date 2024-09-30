#' @export
deg_count = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  names(df) = c('degree', 'count')
  return(df)
}
#' @export
deg_surv = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  df[,2] = 1 - cumsum(df[,2])/sum(df[,2])
  names(df) = c('degree', 'surv')
  return(df)
}
#' @export
deg_dist = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  df[,2] = df[,2]/sum(df[,2])
  names(df) = c('degree', 'prob')
  return(df)
}

