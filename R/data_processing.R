#' Get degree counts
#' Returns the degree counts from a degree sequence
#'
#' @param degs integer vector of degrees
#' 
#' @export
deg_count = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  names(df) = c('degree', 'count')
  return(df)
}


#'Get degree survival
#'
#'Returns empirical survival from degree sequence as a data frame
#'
#' @param degs integer vector of degrees
#'
#' @export
deg_surv = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  df[,2] = 1 - cumsum(df[,2])/sum(df[,2])
  names(df) = c('degree', 'surv')
  return(df)
}


#'Get degree PMF
#'
#'Returns empirical PMF from degree sequence
#'
#'@param degs integer vector of degrees
#'
#' @export
deg_dist = function(degs){
  df = as.data.frame(table(degs))
  df[,1] = as.numeric(names(table(degs)))
  df[,2] = df[,2]/sum(df[,2])
  names(df) = c('degree', 'prob')
  return(df)
}

