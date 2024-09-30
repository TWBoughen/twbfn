pars_to_ptrs2 = function(a,b,n0){
  bma = as.character(b-a)
  pt_str = "(((x+1)<n0) ? (pow(x+1,A)):pow(n0,A)+B*pow(x+1-n0,1))"
  var_names = c('outs', 'ins', 's')
  ptr_list = list(outs='', ins='i', s='')
  for(i in 1:length(var_names)){
    t_str = gsub('x', var_names[i],pt_str)
    t_str= gsub('BMA', bma,t_str)
    t_str = gsub('n0',as.character(n0),t_str)
    t_str = gsub('A', as.character(a), t_str)
    t_str = gsub('B', as.character(b), t_str)
    ptr_list[[i]] =t_str
  }
  return(ptr_list)
}