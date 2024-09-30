#' @export
sample_ppa = function(N,a,b,n0,smplr=NULL, ..., init_size=25, quiet=T){
  if(is.null(smplr)){
    return(sample_ppa_tree(N,a,b,n0))
  }else{
    return(sample_ppa_rnd(N,a,b,n0,smplr,...,init_size=init_size, quiet=quiet))
  }
}
#' @export
sample_ppa_tree = function(N,a,b,n0){
  ptrs = pars_to_ptrs2(a,b,n0)
  
  ctrl = wdnet::rpa_control_preference(ftype='customized',
                                pref=ptrs$s,
                                spref = ptrs$outs,
                                tpref = ptrs$ins) +
    wdnet::rpa_control_scenario(beta.loop = FALSE)
  init = list(edgelist = matrix(c(1,1), nrow = 1), edgeweight = 1, directed = TRUE)
  ret  =wdnet::rpanet(N, control=ctrl, initial.network = init)
  G = wdnet::wdnet_to_igraph(ret)
  return(G)
}
#' @export
sample_ppa_rnd = function(N,a,b,n0,smplr,..., init_size = 25, quiet=T){
  m = init_size
  my_sampler = function(n){
    out = smplr(n,...)
    return(out)
  }
  g = Vectorize(function(x,n0=1,a=1,b=1,eps=0){
    #nom
    # const = n0^(1-a)
    if(x<n0){
      return(x^a + eps)
    }else{
      return(max(0,n0)^a + b*(x-(max(0,n0))) + eps)
    }
  }, vectorize.args = 'x')
  
  
  ptrs = pars_to_ptrs2(a,b,n0)
  
  control <- wdnet::rpa_control_scenario(alpha = 1, beta = 0,gamma=0, source.first = F)+ rpa_control_newedge(
    sampler = my_sampler,
    snode.replace = F,
    tnode.replace = F
  ) + wdnet::rpa_control_preference(ftype='customized',
                             pref=ptrs$s,
                             spref = ptrs$outs,
                             tpref = ptrs$ins)
  
  old_nrows = m
  ret1 <- wdnet::rpanet(nstep = 1, control = control, initial.network = list(edgelist = cbind(c(1:m), c(1:m)), directed=T))
  
  for(i in 1:N){
    if(!quiet) message(i,'/',N)
    new_nrows = nrow(ret1$edgelist)
    el = ret1$edgelist
    change = new_nrows - old_nrows
    old_nrows=new_nrows
    el[(nrow(el)-change+1) : nrow(el),1]  = min(el[(nrow(el)-change+1) : nrow(el),1])
    ret1<- wdnet::rpanet(nstep = 1, control=control, initial.network = list(edgelist = el, directed=T))
  }
  
  G = wdnet::wdnet_to_igraph(ret1)
  
  return(G)
  
}