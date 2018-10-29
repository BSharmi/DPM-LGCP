#-------------function 6: pool samples according to s------------------------------------
pool_list = function(Y_p, n_all, s){
  n1 = length(unique(s))
  expand_s = rep(s,times=n_all)
  new.sample = vector('list', n1)
  for (i in 1:n1){
    new.sample[[i]] = Y_p[expand_s==i,2]
  }
  return(new.sample)
}
