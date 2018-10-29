log_marginal_posterior_lgcp = function(x, Ns=1, T, mesh1, spde1, xs=NULL){
# x: a point process
# Ns=1: 1 if x is a single point process, >1 if it is pooled from multiple point processes.
# xs: can be a list of point processes.  
 
if (length(xs) > 0){
    #l = data_likelihood_betabernoulli([x;xs], G0alphabeta, ones(1,N+1)) - ...
    #    data_likelihood_betabernoulli(   xs , G0alphabeta, ones(1,N  ))
  } else{
    l = data_likelihood_lgcp(x, Ns, T, mesh1, spde1, 1) # done.	
 }
 
 return(l)
}