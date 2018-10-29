#-------------function 2: calculate log intensity for a list of point process ---------------
est_inten_list_v2 = function(samples, Ns, T, mesh1, spde1, mlik1=FALSE, 
                         config1=FALSE, zgrid.eval=FALSE, tgrid, save_result=FALSE, take_sample=FALSE, verbose=FALSE){
# in _v2, some errors were corrected in E_point_process, and lam's estimation.
# Now the samples taken are in the correct scale. 
						 
  # This function estimates the log intensity of a list of point processes. 
  # Input:
  # samples: a list containing the multiple point processes.
  # Ns: a vector containing the number of point processes in each of of list component in samples. Default=1.
  #     If Ns!=2, that means, the point process is pooled from two separate point processes.
  # mesh1, spde1: these are items needed in inla.
  # zgrid.eval: whether evaluate the z(s) on tgrid or not. 
  # take_sample: if TRUE, instead of estimate the posterior mean of log(lam(s)), we will just take a random sample of it.
  # nsample: if take_sample==TRUE, the number of samples to take. 
  
  N = length(samples)
  if (zgrid.eval){
    lg_inten = matrix(NA, nrow=N, ncol=length(tgrid))
  }else{
	lg_inten = matrix(NA, nrow=N, ncol=mesh1$n)
  }
  if (length(Ns)<N){
    Ns = rep(Ns,N)
  }
  if (take_sample == TRUE){
    save_result = TRUE
    config1 = TRUE
  }
  #++++++++++++++
  # mlik1=FALSE 
  # config1=FALSE
  # zgrid.eval=FALSE
  # save_result = TRUE
  # win.graph()
  # plot(mesh1$loc, exp(outi$z.out[,1]), type='l')
  # z.tempi = inla.posterior.sample(n = 1, result = outi$res)[[1]]$latent
  # z.samplei = z.tempi[grep('idx.', rownames(z.tempi))]
  # lines(mesh1$loc, exp(z.samplei), type='l', col='red')
  #+++++++++++++++
  for (i in 1:N){
    # if (i%%5==0){
    # cat(paste('i=', i, '\n'))
    # }
    outi = lgcp_inla_v2(samples[[i]], Ns[i], T, mesh1, spde1, mlik1, config1, zgrid.eval, tgrid, save_result, verbose1=verbose)
    # the v2 version of lgcp_inla is a corrected version.
	
    if (take_sample == TRUE){
      z.tempi = inla.posterior.sample(n = 1, result = outi$res)[[1]]$latent		
      z.samplei = z.tempi[grep('idx.', rownames(z.tempi))]
    } else{
	  if (zgrid.eval){
	    z.samplei = outi$zgrid.out[,1]
	  }else{
		z.samplei = outi$z.out[,1]
	  }
    }
    lg_inten[i,] = z.samplei
  }
  
  return(lg_inten)
}
