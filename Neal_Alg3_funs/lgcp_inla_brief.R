lgcp_inla_brief = function(Y, Ns=1, T, mesh1, spde1, mlik1, verbose1=FALSE){
  # This function is developed based on lgcp_inla.R, it is a brief version---in this version
  # I only take the marginal likelihood, no estimation of lambda (or Z value) is done. 
  
  # Input:
  # Y: a vector containing the location of the events. 
  # Ns: the total number of samples that consistute Y. If Ns is greater than 1, we need to adjust for
  #     E_point_process, and the marginal likelihood. See my derivations. 
  # mesh1: the mesh on which we approximate the GMRF
  # spde1: the spde object. 
  
  nV=mesh1$n
  nData = length(Y)
  IntegrationMatrix = sparseMatrix(i=1:nV,j=1:nV,x=rep(1,nV)) #simple integration scheme #IntegrationMatrix1 = inla.mesh.project(mesh, mesh$loc)$A
  LocationMatrix = inla.mesh.project(mesh1, Y)$A   # same for : A = inla.spde.make.A(mesh, loc=y) # N by p. 
  ObservationMatrix = rBind(IntegrationMatrix,LocationMatrix) #Glue matrices together  
  IntegrationWeights = diag(Ns*(inla.mesh.fem(mesh1)$c0)) # Finn Lindgren said should use this line. 
  # when Ns>1, we need to rescale alpha_tilde to N*alpha_tilde since we need to multiply Ns likelihood functions.  
  E_point_process =c(IntegrationWeights,rep(0,nData)) #Glue everything together
  fake_data = c(rep(0,nV),rep(1,nData)) #Put a zero where there aren't observations and a 1 where there is a point
  
  formula = y ~ -1 + f (idx, model=spde1) #Basic latent model - feel free to add covariates etc
  data = list(y=fake_data, idx = c(1:nV)) #put the data in
  
  #The INLA call. Likelihood is Poisson with Observation Matrix and appropriate value fo E.
  result = inla(formula, data=data, family="poisson", control.predictor=list(A=ObservationMatrix),
                E=E_point_process, control.compute=list(mlik=mlik1, config=FALSE), verbose=verbose1)  
  
  # mean.coef=result$summary.random$idx[,"mean"]-log(Ns)
  # quant.025=result$summary.random$idx[,"0.025quant"]-log(Ns)
  # quant.975=result$summary.random$idx[,"0.975quant"]-log(Ns)
  # z.out = cbind(mean.coef, quant.025, quant.975)
  
  
  # if (zgrid.eval==TRUE){
    # proj = inla.mesh.projector(mesh1, dims = length(tgrid))
    # gridvalue.mean = inla.mesh.project(proj, mean.coef)-log(Ns)
    # gridvalue.025 =  inla.mesh.project(proj, quant.025)-log(Ns)
    # gridvalue.975 =  inla.mesh.project(proj, quant.975)-log(Ns)
    # zgrid.out = cbind(gridvalue.mean, gridvalue.025, gridvalue.975)	
  # }else{
    # zgrid.out = NULL
  # }
  
  
  # mlike = result$mlik[[1]] + Ns*(T[2]-T[1])-nData # here added the |C| to the log likelihood. 
  mlike = result$mlik[[1]] + Ns*(T[2]-T[1])    
  # when Ns>1, we need to rescale |D| to Ns*|D| since we need to multiply Ns likelihood functions.  
  
  
  # if (save_result){
    # res = result    
  # }else{
    # res = NULL
  # }
  # return(list(z.out=z.out, zgrid.out=zgrid.out, mlike = mlike, res=res))
  return(mlike)
}
