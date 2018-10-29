#-------------get_mesh_spde------------------------------------
get_mesh_spde = function(nknots = 50, TT, sigma0, kappa0, tau0){
  # This function get mesh on 1-d domain for use in INLA. 
  # Input:
  # nknots: how many knots to use.
  # T: the boundary of the domain. 
  
  knots = seq(TT[1], TT[2], length = nknots) # prespecified knots that does not contain data locations. 
  mesh0 = inla.mesh.1d(knots, interval = TT, degree = 1, boundary = "free")
  
  spde0 = inla.spde2.matern(mesh0, constr = FALSE,
                            B.tau = cbind(log(tau0)),
                            B.kappa = cbind(log(kappa0))) 
  
  # Note: results will be sensitive to the three parameters here. Need to choose them
  # very carefully. 
  #sigma0 = 1.0 ## Prior median for the standard deviation
  #kappa0 = 1e-4 ## Fixed kappa
  #tau0 = 1/(4*kappa0^3*sigma0^2)^0.5 ## Prior median for tau matching sigma0

  # Note: I will always have inla crashing problem when allow 1 in B.tau. 
  # spde = inla.spde2.matern(mesh, constr=FALSE,
  #                          B.tau = cbind(log(tau0), 1),
  #                          B.kappa = cbind(log(kappa0), 0),
  #                          theta.prior.prec = 1e-6) # theta.prior.prec = 1e-5
  
  # spde = inla.spde2.matern(mesh, constr=FALSE,
  #                          B.tau = cbind(log(tau0), 1),
  #                          B.kappa = cbind(log(kappa0), 0),
  #                          theta.prior.prec = 1e-6) # theta.prior.prec = 1e-5
  
  return(list(mesh = mesh0, spde = spde0))
}
