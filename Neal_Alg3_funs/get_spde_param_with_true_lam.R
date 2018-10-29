get_spde_param_with_true_lam <- function (lam_s, tgrid, d, plt){
  # This is the version of setting prior parameters for the SPDE model when the true lambda(s) is known. 
  # Reference: see test code at script_naive_simu_test_v5.R
  # Input:
  # lam.s: the true lam(s) evaluated on a grid of values. This should be a vector. 
  # tgrid: the time grid on which lam.s is evaluated. 
  # d: d=1 if s is 1-D, d=2 if s is 2-D. 
  # plt: 1/0 whether plot the auto correlations. 
  
  # Output:
  # sigma0: prior median for the marginal standard deviation of the GMRF prior. 
  # kappa0: the spatial scale parameter, smaller kappa0 means larger scale.  
  # tau0: the variance parameter, smaller tau0 means larger marginal variance. 
  
  # determine sig0, marginal standare deviation. 
  # sigma0 = 1 ## Prior median for the standard deviation, the one I used in version 1.
  acf.out1 = acf(log(lam_s), type="covariance", plot=FALSE)
  sigma0 = sqrt(acf.out1$acf[1])  # lag 0 of autocovariance gives marginal variance.
  
  # Determine kappa0, see page 4, "Lindgren and Rue, JSS/Bayesian spatial modelling with R-INLA" para. near(4).
  # kappa0 = 1e-3 ## Fixed kappa, the one I used in version 1.
  acf.out2 = acf(log(lam_s), plot=FALSE) # choose lag=17
  rho = tgrid[acf.out2$lag[which(abs(acf.out2$acf-0.13)==min(abs(acf.out2$acf-0.13)))]+1] # 0.13 is given in Lindgren and Rue JSS.
  
  # Determine tau0, see page 3, (2) of Lindgren and Rue, JSS.
  # tau0 = 1/(4*kappa0^3*sigma0^2)^0.5 ## Prior median for tau matching sigma0
  alpha = 2
  nu = alpha - d/2
  kappa0 = sqrt(8*nu)/rho ## Fixed kappa, if don't know lam.s, using (T[2]-T[1])/5
  tau0 = sqrt(gamma(nu)/(gamma(alpha)*(4*pi)^(d/2)*kappa0^(2*nu)*sigma0^2))
  
  if (plt==1){
    delta = mean(diff(tgrid))
    ngrid = length(tgrid)
    lag = c(0, (1:ngrid)*delta)
    par(mfrow=c(1,2))
    plot(lag[acf.out1$lag+1], acf.out1$acf, type='h', main='auto covariance')
    plot(lag[acf.out2$lag+1], acf.out2$acf, type='h', main='auto correlation')
  }
  
  return(list(sigma0 = sigma0, kappa0 = kappa0, tau0=tau0, range0=rho))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
Auto.cov <- function(dist, kappa, sigma, nu){
  # dist: a vector of distances, starting from 0, ends till T.
  
  auto.cor =  1/(2^(nu-1)*gamma(nu))*(kappa*dist)^(nu)*besselK(kappa*dist,nu)
  auto.cor[1] = 1
  auto.cov = sigma^2*auto.cor
  
  
  return(list(cor=auto.cor, cov=auto.cov))
}