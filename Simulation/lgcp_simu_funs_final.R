gen_intensity = function(tgrid, l, scaling)
  # generate a single intensity function using log-Gaussian process with zero mean and squared exponential covariance
  # input:
  #       tgrid: vector, time grid of integer values, denote (CSM site) positions in real data
  #       l: scalar, characteristic length scale
  #       scaling: scalar, scaling parameter for intensity function
  # output:
  ##       a vector of length(tgrid), an intensity function scaled by scaling parameter
{
  library(MASS)
  p = length(tgrid)
  mu = rep(0, p)
  rmat = matrix(NA, p, p)
  for (i in 1 : p)
  {
    rmat[i, ] = tgrid - tgrid[i]
  }
  Sigma = exp(-rmat ^ 2 / 2 / l ^ 2)
  fGP = mvrnorm(1, mu, Sigma)
  
  return(fint = exp(fGP) / scaling)
}

inhomo_poisson_sampler = function(tgrid, fint)
  # generate a sample of non-homogeneous Poisson process with intensity function
  # input:
  #       tgrid: vector, time grid of integer values, denote (CSM site) positions in real data
  #       fint: vector, an intensity function of length(tgrid)
  # output:
  #       a list with
  #       sample: vector
{
  if (length(tgrid) != length(fint))
  {
    stop("Wrong!")
  }
  # (1) Compute the volumn of the region. 
  omega = max(fint) + 0.01
  TT = range(tgrid)
  vol = TT[2] - TT[1]

  # (2) Generate a poisson number
  num.allevents = rpois(1, omega * vol)
  
  # (3) Conditional on the number, generate uniform event locations.
  sample.allevents = runif(num.allevents, min = min(TT), max = max(TT))

  # (4) Generate GP at the locations.
  Y = approx(tgrid, fint, xout = sample.allevents)$y
  
  # (5) Thinning
  U = runif(num.allevents, min = 0, max = omega)
  ic = (U < Y)
  sample = sample.allevents[ic]

  return(list(sample = sample, sample.allevents = sample.allevents, thin.rate = sum(ic) / num.allevents))
}
