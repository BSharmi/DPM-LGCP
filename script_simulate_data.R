# generate independent samples from 1D non-homogeneous Poisson process with log-Gaussian intensity
# ------------------------------------------------------------------------------------------------
#rm(list = ls())
## change work path as needed. It should include full path upto LGCP_package
#workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/"
source(paste(workpath, "Simulation/lgcp_simu_funs_final.R", sep = "/"))

TT = c(0, 10)
tgrid = seq(TT[1], TT[2], length = 50)

set.seed(21)
f1 = gen_intensity(tgrid, 2, 1.8) # 2 (0.3635->0.0105)
f2 = gen_intensity(tgrid, 2, 1.8)
f3 = gen_intensity(tgrid, 2, 1.8)
f.vec = cbind(f1, f2, f3)

## 20 TFs randomly assigned to three clusters
nsample = 20
clust = c(rep(1, 6), rep(2, 7), rep(3, 7)) ## 

samples = vector('list', nsample)
sample.mat = NULL
sample.len = rep(NA, nsample)
set.seed(5000)
for (i in 1 : nsample){
  samples[[i]] = inhomo_poisson_sampler(tgrid, f.vec[, clust[i]])$sample
  sample.mat = rbind(sample.mat, cbind(i, samples[[i]]))
  sample.len[i] = length(samples[[i]])
}

## save simulation data
save(tgrid, f.vec, samples, clust, sample.mat, sample.len, file = paste(workpath, "Simulation/data/lgcp_simu_data_v0.RData", sep = "/"))

