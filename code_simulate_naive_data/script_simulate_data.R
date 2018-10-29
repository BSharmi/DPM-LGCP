# generate independent samples from 1D non-homogeneous Poisson process with log-Gaissian intensity
# ------------------------------------------------------------------------------------------------
rm(list = ls())
workpath = "/Users/bsharmi6/Documents/R/BI_ChIPSeq/SharmiWork"
source(paste(workpath, "code_simulate_naive_data/lgcp_simu_funs_final.R", sep = "/"))

TT = c(0, 10)
tgrid = seq(TT[1], TT[2], length = 50)

set.seed(21)
f1 = gen_intensity(tgrid, 2, 1.8) # 2 (0.3635->0.0105)
f2 = gen_intensity(tgrid, 2, 1.8)
f3 = gen_intensity(tgrid, 2, 1.8)
f.vec = cbind(f1, f2, f3)
# mu = c(1, 5, 9)
# sd = rep(0.6, 3)
# N_expect = 5
# dens1 = N_expect*dnorm(tgrid, mean = mu[1], sd=sd[1])
# dens2 = N_expect*dnorm(tgrid, mean = mu[2], sd=sd[2])
# dens3 = N_expect*dnorm(tgrid, mean = mu[3], sd=sd[3])

nsample = 20
clust = c(rep(1, 6), rep(2, 7), rep(3, 7))

samples = vector('list', nsample)
sample.mat = NULL
sample.len = rep(NA, nsample)
set.seed(5000)
for (i in 1 : nsample){
  samples[[i]] = inhomo_poisson_sampler(tgrid, f.vec[, clust[i]])$sample
  sample.mat = rbind(sample.mat, cbind(i, samples[[i]]))
  sample.len[i] = length(samples[[i]])
}

save(tgrid, f.vec, samples, clust, sample.mat, sample.len, file = paste(workpath, "data/lgcp_simu_data_v0.RData", sep = "/"))

# ------------------------------------------------------------------
# Junk code:
#col.vec = c("blue", "red", "green")
#pch.vec = c(rep(0, 6), rep(1, 7), rep(2, 7))#c(1 : 6, 1 : 7, 1 : 7)
#cex.vec = c(rep(0.8, 6), rep(1, 7), rep(0.8, 7))
#windows()
#plot(tgrid, f.vec[, 1], type = "l", lwd = 2, col = col.vec[1], ylim = c(-4.2, 5.2), xlab = "Genomic location", ylab = "")
#lines(tgrid, f.vec[, 2], type = "l", lwd = 2, col = col.vec[2])
#lines(tgrid, f.vec[, 3], type = "l", lwd = 2, col = col.vec[3])
#for (i in 1 : nsample)
#{
#  points(samples[[i]], rep(-i * 0.2, sample.len[i]), col = col.vec[clust[i]], pch = pch.vec[i], cex = cex.vec[i]) #pch = paste(pch.vec[i]), cex = 0.6
#}
#