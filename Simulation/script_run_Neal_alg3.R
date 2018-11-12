# run Neal 1998 algorithm 3 for clustering
# ----------------------------------------
#rm(list = ls())
library(INLA)
#INLA:::inla.dynload.workaround() 
library(miceadds)
library(igraph)
library(caTools)

## change work path as needed. It should include full path upto LGCP_package
workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/"
load(paste(workpath, "Simulation/data/lgcp_simu_data_v0.RData", sep = "/"))
source.all(path = paste(workpath, "Neal_Alg3_funs", sep = "/"), print.source = FALSE)


# Get mesh and spde setup
nclust = ncol(f.vec)
param = matrix(NA, nclust, 4)
colnames(param) = c('sigma0', 'kappa0', 'tau0', 'range0')
for (i in 1 : nclust)
{
  outi = get_spde_param_with_true_lam(f.vec[, i], tgrid, d = 1, 0)
  param[i, ] = c(outi$sigma0, outi$kappa0, outi$tau0, outi$range0)
}
sigma0 = min(param[, 1])
kappa0 = min(param[, 2])
tau0 = min(param[, 3])

nknots = 35
mesh_spde = get_mesh_spde(nknots = nknots, range(tgrid), sigma0, kappa0, tau0)
mesh = mesh_spde$mesh
spde = mesh_spde$spde

alpha = 1
n.iter = 10 #200
n.burn = 0
init.method = "random"
init.k = 6
timeupdate = 5
system.time({
  res = Neal_alg3_lgcp(alpha, n.iter, n.burn, samples, sample.mat, sample.len, range(tgrid), init.method, init.k, timeupdate, mesh, spde, verbose = FALSE)					 
}) 

#-------------------------------------------------------------- 
# now look at the clustering result. 
# ----------------------------------------------------------

nsample = length(samples)
nMc = nrow(res$s_record)
pairs = matrix(NA, nrow = nMc, ncol = choose(nsample, 2))
adjmat = array(NA, dim = c(nsample, nsample, nMc)) 
for (i in 1 : nMc)
{
  si = res$s_record[i, ]  
  pairs[i, ] = s_to_vec(si, nsample)
  adjmat[, , i] = vec_to_mat(pairs[i, ], nsample)
}

true.vec = s_to_vec(clust, nsample)
true.mat = vec_to_mat(true.vec, nsample)
g.true = graph.adjacency(true.mat, mode = "undirected")
init.vec = s_to_vec(res$initial$s0, nsample)
init.mat = vec_to_mat(init.vec, nsample)
g.init = graph.adjacency(init.mat, mode = "undirected")
est.vec = colMeans(pairs) > 0.5
est.mat = vec_to_mat(est.vec, nsample)
g.est = graph.adjacency(est.mat, mode = "undirected")

V(g.true)$label.cex = 0.8
V(g.init)$label.cex = 0.8
V(g.est)$label.cex = 1.2
V(g.true)$size = 18
V(g.init)$size = 18
V(g.est)$size = 20
E(g.true)$width = 0.8
E(g.init)$width = 0.8
E(g.est)$width = 1.2

save(param, mesh_spde, res, pairs, adjmat, file = paste(workpath, "Simulation/results/lgcp_simu_result_v0.RData", sep = "/"))

# Plot the (marginal) posterior mean estimate in terms of graphs
par(mfrow = c(1, 3))
plot(g.true, main = "True Cluster")
plot(g.init, main = "Initial Cluster", vertex.color = "white")
plot(g.est, main = "Estimated Cluster", vertex.color = "green")
#dev.off()
# calculate the statistic
print(paste("initial error rate: ", mean(abs(init.vec - true.vec)), sep = ""))
print(paste("estimation error rate: ", mean(abs(est.vec - true.vec)), sep = ""))

# Plot true intensity functions vs. estimated intensity functions within each cluster
# vs. intensity function of the posterior mode

## True intensity
par(mfrow = c(1, 1))
plot(tgrid, f.vec[, 1], type = "l", col = "black", lwd = 2)
lines(tgrid, f.vec[, 2], col = "green", lwd = 2)
lines(tgrid, f.vec[, 3], col = "red", lwd = 2)
legend(6,3,c('lam1', 'lam2','lam3'),col=c('black','green','red'), lty=c(1,1,1),cex=0.5)

# Estimate the intensity for each cluster
clust.est = vec_to_s(est.vec, nsample)
new.ns = table(clust.est)
new.sample = pool_list(sample.mat, sample.len, clust.est)
mean.zvalue = est_inten_list_v2(new.sample, new.ns, range(tgrid), mesh, spde, mlik1 = FALSE, config1 = FALSE, zgrid.eval = TRUE, tgrid = tgrid, save_result = FALSE, take_sample = FALSE, verbose = FALSE)

col.vec = c("blue", "red", "green")
# pch.vec = c(1 : 6, 1 : 7, 1 : 7)
plot(tgrid, exp(mean.zvalue[1, ]), type = "l", lwd = 2, col = col.vec[1],  ylim = c(min(exp(mean.zvalue)), max(exp(mean.zvalue))),
     xlab = "location", ylab = "")
lines(tgrid, exp(mean.zvalue[2, ]), type = "l", lwd = 2, col = col.vec[2])
lines(tgrid, exp(mean.zvalue[3, ]), type = "l", lwd = 2, col = col.vec[3])
legend(6.5, 5, c("Estimated intensity of cluster 1", "Estimated intensity of cluster 2", "Estimated intensity of cluster 3"), col = col.vec, lty = c(1, 1, 1), cex = 0.7)




# Looking at result
like1 = log_likelihood_given_c_lgcp(samples, clust.est, range(tgrid), mesh, spde)
marg.lik.all.est = like1$l # overall marginal likelihood based on estimated clusters
marg.lik.cluster.est = like1$marg_like # marginal likelihood of each cluster

like2 = log_likelihood_given_c_lgcp(samples, clust, range(tgrid), mesh, spde)
marg.lik.all.true = like2$l  # overall marginal likelihood based on true cluster
marg.lik.cluster.true = like2$marg_like

print(paste("overall marginal likelihood: ", marg.lik.all.est, " (estimated); ", marg.lik.all.true, " (true)", sep = ""))
print(paste("marginal likelihood per cluster: ", toString(marg.lik.cluster.est), " (estimated); ", toString(marg.lik.cluster.true), " (true)", sep = ""))

