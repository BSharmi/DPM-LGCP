#rm(list = ls())
#library(INLA)
#INLA:::inla.dynload.workaround() 
library(miceadds)
library(igraph)
library(caTools)
library(foreach)
library(doParallel)

## parallel processing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

##workpath = "/home/bsharmi6/R_scripts/SharmiWork/"

## call function to compute clustering across a state within different randomly selected windows
compute_cluster <- function(Dstate,fpath){
  library(INLA)
  #INLA:::inla.dynload.workaround() 
  library(miceadds)
  library(igraph)
  library(caTools)
  library(foreach)
  library(doParallel)
  
  ## load TFpeaks_domains data for each state 
  file_name=paste0(workpath,"Real/results/NSC_TFpeaks_diHMM_domain_",Dstate,".RData")
  load(file_name)  
  
  ## get TF peak locations. First list element is the merged or shifted peak locations
  TFpeaks_each_domain_coniditonal_k=TFpeaks_each_domain_coniditonal[[1]]
  ## remove TFs with one peak
  TFpeaks_each_domain_coniditonal_k=TFpeaks_each_domain_coniditonal_k[unlist(unlist(lapply(TFpeaks_each_domain_coniditonal_k, function(x) length(x)>1)))]
  
  TF_in_cluster=names(TFpeaks_each_domain_coniditonal_k)
  nsample_k<-length(TF_in_cluster)
  
  ## set parameters
  sigma0 = 0.9809434
  kappa0 = 1.131607
  tau0 = 0.196793
  L = 200
  #L=max(unlist(lapply(TFpeaks_each_domain_coniditonal_k, function(x) length(x))))
  sample.mat = NULL
  sample.len<-c()
  ## scale points
  min_val=min(unlist(TFpeaks_each_domain_coniditonal_k))
  max_val=max(unlist(TFpeaks_each_domain_coniditonal_k))
  samples=lapply(TFpeaks_each_domain_coniditonal_k, function(x) L *(x-min_val)/(max_val-min_val))
  for(ii in 1:length(samples)){
    sample.mat = rbind(sample.mat, cbind(ii, samples[[ii]]))
    sample.len[ii] = length(samples[[ii]]) 
  }
  tgrid = seq(0, L, length = 1000)
  #tgrid = seq(min(unlist(TFpeaks_each_domain_coniditonal_k)), max(unlist(TFpeaks_each_domain_coniditonal_k)), length = 1000)
  nknots = 700 #35
  mesh_spde = get_mesh_spde(nknots = nknots, range(tgrid), sigma0, kappa0, tau0)
  mesh = mesh_spde$mesh
  spde = mesh_spde$spde
  alpha = 1
  n.iter = 20 # 200
  n.burn = 0
  init.method = "random"
  init.k = 3
  timeupdate = 5
  system.time({
    res = Neal_alg3_lgcp(alpha, n.iter, n.burn, samples, sample.mat, sample.len, range(tgrid), init.method, init.k, timeupdate, mesh, spde, verbose = FALSE)           
  }) # 2151 sec
  
  nMc = nrow(res$s_record)
  pairs = matrix(NA, nrow = nMc, ncol = choose(nsample_k, 2))
  adjmat = array(NA, dim = c(nsample_k, nsample_k, nMc)) 
  for (i in 1 : nMc)
  {
    si = res$s_record[i, ]  
    pairs[i, ] = s_to_vec(si, nsample_k)
    adjmat[, , i] = vec_to_mat(pairs[i, ], nsample_k)
  }
  
  init.vec = s_to_vec(res$initial$s0, nsample_k)
  init.mat = vec_to_mat(init.vec, nsample_k)
  #colnames(init.mat)<-names(npos1_new_k[iidx])
  rownames(init.mat)<-TF_in_cluster
  colnames(init.mat)<-TF_in_cluster
  g.init = graph.adjacency(init.mat, mode = "undirected")
  est.vec = colMeans(pairs) > 0.5
  est.mat = vec_to_mat(est.vec, nsample_k)
  rownames(est.mat)<-TF_in_cluster
  colnames(est.mat)<-TF_in_cluster
  #estmat[[k]] = est.mat
  g.est = graph.adjacency(est.mat, mode = "undirected")
  
  #format graph nodes and vertices
  V(g.est)$label.cex = 0.8
  V(g.est)$size = 20
  E(g.est)$width = 0.9
  
  # Estimate the intensity for each cluster
  clust.est = vec_to_s(est.vec, nsample_k)
  new.ns = table(clust.est)
  new.sample = pool_list(sample.mat, sample.len, clust.est)
  mean.zvalue = est_inten_list_v2(new.sample, new.ns, range(tgrid), mesh, spde, mlik1 = FALSE, config1 = FALSE, zgrid.eval = TRUE, tgrid = tgrid, save_result = FALSE, take_sample = FALSE, verbose = FALSE)
  
  # save estamt, g.est
  gsymbol_g.est <- list(est.mat, g.est,mean.zvalue,res)  
  #return(est.mat)
  #return(gsymbol_g.est)
  
  ## save results
  save(gsymbol_g.est, file = paste0(workpath,"Real/results/diHMM_NSC_res_TFcluster_",Dstate,".RData"))
}

# load state names
Dstate_all<-get(load(paste0(workpath,"Real/results/diHMM_state_names.RData")))

# run clustering for different chromatin states
foreach (i_state = 1 : length(Dstate_all), .verbose=TRUE) %dopar% {
  compute_cluster(Dstate_all[i_state],fpath)}  
