
Neal_alg3_lgcp <- function(alpha, n.iter=100, n.burn = 10, samples, Y_p, n_all, T, init.method='random', 
                            init.k = 3, timeupdate = 2, mesh1, spde1, verbose=FALSE){
  # This is the Gibbs sampler following Alg3 of Neal 1998 "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
  # This algorithm is described in section 2.2 of Hal Daume III's paper on DPM search. The algorithm 3 does not require sampling
  # theta_star, but it also uses a polya urn sampler on marginal likelihood. 
  
  # Input:
  # n.iter: total number of MCMC iterations.
  # n.burn: number of burnin period to throw away.
  # samples: the point process data in the form of a list.
  # Y_p: the point process data in the form of a matrix, col1: the sample number, col2: the points.
  # n_all: a vector giving the number of points (of each point process) in the samples. 
  # T: a vector of two component, giving the lower and upper bound of the 1-D domain.
  # z.all: the initial values of log(lamb(s)), evaluated on the knots locations. defaul=NULL. 
  # init.method: the method to use to create the initial cluster. default='random'
  # init.k: the initial number of clusters in the initial clustering. 
  # nknots: number of knots (equal distant grid) on which we create mesh,spde, and evaluate z(s) on.
  # timeupdate: the cycle by which we print on screen the iteration id. 
  
  # Output:
  # post: a list that record all posterior clustering, including s (cluster assignment), 
  #       nj (size of each cluster), z.star (unique values of z*(s) for each cluster) 
  # pairs: a vector of size N*(N-1)/2, a vectorized 1/0 vector indicating whether
  #       a pair of (i,j) belongs to the same cluster.
  # nlist: a matrix recording the size of each cluster during the MCMC, the size is sorted from large to small.
  # klist: a vector indicating the number of clsuters found during the MCMC.
  # mesh: the common mesh object used in inla.
  # spde: the common spde object used in inla.
  
  ## add library to fix NA problem
  library(matrixStats)



  N = length(samples)
  
  #create initial clustering s
  if (init.k == 1){
    s = rep(1, N)
    lg_inten1 = NULL
  } else if (init.method=='random'){
    s = sample(1:init.k, N, replace=TRUE, prob=rep(1/init.k,times=init.k))
    if (length(unique(s))<init.k){
      s = relabel_s(s)
      init.k = length(unique(s))      
    }
    lg_inten1 = NULL
  } else if (init.method=='hc'){
    lg_inten1 = est_inten_list_v2(samples, Ns=1, T, mesh1, spde1, mlik1=FALSE, config1=FALSE,
                              zgrid.eval=FALSE, tgrid=NULL, save_result=FALSE)
    s = two_stage_clustering(lg_inten1, init.k)$s
  }else{ # This is the case when z.all is prespecified.
    s = get_unique_rows(z.all)$s
    lg_inten1 = NULL
  }   
 
  s_all = rep(s, times=n_all)  
  initial = list(lginten = lg_inten1, s0 = s, init.method = init.method, knots=mesh1$loc)
  cat(paste('....cluster initialized using method', init.method, '....\n', sep=" "))
  
  # initialize marginals  
  marginals = rep(NA,times=N)
  for (i in 1 : N){
  # zhu: logH(x_n), n=1,...,N.
	marginals[i] = log_marginal_posterior_lgcp(samples[[i]], Ns=1, T, mesh1, spde1, xs=NULL)    
	if (i%%timeupdate==0){
	cat(paste('i=',i),'\n')
	}  
  }
  cat(paste('...the marginal likelihood calculated...'),'\n')  
  
  
  #---  start the Gibbs sampler -----------
  s_record = matrix(NA, nrow = n.iter-n.burn, ncol = N) # this records the cluster assignment at every iteration, the s vector.
  njlist <- matrix(NA, nrow = n.iter-n.burn, ncol = 30)  # record sizes of 30 largest clusters
  klist <- rep(NA, times = n.iter-n.burn)   # record the number of clusters in each simulation
  pairs = matrix(NA, nrow = n.iter-n.burn, ncol = N*(N-1)/2) # vectorized pairwise matrix (upper triangular, by row) of two point-processes belonging to the same cluster.
   
  # create a matrix to save H({x_i, i in I}) for arbitrary set of I, with I
  logH_all = cbind(diag(rep(1,N)),marginals) # the first N columns save whether i in I (1 vs 0),  
  # the last component saves the log(H({x_i, i in I}) )
  
  t0 = proc.time()
  for(iter in 1:n.iter){
    
    if (iter%%timeupdate==0) cat(paste('iter =',iter, '\n'))
    
    for (i in 1:N){
	    if (i%%(round(N/2))==0){
			cat(paste('iter=',iter,'; i=',i),'\n')
		}  
		X_i = samples[[i]]	
		s_minus = s[-i];
		nj=table(s_minus)
		sj=as.numeric(names(nj))
		
		# calculate the conditional Hj = H(xi | x_{k,k!=i,c_k=j})
	    Hj = rep(NA, length(nj)) # all calculated in log scale.
		for (j in 1:length(nj)){
		    list_ij = (s_all==sj[j] & Y_p[,1]!=i) 
			X_k =Y_p[list_ij,2] # all points in X_{k,k!=i,c_k=j}, excluding the ith sample
			
			# the following are for the purpose of saving some previous calculated H(x_{i, i in I}), to save computation. 
			id_k = unique(Y_p[list_ij,1])  
			vec_k = rep(0, N)
			vec_k[id_k] = 1
						
			X_joint =c(X_i,X_k);
			id_joint = c(i, id_k);
			vec_joint = rep(0, N)
			vec_joint[id_joint] = 1
			
			find_joint_id = rowSums(abs(t(vec_joint-t(logH_all[,1:N]))))
			if (any(find_joint_id==0)){
			   Hij_joint = logH_all[which(find_joint_id==0),N+1]
			}else{
			   Hij_joint = lgcp_inla_brief(X_joint, Ns=nj[j]+1, T, mesh1, spde1, mlik1=TRUE, verbose1=FALSE)
			   logH_all = rbind(logH_all,c(vec_joint, Hij_joint)); # save Hij_joint and Hij_k for future use.	
			}
				
			find_k_id = rowSums(abs(t(vec_k-t(logH_all[,1:N]))))
			if (any(find_k_id==0)){
			   Hij_k = logH_all[which(find_k_id==0),N+1]
			}else{
			   Hij_k =  lgcp_inla_brief(X_k, Ns=nj[j], T, mesh1, spde1, mlik1=TRUE, verbose1=FALSE);
			   logH_all = rbind(logH_all,c(vec_k, Hij_k)); # save Hij_joint and Hij_k for future use.	
			}
			Hj[j] = Hij_joint - Hij_k								
		}

    ## change by Man
    logpj <- c(Hj+log(nj),marginals[i]+log(alpha))  
    logsumpj=logSumExp(logpj)
    s_i <- sample(c(sj,max(sj)+1), 1, replace=FALSE, prob = exp(logpj-logsumpj)) 
		 
		#pj <- exp(c(Hj+log(nj),marginals[i]+log(alpha)))                        # p(s[i]=j | ...), j=1..k, k+1
    #s_i <- sample(c(sj,max(sj)+1), 1, replace=FALSE, prob = pj/sum(pj))           # sample s[i]
        # max(sj)+1 means we have created a new cluster. since this new cluster is relative to sj or s[-i], 
		# we may ended up having higher and higher labels (and the lower labels may disapear). so everytime
		# we will need to relabel the s. it won't influence our clustering pattern, only influnce the labeling.
		
		s[i] = s_i 
		s = relabel_s(s)	
        s_all = rep(s, times=n_all)		
	}  
    
    # save related posterior samples
    if (iter > n.burn){
      s_record[iter-n.burn,] = s      
      njlist[iter-n.burn, 1:length(unique(s))] = sort(table(s),decreasing =TRUE)
      klist[iter-n.burn] = length(unique(s))
      pairs[iter-n.burn,] = s_to_vec(s, N)
    }
  }
  time1 = proc.time() - t0 
  
  return(list(s_record=s_record, pairs=pairs, nlist=njlist, klist = klist, mesh=mesh1, spde = spde1, 
              time=time1, initial = initial, marginals=marginals))
}
