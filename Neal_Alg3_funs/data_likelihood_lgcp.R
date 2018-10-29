data_likelihood_lgcp = function(x, Ns=1, T, mesh1, spde1, cc){

if (length(cc)>1){
	#x_sum = colSums(X)
    #n_s = dim(X)[1]
    #l = sum(lgamma(x_sum+G0alphabeta[1,]) + lgamma(n_s - x_sum + G0alphabeta[2,]) - lgamma(G0alphabeta[1,]) - lgamma(G0alphabeta[2,])
    #    + lgamma(colSums(G0alphabeta)) - lgamma(colSums(G0alphabeta)+n_s))      
} else if (length(cc)==1){
	
	if (is.list(x)==FALSE){ # this is the case when x is one point process. 
		l = lgcp_inla_brief(x, Ns, T, mesh1, spde1, mlik1=TRUE, verbose1=FALSE)	
	} else { # this is the case when x is a list of multiple point processes. 
		l = lgcp_inla_brief(x[[1]], Ns, T, mesh1, spde1, mlik1=TRUE, verbose1=FALSE)				
	}
	# if (is.vector(X)){
		# l = sum(lgamma(X + G0alphabeta[1,]) + lgamma(1-X+G0alphabeta[2,])- lgamma(G0alphabeta[1,]) - lgamma(G0alphabeta[2,]) +
			# lgamma(colSums(G0alphabeta)) - lgamma(G0alphabeta[1,]+G0alphabeta[2,]+1)) # done	
	# }else{	
		# l = sum(lgamma(X[1,] + G0alphabeta[1,]) + lgamma(1-X[1,]+G0alphabeta[2,])- lgamma(G0alphabeta[1,]) - lgamma(G0alphabeta[2,]) +
			# lgamma(colSums(G0alphabeta)) - lgamma(G0alphabeta[1,]+G0alphabeta[2,]+1)) # done	
	# }	
}

return(l)
}