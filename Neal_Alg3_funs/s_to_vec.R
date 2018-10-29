#...........function 11: transfer the cluster assignment s to a vector of size (N-1)*N/2.........
# vec.s is a vector of an upper triangular matrix (by rows) indicating whether sample i and j belong to the same
# cluster. 

s_to_vec = function(s, N){
  
  vec.s = rep(0, (N-1)*N/2)
  tag  = 1 
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (s[i]==s[j]){
        vec.s[tag] = 1
      }
      tag = tag + 1
    }
  }
  
  return(vec.s)
}


vec_to_s = function(vec, N){
  
  s = rep(NA, N)
  s[1] = 1
  
  tag  = 0 
  nclust = 1
  for (i in 1:(N-1)){
 # for (i in 1:6){
    for (j in (i+1):N){	
	#for (j in 9:N){	
        tag = tag + 1	     	
	    if (vec[tag]==1){
			if (any(!is.na(s[i])||!is.na(s[j]))){
				#id1= c(which(!is.na(s[i])), which(!is.na(s[j])))
				if (!is.na(s[i])){				 
					s[j] = s[i]
				} else if (!is.na(s[j])){
				    s[i] = s[j]
				}
		    } else {
				s[i] = nclust+1	
				s[j] = nclust+1	
				nclust = nclust+1
			}
		}		
	  }	     
    }  
  s[which(is.na(s))] = nclust+1;  
  
  return(s)
}

#...........function 11: transfer the vec.s in function 11 to a matrix.........
vec_to_mat = function(vec.s, N){
  mat.s = matrix(0, nrow=N, ncol=N)
  tag = 1
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      mat.s[i,j] = vec.s[tag]
      mat.s[j,i] = mat.s[i,j]
      tag = tag + 1
    }
  }
  
  return(mat.s)
}
