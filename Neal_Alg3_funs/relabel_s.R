relabel_s = function(s0){
  #  s0 is 1 by N vector, the cluster labels of the N samples. 
  #  to avoid keep giving higher labels, we will do the relabelling. 
  
  uni.val = unique(s0) # these unique values are ordered according to the original z.all values. 
  n.u = length(uni.val)  
  s.new = rep(NA, length(s0))
  
  for (i in 1:n.u){
    idi = which(s0==uni.val[i])        
    s.new[idi] = i
  }
  
  return(s.new)
}