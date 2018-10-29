
## get TF peaks in different diHMM windows based on a threshold number of peaks
get_TF_peaks_in_states <- function(state_list_ij,combined.TF.regionmeans,combined.TF.chrom,combined.TF.region){
  s_j=state_list_ij$start
  e_j=state_list_ij$end
  chrom_j<-unique(state_list_ij$chr) ## since based on chromosome
  index_region = lapply(combined.TF.region, function(qq) which(apply(qq, 1, function(x) any(x[1] == chrom_j & as.numeric(as.character(x[4])) > s_j & as.numeric(as.character(x[4])) < e_j))))
  index_region<-index_region[which(sapply(index_region,length)>0)]
  iidx<-which(names(combined.TF.regionmeans)%in% names(index_region))
  temp=mapply(function(x, y) x[y],combined.TF.regionmeans[iidx],index_region,SIMPLIFY = FALSE)
  temp_region=mapply(function(x, y) x[y,],combined.TF.region[iidx],index_region,SIMPLIFY = FALSE)
  return (list(temp,state_list_ij, temp_region))
}


######################## main code #############################
library(miceadds)
library(igraph)
library(caTools)
library(foreach)
library(doParallel)

## parallel processing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)


## read chromatin state files by diHMM; remove title
dihmm_file <- read.table(paste0(workpath,"Real/data/NSC_nD30_nB30_domainLevelStatesColor.bed"), sep='\t', quote=NULL, comment='', header=FALSE, skip=1)
## select only required columns
dihmm_file<-dihmm_file[,1:4]
names(dihmm_file)<-c('chr', 'start', 'end', 'domains')
## add a peak length column
dihmm_file$length=dihmm_file$end-dihmm_file$start
## get state names
states=unique(dihmm_file$domains)

## store all state details in a list
temp<-split(x = dihmm_file,f = dihmm_file$domains)

## split domain level states by chromosome
dom_list_chr=lapply(temp, function(l) split(x=l, f=l$chr)) ##each domain state and by chromosome

## update state names to match same order as dom_list_chr
states=names(dom_list_chr)

## save states
save(states, file=paste0(workpath,"Real/results/diHMM_state_names.RData"))

## save state_list 
save(dom_list_chr, file = paste0(workpath,"Real/results/diHMM_domain_states.RData"))


## read TF peak files
filelist<-list.files(paste0(workpath,"Real/data/"),(pattern = "peaks"))
combined.TF.chrom<-vector("list",length(filelist))
combined.TF.regionmeans<-vector("list",length(filelist))
combined.TF.region<-vector("list",length(filelist))
strnames<-unlist(lapply(filelist, function(x) strsplit(x, "_")[[1]][1]))
names(combined.TF.regionmeans)<-strnames
names(combined.TF.chrom)<-strnames
names(combined.TF.region)<-strnames
for (i in 1:length(filelist)){
  TF_final <- read.table(paste0(workpath,"Real/data/",filelist[i]), sep='\t', quote=NULL, comment='', header=FALSE)
  TF_final[,4]=(TF_final[,2]+TF_final[,3])/2
  names(TF_final)<-c("Chromosome","Start","End","Avg")
  TF.regionmeans<-(TF_final$Start+TF_final$End)/2
  combined.TF.chrom[[i]]<-as.character(TF_final$Chromosome)
  combined.TF.regionmeans[[i]]<-TF.regionmeans
  combined.TF.region[[i]]<-TF_final  
}

## no of TFs
nsample<- length(combined.TF.regionmeans)

### get TF_peaks for each state
for (i_state in 1:length(dom_list_chr)){
  ## first get peaks for each TF in 20 chromosomes for each state
  chr_indx=which(unlist(lapply(dom_list_chr[[i_state]], function(x) (dim(x)[1])>0)))
  ## list of TF peaks for all chromosomes
  TFpeaks_each_domain_coniditonal<-foreach (j = 1 : length(chr_indx),.combine = 'list',.multicombine = TRUE,.maxcombine = length(chr_indx),.final = function(x) setNames(x, names(chr_indx))) %dopar% {
    get_TF_peaks_in_states(dom_list_chr[[i_state]][[chr_indx[j]]],combined.TF.regionmeans,combined.TF.chrom,combined.TF.region)
  }
  
  ### now merge TF peaks for all chromosomes for each state
  
  ## peak locations of all TFs in 20 chromosomes
  all_TF_peaks=sapply(TFpeaks_each_domain_coniditonal, '[[', 1) ## eliminates one level of list
  ## remove empty chormosome list
  all_TF_peaks=all_TF_peaks[which(unlist(lapply(all_TF_peaks, function(x) length(x)>0)))]
  
  ## sort peaks locations increasing to decreasing
  all_TF_peaks=lapply(all_TF_peaks, function(x) lapply(x, function(y) sort(y)))
  
  ## peak windows of all TFs in 20 chromosomes
  all_TF_peak_windows=sapply(TFpeaks_each_domain_coniditonal, '[', 3) ## eliminates one level of list
  ## remove empty chromosomes
  all_TF_peak_windows=all_TF_peak_windows[which(unlist(lapply(all_TF_peak_windows, function(x) length(x)>0)))]
  
  ## state_windows of all TFs in 20 chromosomes
  all_state_windows=sapply(TFpeaks_each_domain_coniditonal, '[', 2) ## eliminates one level of list
  ## remove  unmatched from all_TF_peaks
  all_state_windows=all_state_windows[intersect(names(all_state_windows),names(all_TF_peaks))]
  
  ## chr list obtained from mm10 unique chromosomes. 
  chr_list=c(paste0("chr",seq(1:19)), "chrX","chrY")
  
  ## sort TF peaks list by chromosome
  all_TF_peaks=all_TF_peaks[match(chr_list, names(all_TF_peaks),nomatch = 0)]
  
  ## sort TF peak windows list by chromosome
  all_TF_peak_windows=all_TF_peak_windows[match(chr_list, names(all_TF_peak_windows),nomatch = 0)]
  
  ## sort state windows list by chromosome
  all_state_windows=all_state_windows[match(chr_list, names(all_state_windows),nomatch = 0)]
  
  ## add cumulative max value of last peak location form all TFs of previous chromosome to current chromosome
  all_TF_peaks_shifted=list()
  all_TF_peaks_shifted[1]=all_TF_peaks[1]## first chromosome no change
  for (i in 2:length(all_TF_peaks)){
    append_value=max(unlist(all_TF_peaks_shifted[i-1])) ## max peak location on previous chromosome
    all_TF_peaks_shifted[i]= lapply(all_TF_peaks[i], function(x) lapply(x, function(y) y+append_value))
  }
  names(all_TF_peaks_shifted)=names(all_TF_peaks)
  
  ## merge TF peaks for all chromosomes
  TF_names=unique(unlist(c(sapply(all_TF_peaks_shifted,names))))
  all_TF_peaks_merged=vector('list',length(TF_names))
  names(all_TF_peaks_merged)=TF_names
  for(ii in 1:length(TF_names)){
    iidx=which(unlist(lapply(all_TF_peaks_shifted, function(x) any(names(x) %in% names(all_TF_peaks_merged)[ii]))))
    all_TF_peaks_merged[[ii]]=unname(unlist(lapply(all_TF_peaks_shifted[iidx], function(x) x[which(names(x) %in% names(all_TF_peaks_merged)[[ii]])])))
  }
  ## save file. each state has one file. each file is a list of even numbers. odd number is the TFs list, even number is the state window data frame
  TFpeaks_each_domain_coniditonal=list(all_TF_peaks_merged,all_state_windows,all_TF_peak_windows,all_TF_peaks)
  save(TFpeaks_each_domain_coniditonal,file=paste0(workpath,'Real/results/NSC_TFpeaks_diHMM_domain_',names(dom_list_chr)[i_state],'.RData'))  
}

