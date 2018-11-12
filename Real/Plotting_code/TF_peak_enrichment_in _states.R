#################################### code to generate enrichment matrix of TF peaks in chromatin states. 
## rows correspond to states
##columns correspond to Transcription factors

library(miceadds)
library(igraph)
library(caTools)
library(foreach)
library(doParallel)

## parallel processing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

### combine all mm10 bed files converted using Liftover for all TF into a list
filelist<-list.files(paste0(workpath,"Real/data/"),(pattern = ".txt"))
combined.TF.chrom<-vector("list",length(filelist))
combined.TF.regionmeans<-vector("list",length(filelist))
combined.TF.region<-vector("list",length(filelist))
strnames<-unlist(lapply(filelist, function(x) strsplit(x, "_")[[1]][1]))
names(combined.TF.regionmeans)<-strnames
names(combined.TF.chrom)<-strnames

for (i in 1:length(filelist)){
  TF_final <- read.table(paste0(workpath,"Real/data/",filelist[i]), sep='\t', quote=NULL, comment='', header=FALSE)
  names(TF_final)<-c("Chromosome","Start","End")
  TF.regionmeans<-(TF_final$Start+TF_final$End)/2
  combined.TF.chrom[[i]]<-as.character(TF_final$Chromosome)
  combined.TF.regionmeans[[i]]<-TF.regionmeans
  combined.TF.region[[i]]<-TF_final  
}

nsample<- length(combined.TF.regionmeans)

## dihmm
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
state_list<-list("Broad Promoter"=rbind(temp$D30, temp$D23 , temp$D13 , temp$D25, temp$D28 , temp$D6 , temp$D22), "Bivalent Promoter" = rbind(temp$D16 , temp$D24 , temp$D10 , temp$D27 , temp$D26), 
                 "Super Enhancer"=rbind(temp$D18), "Poised Enhancer"=temp$D3, "Upstream Enhancer" = rbind(temp$D2,temp$D17,temp$D12), "Transcribed" = rbind(temp$D4 , temp$D15 , temp$D21), "Boundary"=temp$D14,
                 "Low Signal" = rbind(temp$D1 ,temp$D5),"Low Coverage" = rbind(temp$D7 , temp$D29 , temp$D8, temp$D20), "Polycomb Repressed" = rbind(temp$D9 , temp$D11 , temp$D19))

## get TF peaks across states function
get_TF_peaks_in_states <- function(state_list_j,combined.TF.regionmeans,combined.TF.chrom,nsample,state_list_j_names){
  temp=vector('list',nsample)
  for (j in 1:dim(as.data.frame(state_list_j))[1]){
    s_j=state_list_j$start[j]
    e_j=state_list_j$end[j]
    chrom_j<-state_list_j$chr[j]
    ## split to get chromosome and region means 
    index_region = lapply(combined.TF.regionmeans,function(x) which(unlist(x)>s_j & unlist(x)< e_j))
    index_chr<-lapply(combined.TF.chrom,function(x) which(x==chrom_j))
    index_common=(mapply(intersect, index_chr, index_region))
    if(any(sapply(index_common,length)>0)){
      index_common<-index_common[which(sapply(index_common,length)>0)]
      iidx<-which(names(combined.TF.regionmeans)%in% names(index_common))
      for(kk in 1:length(iidx)){
        temp[[iidx[kk]]]<-c(temp[[iidx[kk]]],combined.TF.regionmeans[[iidx[kk]]] [as.numeric(unlist(index_common[kk]))])
      }
    }
  }
  names(temp)<-strnames
  return (temp)
}

## get TF_peaks for each state
TFpeaks_domains<-vector('list',length(state_list))
TFpeaks_domains<-foreach (j = 1:length (state_list),.combine = 'list',.multicombine = TRUE,.maxcombine = length(state_list),.final = function(x) setNames(x, names(state_list))) %dopar% {
  get_TF_peaks_in_states(state_list[[j]],combined.TF.regionmeans,combined.TF.chrom,nsample, names(state_list)[j])
}

## save states
save(TFpeaks_domains,file = paste0(workpath,"Real/results/TF_peak_enrichment_for_heatmap.RData"))

########################### code to generate heatmap ###################################
## generate table for TF peaks in each state
temp=lapply(TFpeaks_domains, function(x) lapply(x, function(y) length(y)))
temp=as.data.frame(do.call('rbind',temp))
## total number of TF peaks in all states
apply(temp, 2, function(x) sum(unlist(x)))

## enrichment matrix code
TF_D_matrix=t(data.matrix(temp, rownames.force = NA))
D_n=length(state_list)
###denominator for enrichment: M/N
## get window length for all states
all_width=unlist(lapply(state_list, function(x) sum(x[,3]-x[,2]+1)))
den=all_width/sum(all_width) #M/N
num=TF_D_matrix/rowSums(TF_D_matrix) #m/n
TF_norm=t(t(num)/den)
TF_norm_log=10*log10(TF_norm+0.1)
## using default pheatmap
annotation_col = data.frame(State=colnames(TF_norm_log))
rownames(annotation_col)<-colnames(TF_norm_log)
Var1= c("firebrick1", "deeppink1","gold1","orange3","yellow","chartreuse4","lightskyblue","gray48","mintcream","navajowhite1","seagreen1")
names(Var1) <- rownames(annotation_col)
anno_colors <- list(Var1 = Var1)
#colnames(annotation_col)<-'State'
library(RColorBrewer)
library(pheatmap)
pheatmap(TF_norm_log, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq_len(10))),annotation_colors=anno_colors,annotation_col=data.frame(annotation_col),clustering_distance_cols = "euclidean",clustering_method = "complete",
         main='Log enrichment of TFs in states',show_colnames = F,legend_labels = "Enrichment Value")

## using complex heatmap
library(ComplexHeatmap)
library(pvclust)
library(MASS)
library(circlize)
type = colnames(TF_norm_log)
ha = HeatmapAnnotation(df = data.frame(type = type),col = list(type = c("Broad Promoter"="red3", "Bivalent Promoter"="orchid1","Super Enhancer" ="gold1","Poised Enhancer"="orange3","Upstream Enhancer"="yellow",
                                                                        "Transcribed" = "chartreuse4","Boundary" ="lightskyblue","Low Signal" ="mintcream", "Polycomb Repressed" ="gray48","Low Coverage" ="seagreen1")))
Heatmap(TF_norm_log, name = "enrichment", col = colorRamp2(c(min(TF_norm_log),0,max(TF_norm_log)), c("cornflowerblue", "yellow", "red")),top_annotation = ha, top_annotation_height = unit(4, "mm"), show_row_names = TRUE, 
        show_column_names = FALSE,column_title = 'Log enrichment of TFs in states',column_title_gp = gpar(fontsize = 20, fontface = "bold"), clustering_distance_columns = "spearman",clustering_method_columns = "single")

