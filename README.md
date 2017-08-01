# DPM-LGCP
Identifying Transcriptional Regulatory Modules among Different Chromatin States 
# Overview
DPM-LGCP stands for Dirichlet Prior mixture for Log Gaussian Cox Process. It is a non-parametric Bayesian clustering technique based on non-homogeneous Poisson process which has been used to identify transcriptional regulatory modules in distinct chromatin states. The algorithm clusters transcription factors sharing similar intensity function (binding patterns). It has broader applications to cases requiring un-supervised clustering of large datasets with prior information on clusters.
# System Requirements
DP-LGCP is independent of operating systems as it is written in R. Basic requirements for running the pipeline include installing R and the following libraries- 'INLA'(INLA 0.0-1468872408, dated 2016-07-18 (14:43:05+0100)), 'miceadds', 'foreach', 'doParallel', 'igraph', 'caTools'.
# Usage
Unzip the package. Change the current directory in R to the 'GitHub_code' folder containing the code and data organized into subfolders - 'data' contains diHMM generated domain level bed files and mm10 three column peak files of transcription factors to be used in clustering step;  'main' contains all R scripts ; 'results' contains RData files.

The following steps need to be followed to obtain TF regulatory modules in different chromatin states - 

1. In R console, run the command below where 'workpath' is the full file path to the 'GitHub_code' directory, e.g.


```
> workpath = "/Users/sharmibanerjee/Documents/Summer2017/GitHub_code/"
> source(paste0(workpath,"/main/script_generate_TF_peaks_in_diHMM_states.R"))
> get_TF_peaks_in_states_conditional_threshold(workpath)
```
This step reads diHMM and peak files and selects those TFs in each chromatin state window that contain peaks greater than the threshold. After step 1 completes, RData files should be generated in the "results" folder.

2. Run the following commands on R console to generate clustering results
```
> source(paste0(workpath,"/main/script_get_clusters_for_each_diHMM_state.R"))
> get_clusters_for_each_chromatin_state(workpath)
```

The RData files from clustering step should be generated in the 'results' folder with the names 'diHMM_celltype_res_TFcluster_threshold_chromatin_state_name'
