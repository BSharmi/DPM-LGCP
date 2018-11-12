# Chromatin state specific transcriptional regulatory network identification

Version	1.0

Date		11/12/2018

Description	DPM-LGCP stands for Dirichlet Prior mixture for Log Gaussian Cox Process. It is a non-parametric Bayesian clustering technique based on non-homogeneous Poisson process which has been used to identify transcriptional regulatory modules in distinct chromatin states. The algorithm clusters transcription factors sharing similar intensity function (binding patterns). It has broader applications to cases requiring un-supervised clustering of large datasets with prior information on clusters.	
Author 	Sharmi Banerjee, Honxiao Zhu, Man Tang, Xiaowei Wu, Wu Feng, David Xie	

## System Requirements
DP-LGCP is independent of operating systems as it is written in R. Basic requirements for running the pipeline include installing R and the following libraries- ‘miceadds’, ‘foreach’, ‘doParallel’, ‘igraph’, ‘caTools’, ‘ComplexHeatmap’, ‘pvclust’, ‘MASS’, ‘circlize’, ‘RColorBrewer’, ‘pheatmap’. INLA package in provided here. Please see description below on how to add the package to the R libraries.

## Description of the folders

- INLA contains the R INLA library
- Real contains script, data and some results for neural stem cells used in the paper 
  
  -Data – peak files for twenty-one transcription factors, diHMM generated domain bed files for thirty domains and mm10 refseq genes. The basic requirements towards using the peipeline are 
   
   1. Chromatin state identified file generate by a genome segmentation software. In this case we are using Domain-level state files generated by diHMM. The file name is under Real/Data/NSC_nD30_nB30_domainLevelStatesColor.bed. The first 4 columns from this file are read which are: chromosome, start and end positions and the label of the Domain. On usage of diHMM tool, please see - https://github.com/gcyuan/diHMM
   
   2. Peak locations of the transcription factor. These files are 3-column bed/txt files with three columns: chromosome, start and end positions
  
  -Main_Code – scripts to generate clusters in chromatin states
  
  -Plotting_Code -  scripts to generate enrichment maps of TF peaks and TF modules in chromatin states
  
  -Results – clustering results (RData files) on chromatin states

- Simulation has code and results for simulation data

  -code – scripts to simulate transcription factors binding sites with twenty TFs and three clusters and evaluate results by applying clustering algorithm using non-homogeneous Poisson point process
  
  -data – simulated data generated from scripts in ’code’ folder
  
  -results – clustering results
  
  -Neal_Alg3_funs’ has the main algorithm clustering functions.



### Setting INLA path 
Download and copy INLA folder under LGCP-package to the R library of the system. For example for a R 3.3 version, the path is  - /Library/Frameworks/R.framework/Versions/3.3/Resources/library

### Clustering on simulation data
1.	Set working directory in R to the path up to the Simulation directory in LGCP package. E.g. -
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Simulation/')
2.	Set a variable 'workpath' in R console. E.g. on R console -
workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Simulation/code/”
3. Source the functions needed
source.all('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Neal_Alg3_funs/')
4.	Run simulation script -
source("script_simulate_data.R")
5.	Evaluate clustering result -
source("script_run_Neal_alg3.R")

### Simulation case study – 
•	Three 1D non-homogeneous Poisson processes with Log-Gaussian intensity were generated. Next, binding site locations of twenty transcription factors were simulated by drawing independent samples from these non-homogeneous Poisson processes. 

•	Using the simulated binding site locations, we applied the clustering algorithm. The TFs were randomly assigned to any cluster at initializations. The results are shown after ten iterations –
  - initial error rate: 0.39
  - estimation error rate: 0
  - overall true marginal likelihood: 76.73
  - overall estimated marginal likelihood: 76.73
  -	true marginal likelihood per cluster: 27.27, 14.95, 34.49
  -	estimated marginal likelihood per cluster: 27.27, 14.95, 34.49

### Clustering on real data
1.	Set working directory in R to the path up to the Real directory in LGCP package. E.g. -\
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Real/Main_Code/')
2.	Set a variable 'workpath' in R console. E.g. on R console-\
workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/”
3. Source the functions needed\
source.all('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Neal_Alg3_funs/')
4.	Run the following scripts –

    a. source("script_generate_TF_peaks_in_diHMM_states.R")\
This step generates a results (in .RData) format containing the TF peaks contained within each chromatin state

    b. source("script_get_clusters_for_each_diHMM_state.R")\
This step generates clustering results within each chromatin state using results from step (a)

5.	Files generated 

    a.	Clustering results will be stored in the ‘results’ folder under ‘Real’ with the file names diHMM_celltype_cluster_domainname.RData
    
    b.	TF peaks in each chromatin state will be generated with names celltype_TFpeaks_diHMM_chromatinstate.RData

### Plotting enrichment heat-maps
1.	Set working directory in R to the path up to the Real directory in LGCP package. E.g. -
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Real/Plotting_Code/')
2.	Set work path in R console. E.g. -
workpath = “/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/”
3.	Run the following scripts to generate heat map of TF peaks within different chromatin states
source("TF_peak_enrichment_in _states.R")


## Contact
bsharmi6@vt.edu
