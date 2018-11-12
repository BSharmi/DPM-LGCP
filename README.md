# Chromatin state specific transcriptional regulatory network identification

Version	1.0

Date		11/12/2018

Description	DPM-LGCP stands for Dirichlet Prior mixture for Log Gaussian Cox Process. It is a non-parametric Bayesian clustering technique based on non-homogeneous Poisson process which has been used to identify transcriptional regulatory modules in distinct chromatin states. The algorithm clusters transcription factors sharing similar intensity function (binding patterns). It has broader applications to cases requiring un-supervised clustering of large datasets with prior information on clusters.	
Author 	Sharmi Banerjee, Honxiao Zhu, Man Tang, Xiaowei Wu, Wu Feng, David Xie	

## System Requirements
DP-LGCP is independent of operating systems as it is written in R. Basic requirements for running the pipeline include installing R and the following libraries- 'INLA'(INLA 0.0-1468872408, dated 2016-07-18 (14:43:05+0100)), ‘miceadds’, ‘foreach’, ‘doParallel’, ‘igraph’, ‘caTools’, ‘ComplexHeatmap’, ‘pvclust’, ‘MASS’, ‘circlize’, ‘RColorBrewer’, ‘pheatmap’.

## Usage
Unzip the package. Change the current directory in R to the ‘LGCP_package' folder containing the code and data organized into subfolders 	

- INLA contains the R package
- Real contains script, data and some results for neural stem cells used in the paper 
  
  -Data – peak files for twenty-one transcription factors, diHMM generated domain bed files for thirty domains and mm10 refseq genes. 
  
  -Main_Code – scripts to generate clusters in chromatin states
  
  -Plotting_Code -  scripts to generate enrichment maps of TF peaks and TF modules in chromatin states
  
  -Results – clustering results (RData files) on chromatin states

- Simulation has code and results for simulation data

  -code – scripts to simulate transcription factors binding sites with twenty TFs and three clusters and evaluate results by applying clustering algorithm using non-homogeneous Poisson point process
  
  -data – simulated data generated from scripts in ’code’ folder
  
  -results – clustering results
  
  -Neal_Alg3_funs’ has the main algorithm clustering functions.



### Setting INLA path 
Extract INLA compressed file and copy INLA folder under LGCP-package to R library located in - /Library/Frameworks/R.framework/Versions/3.3/Resources/library

### Clustering on simulation data
1.	Set working directory in R to the path up to the Simulation directory in LGCP package. E.g. -
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Simulation/')
2.	Set work path in R console. E.g. -
workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Simulation/code/”
3.	Run simulation script -
source("script_simulate_data.R")
4.	Evaluate clustering result -
source("script_run_Neal_alg3.R")

### Simulation case study – 
•	Three 1D non-homogeneous Poisson processes with Log-Gaussian intensity were generated. Next, binding site locations of twenty transcription factors were simulated by drawing independent samples from these non-homogeneous Poisson processes. In Figure 1A, the solid curves denote the true intensities. 
•	Using the simulated binding site locations, we applied the clustering algorithm. The TFs were randomly assigned to any cluster at initializations. The results are shown after ten iterations –
o	initial error rate: 0.39
o	estimation error rate: 0
o	overall true marginal likelihood: 76.73
o	overall estimated marginal likelihood: 76.73
o	true marginal likelihood per cluster: 27.27, 14.95, 34.49
o	estimated marginal likelihood per cluster: 27.27, 14.95, 34.49

### Clustering on real data
1.	Set working directory in R to the path up to the Real directory in LGCP package. E.g. -
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Real/Main_Code/')
2.	Set work path in R console. E.g. -
workpath = "/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/”
3.	Run the following scripts –
source("script_generate_TF_peaks_in_diHMM_states.R")
source("script_get_clusters_for_each_diHMM_state.R")
4.	Files generated 
a.	Clustering results will be stored in the ‘results’ folder under ‘Real’ with the file names diHMM_celltype_cluster_domainname.RData
b.	TF peaks in each chromatin state will be generated with names celltype_TFpeaks_diHMM_chromatinstate.RData

### Plotting enrichment heat-maps
1.	Set working directory in R to the path up to the Real directory in LGCP package. E.g. -
setwd('/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/Real/Plotting_Code/')
2.	Set work path in R console. E.g. -
workpath = “/Users/sharmibanerjee/Documents/Summer2017/LGCP_package/”
3.	Run the following scripts to generate heat map in Figure 2(a) in main paper
source("TF_peak_enrichment_in _states.R")


## Contact
bsharmi6@vt.edu
