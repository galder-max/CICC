---
author:
- Maxime Tarabichi
title: 'ClusterID-based consensus clustering'
---
# CICC
## Description
A description of the concepts is available in the pdf and the description folder (corresponding knitr Rnw).
## Run - step1
>run _Step1.submitALL.R_
This pipeline goes through the PCAWG sample IDs and submit one job per sample on a slurm-based cluster to run CICC.  
The job will run _runCICC.R_, which loads the required data using utility functions from _loadData.R_, where paths are encoded, and then runs CICC and saves the output in a consensus format using functions from _loadData.R_.
## Run - step2
>run _Step2.plotClusterCCF.R_
This is useful for visualisation of the results. It plots histograms of CCF and the cluster positions for each method and for the consensus.
