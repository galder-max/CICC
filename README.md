---
author:
- Maxime Tarabichi
title: 'ClusterID-based consensus clustering'
---
# CICC
## Description and example run
A description of the concepts and example runs with dummy data is available in the pdf and the description folder (corresponding knitr Rnw).
## Installation 
There is no installation required, the R scripts that can be run on their own but do have some dependencies.
Compile the C code in the scripts directory with the following command:
> R CMD SHLIB scoringlite.c

This will generate the dynamic library to make co-clustering matrices from hard assignments vectors.
### Dependencies
R packages:
BiocGenerics
S4Vectors
IRanges
GenomeInfoDb
GenomicRanges
## Run - step1
run _Step1.submitALL.R_ with:
>Rscript Step1.submitALL.R

This pipeline goes through the PCAWG sample IDs and submit one job per sample on a slurm-based cluster to run CICC.  
The job will run _runCICC.R_, which loads the required data using utility functions from _loadData.R_, where paths are encoded, and then runs CICC and saves the output in a consensus format using functions from _loadData.R_. 
## Run - step2
run _Step2.plotClusterCCF.R_ with:
> Rscript Step2.plotClusterCCF.R

This is useful for visualisation of the results. It plots histograms of CCF and the cluster positions for each method and for the consensus.
