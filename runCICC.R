###############################################
## CICC pipeline for one sample
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################

## #####################################
## arguments to scripts
args <- commandArgs(TRUE)
nbOutlierstoRemove <- toString(args[1]) ## number of outlier methods to remove (set to 1)
IDS <- toString(args[2]) ## PCAWG/Simulation ID of the sample to run on
DIRNAME <- toString(args[3]) ## output directory (will be created through mkdir if does not exist)
## #####################################
## libraries to load
library(BiocGenerics,lib="~/R/library/")
library(S4Vectors,lib="~/R/library/")
library(IRanges,lib="~/R/library/")
library(GenomeInfoDb,lib="~/R/library/")
library(GenomicRanges,lib="~/R/library/")
## #####################################
## set working directory
setwd("~/")
source("CICC.PAC.R") ## source CICC functions
source("loadData.R") ## source functions for input and output processing
dyn.load("scoringlite.so") ## source dynamic C library for creation of hard-assignment matrix
## #####################################
## tries to create output directory if does not exist
try(system(paste0("mkdir ",DIRNAME)),silent=T)
## #####################################


## #####################################
## names of the methods (loadData.R contains reading functions with hard coded paths to the data and takes the method names as input)
methods=c("cloneHD",
          "DPClust",
          "phylogic",
          "CCube",
          "pyclone",
          "sclust",
          "CliP",
          "CTPsingle",
          "BayClone",
          "phylowgs",
          "svclone"
          )
## #####################################


## #####################################
## loading and saving data for all methods to run on
print("loading data")
lAA <- lapply(IDS,function(x) try(loadAllMethods(x,methods=methods),silent=T))
save(lAA,file=paste0(DIRNAME,"/",
                     "lAA.nontransformed.",nbOutlierstoRemove,".Rda"))
## remove outliers (here: 1 outlier max and 0 outlier min based on how many an the fraction of mutations they report on
print("removing outliers")
lAA2 <- transformlAA(lAA,
                     nbOutliers=nbOutlierstoRemove,
                     downsamplers=NULL
                     )
names(lAA2) <- IDS
save(lAA2,file=paste0(DIRNAME,"/","lAA.transformed.",nbOutlierstoRemove,".Rda"))
## #####################################


## #####################################
print("running CICC")
## #####################################
## run CICC on preprocessed input and save results
system.time(allRes <- lapply((1:length(lAA2)),function(x)
{
    print(x)
    cicc <- try(consensusMatrix(lAA2[[x]],
                                pMethods=c(1),
                                repeats=100,
                                x),silent=T)
    cicc
}))
names(allRes) <- names(lAA2)
save(allRes,file=paste0(DIRNAME,"/allResClusts.out",nbOutlierstoRemove,".Rda"))
## #####################################
## merge clusters that are "too close" to each other
## this step is not performed anymore: stands as a dummy for potential later re-use
print("merging step")
allRes2 <- lapply(1:length(allRes),function(x)
{
    try(if(F){clusts <- allRes[[x]]$clusts
            mergedClusts <- mergeClusters(lAA2,clusts,x,p=0.05)
            return(append(allRes[[x]],list(mergedClusts=mergedClusts)))},silent=T)
    return(try(append(allRes[[x]],list(mergedClusts=allRes[[x]]$clusts)),silent=T))
})
names(allRes2) <- names(lAA2)
save(allRes2,file=paste0(DIRNAME,"/allResClusts2.out",nbOutlierstoRemove,".Rda"))
## #####################################
## writes the outputs: mutation assignments and subclonal structures in CP and CCF
print("write results")
writeResultsMA(allRes2,lAA2,DIRNAME)
writeResultsClusterCCF(allRes2,lAA2,DIRNAME)
writeResultsClusterCCF2(allRes2,lAA2,DIRNAME)
## #####################################


## #####################################
q(save="no")
## #####################################
