###############################################
## pipeline to submit all CICC runs for PCAWG
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################


###############################################
## reads in a file with all PCAWG sample IDs
IDS <- as.character(read.table("~/CICC/consensus.20170119.purity.ploidy.txt",header=T)[,1])
## output directory where output directories for each sample will be created
DIRNAME <- "/srv/shared/vanloo/ICGC-consensus-clustering/CICC-final/"
## number of outliers to remove before running CICC
NBOUT <- 1
###############################################


###############################################
## libraries to load
library(BiocGenerics,lib="~/R/library/")
library(S4Vectors,lib="~/R/library/")
library(IRanges,lib="~/R/library/")
library(GenomeInfoDb,lib="~/R/library/")
library(GenomicRanges,lib="~/R/library/")
###############################################




######### create submit file for one sample for SLURM
generateSubmitFile<-function(nbOut, sampID, dirout,file)
{
    cat("#!/bin/bash",file=file)
    cat("\n",file=file,append=T)
    cat("## a cluster job for CICC",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH -n 1",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH -t 12:00:00 ",file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("#SBATCH --job-name=",sampID,sep=""),file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH --mem=14G",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("## execute the CICC pipeline",file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("module load R"),
        file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("module load R-bundle-Bioconductor/3.2-foss-2016a-R-3.2.3"),
        file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("Rscript runCICC.R",
              nbOut,#1
	      sampID, #2
	      dirout, #3
              sep=" "),
        file=file,append=T)
    cat("\n",file=file,append=T)
    system(paste("chmod u+rwx ",file,sep=""))
}

callBash<-function(file)
{
    print(file)
    print(system(paste("sbatch ",file," -e ~/err.txt &",sep="")))
}


submitJob <- function(id)
{
    submitFile <- paste0("~/submitFiles/",id,".final.bash")
    generateSubmitFile(nbOut=NBOUT,
                       sampID=id,
                       dirout=paste0(DIRNAME,"/",id),
                       file=submitFile)
    callBash(submitFile)
}
###############################################


###############################################
## submits one job per sample
kNull <- lapply(IDS,submitJob)
###############################################
q(save="no")
###############################################



