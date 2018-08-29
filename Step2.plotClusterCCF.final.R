######################################################
## plotting script for furhter manual inspections
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################

###############################################
args <- commandArgs(TRUE) ## no args
######################################################
DIRNAME <- paste0("/srv/shared/vanloo/ICGC-consensus-clustering/CICC-final/")
setwd(DIRNAME)
######################################################
IDS <- dir(DIRNAME)
######################################################
METHODS <- c("DPClust",
             "phylogic",
             "CCube",
             "pyclone",
             "sclust",
             "CliP",
             "svclone",
             "phylowgs",
             "CTPsingle",
             "BayClone")
######################################################


######################################################
source(paste0("loadData.R"))
######################################################
library(parallel)
library(BiocGenerics,lib="~/R/library/")
library(S4Vectors,lib="~/R/library/")
library(IRanges,lib="~/R/library/")
library(GenomeInfoDb,lib="~/R/library/")
library(GenomicRanges,lib="~/R/library/")
## #####################################


######################################################
getCCFs <- function(sampID,ids,method)
{
    cn <- read.table(paste0("~/ICGC/CNA/",
                            sampID,
                            ".consensus.20170119.somatic.cna.annotated.txt"),
                     header=T)
    purities <- read.table(paste0("~/CICC/",
                                  "consensus.20170119.purity.ploidy.txt"),
                           header=T)
    purity <- purities[purities[,1]==sampID,"purity"]
    snv <- try(read.table(paste0("~/ICGC/snv_mnv/",
                             sampID,
                             ".consensus.20160830.somatic.snv_mnv.vcf.gz"),
                      header=F),silent=T)
    if(inherits(snv,"try-error"))
        snv <- read.table(paste0("/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/graylist/snv_mnv/",
                        sampID,
                        ".consensus.20160830.somatic.snv_mnv.vcf.gz"),
                        header=F)
                       rns<-paste(snv[,1],snv[,2],sep=":")
                        rem<-!duplicated(rns)
                        snv<-snv[rem,]
    rownames(snv) <- paste(snv[,1],snv[,2],sep=":")
    snv <- snv[ids,]
    alt_count <- as.numeric(gsub("(.*);t_alt_count=(.*?);(.*)","\\2",as.character(snv[,8])))
    ref_count <- as.numeric(gsub("(.*);t_ref_count=(.*?);(.*)","\\2",as.character(snv[,8])))
    tumor_f <- alt_count/(alt_count+ref_count)
    mult <- getMultiplicity(sim=sampID,
                            methods=c(method,method,method),
                            ids=ids)
    cn_cl <- getCN(cn,snv)
    CCFs <- tumor_f/purity*(purity*cn_cl+(1-purity)*2)/mult
    CCFs[is.infinite(CCFs)]<-NA
    return(list(CCFs=CCFs,purity=purity,cn=cn))
}


getCN <- function(cn,snv)
{
    snvOrd <- cbind(as.character(snv[,1]),as.character(snv[,2]))
    snvOrd[,1][snvOrd[,1]=="X"] <- "23"
    snvOrd[,1][snvOrd[,1]=="Y"] <- "24"
    snvOrd <- apply(snvOrd,2,as.numeric)
    orderSNV <- order(snvOrd[,1],snvOrd[,2],decreasing=F)
    snv <- snv[orderSNV,]
    ## ##############################
    require(GenomicRanges)
    gr1 <- GRanges(snv[,1],IRanges(snv[,2],snv[,2]))
    gr2 <- GRanges(cn[,1],IRanges(cn[,2],cn[,3]))
    snvCN <- subsetByOverlaps(gr1,
                              gr2)
    keepSNV <- (gr1%in%snvCN)
    CN <-countOverlaps(gr2,
                       snvCN)
    snv <- snv[keepSNV,]
    nMaj <- inverse.rle(list(lengths=CN,values=cn[,"major_cn"]))
    nMin <- inverse.rle(list(lengths=CN,values=cn[,"minor_cn"]))
    return(nMaj+nMin)
}


plotClusterCCF <- function(sampID,method,clusts)
{
    ccfs <- getCCFs(sampID,names(clusts),method)
    ccfsC <- sapply(unique(clusts[!is.na(clusts)]),function(x)
    {
        median(ccfs$CCFs[clusts==x],na.rm=T)
    })
    nbs <- sapply(unique(clusts[!is.na(clusts)]),function(x)
    {
        sum(clusts==x,na.rm=T)
    })
    brks <- seq(0,max(ccfs$CCFs,na.rm=T)+0.02,0.01)
    par(mar=c(3,3,3,1))
    hist(ccfs$CCFs,breaks=brks,main=paste(method),xlab="CCF distribution",cex.main=.7)
    cols <- rainbow(length(ccfsC))
    abline(v=ccfsC,lwd=3,col=cols)
    mtext(side=3,paste0("purity=",ccfs$purity),cex=.5)
    legend("topleft",
           col=cols,
           cex=.5,
           pch=rep(19,length(ccfsC)),
           as.character(paste(signif(ccfsC,3),nbs,sep="; ")))
}

plotCCALL<-function(sampID)
{
    load(paste0(DIRNAME,"/",sampID,"/lAA.transformed.1.Rda"))
    par(mfcol=c(3,3))
    retV=lapply(names(lAA2[[1]]),function(method)
        if(!method%in%c("cloneHD","phylowgs")) try(plotClusterCCF(sampID,method,lAA2[[1]][[method]]),silent=T))
}
######################################################



######################################################
mclapply(IDS,function(sampID)
{
    pdf(paste0(DIRNAME,"/",sampID,"/",sampID,".hist.clustCCF.all.pdf"))
    k<-try(plotCCALL(sampID),silent=T)
    dev.off()
},mc.cores=1)
######################################################


######################################################
q(save="no")
######################################################
