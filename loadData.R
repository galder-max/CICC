###############################################
## CICC helper
## set of functions to load the input data from disk
## and write the output results on disk
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################


###############################################
load. <- function(sim,method="truth",header=T)
{
    path <- "/srv/shared/vanloo/ICGC-clustering/final_clustering/"
    if(method=="BayClone")
        path <- paste0(path,"BayClone/Ji_ICGC_output_2721_2017-03-02_final_run/mutation_assignment/",
                       sim,"_mutation_assignments.txt")
    else if(method=="CCube")
        path <- paste0(path,"CCube/lustre/fmlab/yuan03/PCAWG/db_consensus_final_012017/results/ccube_v0.3_final_new/",sim,"/",
                       sim,"_mutation_assignments.txt.gz")
    else if(method=="DPClust")
        path <- paste0(path,"dpclust/2_subclones/",
                       sim,"_mutation_assignments.txt.gz")
    else if(method=="pyclone")
        path <- paste0(path,"pyclone/lustre/fmlab/yuan03/PCAWG/db_consensus_final_012017/results/pyclone_v0.4_final/",
                       sim,"/",sim,"_mutation_assignments.txt.gz")
    else if(method=="cloneHD")
        path <- paste0(path,"cloneHD/",
                       sim,".assignment_probability_table.txt")
    else if(method=="phylogic")
        path <- paste0(path,"phylogic/broad_snv_clustering_submission/",
                       sim,"_mutation_assignments.txt")
    else if(method=="phylowgs")
        path <- paste0(path,"phylowgs/consensus.pwgs.runf.combined/",
                       sim,"_mutation_assignments.txt.gz")
    else if(method=="CliP")
        path <- paste0(path,"CliP/consensus_format/",
                       sim,"_mutation_assignments.txt.gz")
    else if(method=="CTPsingle")
        path <- paste0(path,"CTPsingle/summary_final_run_CTPsingle/",sim,"/",
                       "",
                       sim,"_mutation_assignments.txt.gz")
    else if(method=="sclust")
        path <- paste0(path,"sclust/Sclust_Mar_2017_final/",sim,"_mutation_assignments.txt")
    else if(method=="svclone")
        path <- paste0(path,"svclone2/final_v2/SNVs/",sim,"_assignment_probability_table.txt")
    else
        stop("Methods not in table")
    return(read.table(path,header=header))
}

load.multiplicities <- function(sim,method="truth",header=T)
{
    path <- "/srv/shared/vanloo/ICGC-clustering/final_clustering/"
    if(method=="BayClone")
        path <- paste0(path,"BayClone/Ji_ICGC_output_2721_2017-03-02_final_run/multiplicity/",
                       sim,"_multiplicity.txt")
    else if(method=="CCube")
        path <- paste0(path,"CCube/lustre/fmlab/yuan03/PCAWG/db_consensus_final_012017/results/ccube_v0.3_final_new/",sim,"/",
                       sim,"_multiplicity.txt.gz")
    else if(method=="DPClust")
        path <- paste0(path,"dpclust/0_multiplicity/",
                       sim,"_multiplicity.txt.gz")
    else if(method=="pyclone")
        path <- paste0(path,"pyclone/lustre/fmlab/yuan03/PCAWG/db_consensus_final_012017/results/pyclone_v0.4_final/",
                       sim,"/",sim,"_multiplicity.txt.gz")
    else if(method=="cloneHD")
        path <- paste0(path,"cloneHD/",
                       sim,".assignment_probability_table.txt")
    else if(method=="phylogic")
        path <- paste0(path,"phylogic/broad_snv_clustering_submission/",
                       sim,"_multiplicity.txt")
    else if(method=="phylowgs")
        path <- paste0(path,"phylowgs/consensus.pwgs.runf.combined/",
                       sim,"_multiplicity.txt.gz")
    else if(method=="CliP")
        path <- paste0(path,"CliP/consensus_format/",
                       sim,"_multiplicity.txt.gz")
    else if(method=="CTPsingle")
        path <- paste0(path,"CTPsingle/summary_final_run_CTPsingle/",sim,"/",
                       "",
                       sim,"_multiplicity.txt.gz")
    else if(method=="sclust")
        path <- paste0(path,"sclust/Sclust_Mar_2017_final/",sim,"_multiplicity.txt")
    else if(method=="svclone")
        path <- paste0(path,"svclone2/final_v2/SNVs/",sim,"_multiplicity.txt")
    else
        stop("Methods not in table")
    return(suppressWarnings(read.table(path,header=header)))
}

getMultiplicity <- function(sim, methods, ids)
{
    allM <- lapply(methods[methods!="cloneHD"],function(x)
    {
        k <- try(load.multiplicities(sim,x),silent=T)
        if(x=="CliP") colnames(k) <- c("chr","pos","tumour_copynumber","multiplicity")
        if(inherits(k, "try-error")) return(NULL)
        if(is.null(k)) return(NULL)
        k
    })
    allM <- sapply(allM,function(x)
    {
        if(is.null(x)) return(rep(NA,length(ids)))
        nms <- paste(x[,"chr"],x[,"pos"],sep=":")
        ww <- (!duplicated(nms))
        x <- x[ww,]
        nms <- paste(x[,"chr"],x[,"pos"],sep=":")
        rownames(x) <- nms
        x[ids,"multiplicity"]
    })
    floor(apply(allM,1,median,na.rm=T))
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
    nMaj <- inverse.rle(list(lengths=CN,values=cn[,"major_cn"]))
    nMin <- inverse.rle(list(lengths=CN,values=cn[,"minor_cn"]))
    all <- rep(NA,length(keepSNV))
    all[keepSNV] <- nMaj+nMin
    return(all)
}

getCCFs <- function(lAA)
{
    sampID <- names(lAA)[1]
    cn <- read.table(paste0("~/ICGC/CNA/",
                            sampID,
                            ".consensus.20170119.somatic.cna.annotated.txt"),
                     header=T)
    purities <- read.table(paste0("~/CICC/",
                                  "consensus.20170217.purity.ploidy.txt"),
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
    ids <- names(lAA[[1]][[1]])
    snv <- snv[ids,]
    mult <- getMultiplicity(sim=sampID,
                            methods=names(lAA[[1]]),
                            ids=ids)
    cn_cl <- getCN(cn,snv)
    alt_count <- as.numeric(gsub("(.*);t_alt_count=(.*?);(.*)","\\2",as.character(snv[,8])))
    ref_count <- as.numeric(gsub("(.*);t_ref_count=(.*?);(.*)","\\2",as.character(snv[,8])))
    tumor_f <- alt_count/(alt_count+ref_count)
    CCFs <- tumor_f/purity*(purity*cn_cl+(1-purity)*2)/mult
    return(list(CCFs=CCFs,purity=purity))
}

loadTable <- function(sim,method)
{
    if(method=="cloneHD")
    {
        t <- load.(sim,method,header=F)
        ##print(head(t))
        if(ncol(t)==3)
            ass <- rep(1,nrow(t))
        else
            ass <- apply(t[,-c(1,2)],1,which.max)
        names(ass) <- paste(t[,1],t[,2],sep=":")
    }
    else if(method=="svclone")
    {
        t <- load.(sim,method,header=T)
        ##print(head(t))
        if(ncol(t)==3)
            ass <- rep(1,nrow(t))
        else
            ass <- apply(t[,-c(1,2)],1,which.max)
        names(ass) <- paste(t[,1],t[,2],sep=":")
    }
    else if(method%in%c("CCube","DPClust","pyclone",
                        "phylogic","phylosub","phylowgs",
                        "sclust",
                        "truth",
                        "CliP","BayClone","CTPsingle"))
    {
        t <- load.(sim,method,header=T)
        ass <- t[,3]
        names(ass) <- paste(t[,1],t[,2],sep=":")
    }
    return(ass)
}

loadAllMethods <- function(sim,
                           methods=c("cloneHD",
                                     "DPClust",
                                     "phylogic",
                                     "CCube",
                                     "pyclone",
                                     "phylowgs",
                                     "svclone",
                                     "sclust",
                                     "CliP",
                                     "CTPsingle",
                                     "BayClone"
                                     ))
{
    allass <- lapply(methods,function(x)
    {
        tab <- try(loadTable(sim,x),silent=T)
        if(inherits(tab,"try-error"))
            return(NULL)
        tab
    })
    wNull <- which(!sapply(allass,is.null))
    allass <- lapply(wNull,function(x) allass[[x]])
    names(allass) <- methods[wNull]
    nms <- unique(unlist(lapply(allass,names)))
    aa <- lapply(allass,function(x)
    {
        x <- x[nms]
        names(x) <- nms
        x
    })
    passNAcounts <- rowSums(!sapply(aa,is.na))>=round(length(methods)/2+1)
    aa <- lapply(aa,function(x) x[passNAcounts])
    return(aa)
}

writeResultsMA <- function(allRes,lAA,dirName)
{
    for(x in names(lAA))
    {
        try(if(T)
            {
                t.2c <- cbind(gsub("(.*):(.*)","\\1",
                                   names(lAA[[x]][[1]])),
                              gsub("(.*):(.*)","\\2",
                                   names(lAA[[x]][[1]])),
                              allRes[[x]]$mergedClusts)
                colnames(t.2c) <- c("chr","pos","cluster")
                write.table(t.2c,
                            file=gzfile(paste0(dirName,"/",
                                               x,".mutation.assignment.txt.gz")),
                            quote=F,col.names=T,row.names=F)
            },silent=T)
    }
}

writeResultsClusterCCF <- function(allRes,lAA, dirName)
{
    ii <- names(lAA)[1]
    {try(if(T){
        ma <- allRes[[1]]$mergedClusts
        ccfs <- getCCFs(lAA)
        ccfs <- sapply(unique(ma),function(x)
        {
            median(ccfs$CCFs[ma==x],na.rm=T)
        })
        nbs <- sapply(unique(ma),function(x)
        {
            sum(ma==x,na.rm=T)
        })
        tt <- cbind(unique(ma),nbs,ccfs)
        colnames(tt) <- c("cluster","n_ssms","proportion")
        tt <- tt[order(tt[,3],decreasing=T),]
        tt <- rbind(tt,c(NULL,NULL,NULL))
        write.table(tt,file=gzfile(paste0(dirName,"/",
                                          ii,".cluster.ccfs.txt.gz")),
                    quote=F,col.names=T,row.names=F)},silent=T)
    }
}

writeResultsClusterCCF2 <- function(allRes,lAA, dirName)
{
    ii <- names(lAA)[1]
    {try(if(T){
        ma <- allRes[[1]]$mergedClusts
        ccfs <- getCCFs(lAA)
        ccfs <- sapply(unique(ma),function(x)
        {
            median(ccfs$CCFs[ma==x]*ccfs$purity,na.rm=T)
        })
        nbs <- sapply(unique(ma),function(x)
        {
            sum(ma==x,na.rm=T)
        })
        tt <- cbind(unique(ma),nbs,ccfs)
        colnames(tt) <- c("cluster","n_ssms","proportion")
        tt <- tt[order(tt[,3],decreasing=T),]
        tt <- rbind(tt,c(NULL,NULL,NULL))
        write.table(tt,file=gzfile(paste0(dirName,"/",
        ii,".cluster.cps.txt.gz")),
        quote=F,col.names=T,row.names=F)},silent=T)
    }
}
###############################################



