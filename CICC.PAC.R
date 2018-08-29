###############################################
## CICC functions
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################


## pairwise similarity between two strings of assignment to clusters
## each string corresponds to one mutation's vector of assignment
## each cluster in the string corresponds to one method's cluster assignment
## NOTE: NA assignments count as dissimilarities within a method.
votedist <- function(g1,g2)
{
    suppressWarnings(if(T){
    g1 <- as.numeric(strsplit(g1,split="-")[[1]])
    g2 <- as.numeric(strsplit(g2,split="-")[[1]])
    return(sum(g1==g2,na.rm=T))})
}

## total distance matrix between unique vector of assignments
## total number of vector of methods (nMethods) - similarity matrix (m)
votedist.matrix <- function(allgs,nMethods)
{
    m <- matrix(0,length(allgs),length(allgs))
    for(i in 1:(length(allgs)-1))
        for(j in (i+1):length(allgs))
        {
            m[i,j] <- votedist(allgs[i],allgs[j])
            m[j,i] <- votedist(allgs[i],allgs[j])
        }
    rownames(m) <- allgs
    colnames(m) <- allgs
    nMethods-m
}

## distance matrix in mutation number : should have no impact on the final clustering
votedist.matrix.nbmut <- function(allgs,tA)
{
    m <- matrix(NA,length(allgs),length(allgs))
    for(i in 1:(length(allgs)-1))
        for(j in (i+1):length(allgs))
            m[j,i] <- m[i,j] <- -max(tA[allgs[i]],tA[allgs[j]])
    rownames(m) <- allgs
    colnames(m) <- allgs
    m
}

## concatenate cluster votes
getClsA <- function(allA)
{
    matA <- t(sapply(allA,function(x) x))
    clsA <- NULL
    for(i in 1:nrow(matA))
        clsA <- paste(clsA,matA[i,],sep=if(i==1) "" else "-")
    clsA
}

## performs the consensus clustering
## allA is a list of vectors of cluster assignment
## there is one vector of assignment per method
## they have the same length and share the same order of mutations
fastConsensusClustering <- function(allA,
                                    nbClusters,
                                    keepnames)
{
    nMut <- length(allA[[1]])
    nMethods <- length(allA)
    matA <- t(sapply(allA,function(x) x))
    clsA <- NULL
    for(i in 1:nrow(matA))
        clsA <- paste(clsA,matA[i,],sep=if(i==1) "" else "-")
    vuA <- sort(table(clsA),decreasing=T)
    dist1 <- votedist.matrix(names(vuA),nMethods)*max(vuA)*10
    dist2 <- votedist.matrix.nbmut(names(vuA),vuA)
    distF <- as.dist(dist1+dist2)
    hc <- hclust(distF,method="ward.D")
    lClusts <- lapply(nbClusters,function(nC) cutree(hc,k=nC))
    return(lapply(lClusts,function(x){
        clusts <- x[clsA]
        names(clusts) <- keepnames
        clusts
    }))
}

## derives a co-clustering matrix from a hard-assignment vector (v) in C.
makeHardAss <- function(v)
{
    matrix(.C("hardass",
              as.integer(length(v)),
              as.integer(v),
              as.integer(rep(0,length(v)*length(v))))[[3]],
           length(v),length(v))
}

## add one permutation to the consensus co-clustering matrix
CCM <- function(ccm,
                allCC,
                ii,
                jj,
                repeats,
                i,
                representants)
{
    for(j in 1:repeats)
    {
        v <- allCC[[ii]][[jj]][[j]][[i]][representants]
        v[is.na(v)] <- (max(v,na.rm=T)+1):(max(v,na.rm=T)+sum(is.na(v)))
        ccm <- ccm+makeHardAss(v)
        gc()
    }
    ccm
}

## choose optimal K given the consensus clustering matrices
## ccm is a list of consensus clustering matrices from which the PAC is derived
chooseOptimalK <- function(ccm,medianK,weights)
{
    if(length(ccm)==2) return(2)
    ## from Dr. Yasin Şenbabaoğlu:
    ## shenbaba.weebly.com/blog/how-to-use-the-pac-measure-in-consensus-clustering
    Kvec = 2:medianK
    x1 = 0.1; x2 = 0.9 ## threshold defining the intermediate sub-interval
    PAC = rep(NA,length(Kvec))
    names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
    weights <- sapply(weights,function(x) sapply(weights,function(y) x+y))
    weights <- weights[lower.tri(weights)]
    for(i in Kvec){
        M = ccm[[i]]
        vecX <- M[lower.tri(M)]
        ord <- order(vecX)
        vecX <- inverse.rle(list(lengths=weights[ord],values=vecX[ord]))
        Fn = stats:::ecdf(vecX)
        PAC[i-1] = Fn(x2) - Fn(x1)
    }
    ## The optimal K
    ## print(PAC)
    optK = Kvec[which.min(PAC)]
    return(optK)
}

## main method
## takes the median of number of clusters as starting point (median K)
## if medianK>2 then derives consensus clustering matrices and optimises PAC
## to find optimal K. Then runs fastConsensusClutering.
## allA is a list of vectors of assignments.
## pMethods should be set to 1, as it can crash the hierachical clustering if random sampling leads to too few unique methods
## repeats should be set to 100 to derive consensus matrices based on 100 random sampling of the methods with replacement
consensusMatrix <- function(allA,
                            pMethods=c(0.5,0.75,0.8,1,1.1,1.2),
                            repeats=2)
{
    nbClustersMethods <-sapply(allA,function(x) length(unique(x[!is.na(x)])))
    medianK <- floor(median(nbClustersMethods))
    if(medianK==1) return(list(clusts=rep(1,length(allA[[1]]))))
    nbClusters <- max(nbClustersMethods,na.rm=T)
    nbClusters <- medianK
    nbMethods <- length(allA)
    nbFeatures <- length(allA[[1]])
    if(medianK>2)
    {
        system.time(allCC <- lapply(pMethods,function(x)
            lapply(1,function(y)
                lapply(1:repeats,function(smp)
                {
                    cat(".")
                    keepMethod <- sample(1:nbMethods,round(nbMethods*x),rep=T)
                    if(nbMethods>4)
                        while(sum(!duplicated(keepMethod))<4)
                            keepMethod <- sample(1:nbMethods,round(nbMethods*x),rep=T)
                    lAa <- lapply(keepMethod,function(a)
                    {
                        allA[[a]]
                    })
                    fastConsensusClustering(lAa,
                                            1:nbClusters,
                                            keepnames=paste("snv",
                                                            1:nbFeatures,
                                                            sep=":"))
                }))))
        clsA <- getClsA(allA)
        l <- length(unique(clsA))
        representants <- which(!duplicated(clsA))
        ccm <- list()
        if(!all(clsA[representants]==unique(clsA))) stop("internal mismatch problem")
        for(i in 2:nbClusters)
        {
            print(i)
            ccm[[i]] <- matrix(0,l,l)
            rownames(ccm[[i]]) <- colnames(ccm[[i]]) <- unique(clsA)
            for(ii in 1:length(pMethods))
            {
                print(paste("ii",ii))
                for(jj in 1)
                {
                    ccm[[i]] <- CCM(ccm[[i]],allCC,ii,jj,repeats,i,representants)
                }
            }
        }
        ccmres <- list(ccm=ccm,clsA=clsA)
        K <- chooseOptimalK(lapply(ccm,function(x)
        {
            x/repeats
        }),
        medianK,
        table(clsA)[rownames(ccm[[2]])])
        clusts <- fastConsensusClustering(allA,
                                          K,
                                          keepnames=paste("snv",
                                                          1:nbFeatures,
                                                          sep=":"))[[1]]
    return(list(K=K,ccmres=ccmres,clusts=clusts,medianK=medianK))
    }
    else
    {
        K <- 2
        clusts <- fastConsensusClustering(allA,
                                          2,
                                          keepnames=paste("snv",
                                                          1:nbFeatures,
                                                          sep=":"))[[1]]
        return(list(K=K,ccmres=NULL,clusts=clusts,medianK=medianK))
    }
    return(list(K=K,ccmres=NULL,clusts=clusts,medianK=medianK))
}

## performs hierarchical clustering and cuts the tree to identify the consensu clusters
findClusts <- function(ccm,nbClust)
{
    hc <- hclust(as.dist(100-ccm),method="ward.D")
    clusts <- cutree(hc,k=nbClust)
}

## general method to identify and remove outliers
findOutliers <- function(lAA,
                         nbOutliers=NULL,
                         downsamplers=NULL,
                         downsampleTo=5000,
                         minMutCoverage=.7)
{
    methods <- names(lAA)
    nbMut <- length(lAA[[1]])
    remNotEnoughMut <- sapply(lAA,function(x) sum(!is.na(x)))
    remNotEnoughMut <- which(!(remNotEnoughMut/nbMut>minMutCoverage | remNotEnoughMut>downsampleTo*minMutCoverage))
    if(nbMut>downsampleTo*2)
        if(!is.null(downsamplers))
        {
            remNotEnoughMut <- unique(c(remNotEnoughMut,which(methods%in%downsamplers)))
        }
    if(length(remNotEnoughMut)>0)
    {
        lAA <- lapply((1:length(lAA))[-c(remNotEnoughMut)],function(x) lAA[[x]])
        names(lAA) <- methods[-c(remNotEnoughMut)]
        return(lAA)
    }
    lAA
}

## transforms lAA after removing a given number of outliers
## minNbMethods is the minimum number of methods: the function wont remove outliers if the input as less than that
transformlAA <- function(lAA,
                         nbOutliers=2,
                         minNbMethods=8)
{
    lAA. <- lapply(1:length(lAA),function(x)
    {
        if(length(lAA[[x]])<minNbMethods) return(lAA[[x]])
        try(findOutliers(lAA[[x]],nbOutliers=nbOutliers),silent=T)
    })
    return(lAA.)
}

## plot histogram of ccfs with clusters coloured by consensus clustering
## clusts is a vector of consensus cluster assignment of the mutations
## ccfs are the cancer cell fractions of the mutations
plotAgreement <- function(clusts,
                          ccfs,
                          ...)
{
    brks <- seq(0,max(ccfs)+0.025,.025)
    hist(ccfs,breaks=brks,col=rgb(.5,.5,.5,.5),...)
    cols <- sample(rainbow(length(unique(clusts))))
    l <- length(unique(clusts))
    for(i in l:1)
    {
        hist(ccfs[clusts%in%unique(clusts)[i:1]],
             breaks=brks,col=cols[i],add=T,...)
    }
    abline(v=sapply(unique(clusts),function(x) mean(ccfs[clusts==x])),lwd=2,col=cols)
}
