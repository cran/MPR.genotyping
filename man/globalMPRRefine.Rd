\name{globalMPRRefine}
\alias{globalMPRRefine}
\title{
refine SNPs by resampling
}
\description{
The function to refine SNPs by resampling and Bayesian 
inference based on results of globalMPRByMarkers. 
}
\usage{
globalMPRRefine(baseData, markers = NULL, alleleA = NULL, numGroup = ncol(baseData), 
	groupSort = FALSE, numPerm = 10, numTry = 3, numBaseStep = 50, 
	numBaseCandidateStep = numBaseStep * 2, numKnownStep = numBaseStep/2, 
	numKnownCandidateStep = numBaseStep * 2, useMedianToFindKnown = TRUE, 
	maxNStep = 3, scoreMin = 0.8, useOnlyKnownToType = FALSE, useBayes = FALSE, 
	errorRate = 5e-04, saveMidData = FALSE, verbose = FALSE, strSTART = "\r", 
	strEND = "", ...)
}
\arguments{
  \item{baseData}{
Your SNP matrix
}
  \item{markers}{
Your markers data
}
  \item{alleleA}{
We can use the results of globalMPRByMarkers() to be something like markers
}
  \item{numGroup}{
The number of group to run in one permutation.
}
  \item{groupSort}{
Your groups (or RILs) will sort by coverage.(higher coverage will make the result in local region better)
}
  \item{numPerm}{
The number of permutations
}
  \item{numTry}{
The maximum number of the times of using one SNP from one group (or RIL).
}
  \item{numBaseStep}{
The number of SNP to run localMPR().
}
  \item{numBaseCandidateStep}{
The number of SNP candidate in one step
}
  \item{numKnownStep}{
The number of makers to run localMPR().
}
  \item{numKnownCandidateStep}{
The number of makers candidate in one step
}
  \item{useMedianToFindKnown}{
In one local genomic region (window), we will choose the nearest some makers to this region. Median of the region will be the center.
}
  \item{maxNStep}{
The parameter of localMPR()
}
  \item{scoreMin}{
The score is used to control the accuracy of MPR in one local genomic region (window).
}
  \item{useOnlyKnownToType}{
Using different scoring system 
}
  \item{useBayes}{
"useBayes=TRUE" will make this routine integrate with Bayesian inference
}
  \item{errorRate}{
the parameter of genotypeCallsBayes() if you set "useBayes=TRUE".
}
  \item{saveMidData}{
If you want to check the effect of the number of permutations, "saveMidData = TRUE" will help you.
}
  \item{verbose}{
Report verbose progress.
}
  \item{strSTART}{
The part of displaying format (verbose)
}
  \item{strEND}{
The part of displaying format (verbose)
}
  \item{\dots}{
Arguments to be passed to other methods.
}
}
\seealso{
\code{\link{globalMPRByMarkers}}
}
\examples{
\donttest{
data(snpData)
data(markerData)

## select 30 markers randomly
set.seed(123);markers <- sample(names(markerData)[10:50],20)

## select SNP sites which contain the 30 markers
ids <- match(markers,rownames(snpData))
str(myBaseData <- snpData[min(ids):max(ids),])

## global MPR aiding with marker data
allele.MPR <- globalMPRByMarkers(myBaseData,markers=markerData,numTry=3,
			numBaseStep=50,numBaseCandidateStep=100,
            numMarkerStep=10,useMedianToFindKnown=TRUE,
			maxNStep=3,scoreMin=0.8,verbose=TRUE)
			 
## then you need to refine the MPR results
set.seed(123);system.time(all.res <- globalMPRRefine(myBaseData,alleleA=na.omit(
            allele.MPR[,1]),numGroup=238,groupSort=TRUE,numPerm=1,numTry=3,
            numBaseStep=50,numBaseCandidateStep=100,numKnownStep=30,
            numKnownCandidateStep=50,useMedianToFindKnown=TRUE,
            maxNStep=3,scoreMin=0.8,saveMidData=TRUE,verbose=TRUE))
}
}
\keyword{ MPR }
