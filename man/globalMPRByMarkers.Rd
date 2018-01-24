\name{globalMPRByMarkers}
\alias{globalMPRByMarkers}
\title{
MPR inference in whole chromosome
}
\description{
The function to do MPR inference in whole chromosome 
by using localMPR to infer parental genotypes in hundreds of local regions 
and assemble them aiding with low-coverage sequences of one parent or known 
markers.
}
\usage{
globalMPRByMarkers(baseData, markers = NULL, alleleA = NULL, numTry = 3,
	numBaseStep = 50, numBaseCandidateStep = numBaseStep * 2,
	numKnownStep = pmax(numBaseStep/5, 10),
	numKnownCandidateStep = numKnownStep * 1.5,
	useMedianToFindKnown = TRUE, maxNStep = 3, scoreMin = 0.8, verbose = FALSE,
	strSTART = "\r", strEND = "", ...)
}
\arguments{
  \item{baseData}{
(Necessary input) character matrix of SNP dataset
}
  \item{markers}{
character vector of markers data with SNP position names
}
  \item{alleleA}{
(Necessary input) character vector of one parent allele
}
  \item{numTry}{
maximum number of the times of using one SNP from one group (or RIL).
}
  \item{numBaseStep}{
number of SNP to run localMPR().
}
  \item{numBaseCandidateStep}{
number of SNP candidate in one step
}
  \item{numKnownStep}{
number of makers to run localMPR().
}
  \item{numKnownCandidateStep}{
number of makers candidate in one step
}
  \item{useMedianToFindKnown}{
In one local genomic region (window), we will choose the nearest some makers 
		to this region. Median of the region will be the center.
}
  \item{maxNStep}{
parameter of localMPR()
}
  \item{scoreMin}{
be used to control the accuracy of MPR in one local genomic region (one window).
}
  \item{verbose}{
report verbose progress
}
  \item{strSTART}{
part of displaying format (verbose)
}
  \item{strEND}{
part of displaying format (verbose)
}
  \item{\dots}{
arguments to be passed to other methods.
}
}
\seealso{
\code{\link{localMPR}}
}
\examples{
\donttest{
## load sample dataset
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
}
}
\keyword{ MPR }