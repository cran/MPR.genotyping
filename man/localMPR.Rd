\name{localMPR}
\alias{localMPR}
\title{
infer parental genotypes by minimizing the number of recombination events
}
\description{
The core function to infer parental genotypes (MPR inference) in 
a local chromosome region by minimizing the number of recombination events in 
the population. 
}
\usage{
localMPR(baseData, allele.matrix = NULL, maxNStep = 5,
	returnNumRecomEvents = FALSE, verbose = FALSE, strEND = "\n", ...)
}
\arguments{
  \item{baseData}{
a character matrix of SNP
}
  \item{allele.matrix}{
a character matrix of alleles
}
  \item{returnNumRecomEvents}{
whether report the number of recombination events in the final result
}
  \item{maxNStep}{
the number of SNP in an exchange of the alleles.
}
  \item{verbose}{
report verbose progress.
}
  \item{strEND}{
the part of displaying format (verbose)
}
  \item{\dots}{
arguments to be passed to other methods.
}
}
\details{
Several factors may affect the accuracy of MPR inference: 
(a). The number of SNPs processed each time (window size);
(b). The density of putative SNPs (the distance between SNP sites);
(c). The maximum step size of the heuristic perturbation (the parameter of 
maxNStep);
(d). The number of RILs.
And we wrote core code in C to improve computational speed.
}
\value{
the result of local MPR genotyping allele
}
\seealso{
\code{\link{NumRecomEvents},\link{globalMPRByMarkers},\link{globalMPRRefine},}
}
\examples{
set.seed(123)
## select 50 SNP sites to test inference of parental genotypes from snpData.rda
data(myBaseData)

## random assignments of parental genotypes to alleles will result in
## a big number of recombinations 
allele.random <- base2Allele(myBaseData)

## a big number of recombinations
NumRecomEvents(myBaseData,allele.random)
## 162

## MPR inference with maximum step size of 5
allele.MPR <- localMPR(baseData=myBaseData,maxIterate=50,maxNStep=5,showDetail=TRUE)

## should be a small number compared with random assignments above 
NumRecomEvents(myBaseData,allele.MPR)
## 33
}
\keyword{ MPR }