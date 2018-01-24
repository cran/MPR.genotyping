\name{NumRecomEvents}
\alias{NumRecomEvents}
\title{
calculate the number of recombination events
}
\description{
As the parental genotypes we assumed, this function can calculate how many recombination events would be needed to produce genotypes of RILs.
}
\usage{
NumRecomEvents(baseData, allele.matrix, genoData = NULL)
}
\arguments{
  \item{baseData}{
matirx of SNP
}
  \item{allele.matrix}{
matrix: alleles of parental genotypes
}
  \item{genoData}{
matirx of SNP which has been formatted. 
}
}
\details{
We wrote core code for calculating the number of recombination events in C to improve computational speed.
}
\value{
\item{R}{the number of recombination events}
}
\seealso{
\code{\link{localMPR}}
}
\examples{
set.seed(123)
data(myBaseData)
allele.random <- base2Allele(myBaseData)
NumRecomEvents(myBaseData,allele.random)
# 162
}
\keyword{ MPR }
