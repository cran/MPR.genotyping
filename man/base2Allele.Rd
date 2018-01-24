\name{base2Allele}
\alias{base2Allele}
\title{
Get Initial Allele from SNP Matrix
}
\description{
Get initial parental genotypes from SNP matrix randomly.
}
\usage{
base2Allele(baseData = NULL)
}
\arguments{
  \item{baseData}{
matrix of SNP
}
}
\details{
base : matrix of SNP
}
\value{
matrix :2 columns and nrows of SNP data
}
\seealso{
\code{\link{globalMPRByMarkers}}
}
\examples{
base=matrix(c("A",NA,"G",
		NA,"T","A",
		"C","G",NA),3,3,byrow=TRUE)
base2Allele(base)
}
\keyword{ formatting }
