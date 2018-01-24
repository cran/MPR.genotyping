\name{base2Geno_ori}
\alias{base2Geno_ori}
\title{
Formatting the Matrix of SNP to Int Type(NA)
}
\description{
formatting the matrix of SNP to int type(NA)
}
\usage{
base2Geno_ori(baseData = NULL, allele.matrix = NULL)
}
\arguments{
  \item{baseData}{
the matrix of snp
}
  \item{allele.matrix}{
the allele of each snp site
}
}
\details{
the snp from parent1 will become 0,
the snp is not detected by sequencing will become NA,
the snp from parent2 will become 1.
}
\value{
returned a matrix which have same size of matrix of snp, but "A","G","C","T" become "0","NA","1". 
}
\seealso{
\code{\link{base2Geno}}
}
\examples{
base=matrix(c("A",NA,"G",
		NA,"T","A",
		"C","G",NA),3,3,byrow=TRUE)
allele=matrix(c("A","G",
			"T","A",
			"C","G"),3,2,byrow=TRUE)
base2Geno_ori(base,allele)
}
\keyword{ formatting }
