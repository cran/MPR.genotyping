\name{base2Geno}
\alias{base2Geno}
\title{
Formatting the Matrix of SNP to Int Type
}
\description{
formatting the matrix of snp to int type
}
\usage{
base2Geno(baseData = NULL, allele.matrix = NULL)
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
the snp from parent1 will become 1,
the snp is not detected by sequencing will become 0,
the snp from parent2 will become -1.
}
\value{
returned a matrix which have same size of matrix of snp, but "A","G","C","T" become "1","0","-1". 
}
\seealso{
\code{\link{globalMPRByMarkers}}
}
\examples{
base=matrix(c("A",NA,"G",
		NA,"T","A",
		"C","G",NA),3,3,byrow=TRUE)
allele=matrix(c("A","G",
			"T","A",
			"C","G"),3,2,byrow=TRUE)
base2Geno(base,allele)
}
\keyword{ formatting }
