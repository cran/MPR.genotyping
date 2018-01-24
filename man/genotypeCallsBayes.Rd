\name{genotypeCallsBayes}
\alias{genotypeCallsBayes}
\title{
genotype Calls by Bayes
}
\description{
After globalMPRRefine, you should use Bayes to call final genotype. 
}
\usage{
genotypeCallsBayes(ALLELE.num, errorRate = 5e-04, eps = 1e-10, 
		maxIterate = 100, verbose = FALSE)
}
\arguments{
  \item{ALLELE.num}{
one result of globalMPRRefine
}
  \item{errorRate}{
errorRate
}
  \item{eps}{
eps
}
  \item{maxIterate}{
maxIterate to control running time
}
  \item{verbose}{
choose to show details
}
}
\seealso{
\code{\link{globalMPRRefine}}
}
\examples{
ALLELE.num=matrix(c(100,1,
					50,55,
					89,2,
					1,101,
					0,78),5,3,byrow=TRUE)
genotypeCallsBayes(ALLELE.num)
}
\keyword{ Bayes }
