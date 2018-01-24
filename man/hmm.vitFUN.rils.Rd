\name{hmm.vitFUN.rils}
\alias{hmm.vitFUN.rils}
\title{
Genotyping by HMM
}
\description{
correct genotype by Hidden Markov Model
}
\usage{
hmm.vitFUN.rils(geno, position, geno.probability, transitionFUN = 
	phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), ...)
}
\arguments{
  \item{geno}{
a numeric matrix of SNP genotype formatted
}
  \item{position}{
vector of SNP positions
}
  \item{geno.probability}{
vector of start probability
}
  \item{transitionFUN}{
a function about transition
}
  \item{emissionFUN}{
a function about emission
}
  \item{\dots}{
arguments to be passed to other methods.
}
}
\examples{
O <- c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,
     2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1,2,2,1,2,1,2,1,2,1,2,1,2,1,1,1,1,1,1,1,2,1,1,1,1,
     1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1)
O.pos <- 30e3*1:length(O)

O.cr <- hmm.vitFUN.rils(geno=O,position=O.pos,geno.probability=c(0.4975, 0.4975,0.005),
				transitionFUN = phy2get.haldane.rils,
				emissionFUN = makeEmissionFUN(errorRate = 0.01))

rbind(raw=paste(O,collapse=''),correct=paste(O.cr,collapse=''))
##         [,1]
## raw     "11111112111111211111112222222212222222212222212222222222
##				22221211221212121212111111121111111111121111111211111"
## correct "11111111111111111111112222222222222222222222222222222222
##				22223333333333333333111111111111111111111111111111111"
}
\keyword{ HMM }
