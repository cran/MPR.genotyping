\name{phenoData}
\alias{phenoData}
\docType{data}
\title{
data of phenotype
}
\description{
grain width of RILs
}
\usage{data("phenoData")}
\format{
  The format is:
 num [1:241, 1] 2.85 2.7 2.8 2.9 3.2 2.8 2.75 2.75 2.9 3.35 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:241] "MZ001" "MZ002" "MZ003" "MZ004" ...
  ..$ : chr "Grain_width"
}
\details{
Grain width is quantitative phenotype
}
\source{
http://www.ncpgr.cn/supplements/MPR_genotyping/MPR_genotyping.tar.gz
}
\examples{
#load data
data(phenoData)
}
\keyword{datasets}
