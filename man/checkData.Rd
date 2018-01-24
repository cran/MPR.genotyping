\name{checkData}
\alias{checkData}
\docType{data}
\title{
Data for check
}
\description{
this is used to check up the genotype results in my example.
}
\usage{data("checkData")}
\format{
  The format is:
 chr [1:11948, 1:2] "A" "A" "C" "T" "A" ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:11948] "0500000526A" "0500000556A" "0500000559G" "0500000591G" ...
  ..$ : chr [1:2] "ZS97" "MH63"
}
\details{
you can use "table(checkData[ids,1]==alleleA)"
}
\source{
http://www.ncpgr.cn/supplements/MPR_genotyping/MPR_genotyping.tar.gz
}
\examples{
#load data
data(checkData)
}
\keyword{datasets}
