\name{snpData}
\alias{snpData}
\docType{data}
\title{
Example SNP data
}
\description{
A total of 238 RILs developed from the cross be-
tween Zhenshan 97 and Minghui 63 were resequenced with an
Illumina Genome Analyzer II using the bar-coded multiplexed
sequencing approach
}
\usage{data("snpData")}
\format{
  The format is:
 chr [1:15795, 1:238] "A" NA NA NA NA NA NA "G" NA NA NA NA NA ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:15795] "0500000526A" "0500000556A" "0500000559G" "0500000591G" ...
  ..$ : chr [1:238] "MZ001" "MZ002" "MZ003" "MZ004" ...
}
\details{
SNP data format as it
}
\source{
http://www.ncpgr.cn/supplements/MPR_genotyping/MPR_genotyping.tar.gz
}
\examples{
\donttest{
# load data
data(snpData)
}
}
\keyword{datasets}
