\name{myBaseData}
\alias{myBaseData}
\docType{data}
\title{
Example SNP data(smaller)
}
\description{
# set.seed(123);myBaseData <- snpData[sample(200,50),]
}
\usage{data("myBaseData")}
\format{
  The format is:
 chr [1:50, 1:238] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA "G" NA NA NA NA NA NA NA NA NA NA NA NA NA ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:50] "0500069986C" "0500170322G" "0500097599A" "0500174545T" ...
  ..$ : chr [1:238] "MZ001" "MZ002" "MZ003" "MZ004" ...
}
\details{
smaller
}
\source{
snpData
}
\examples{
## load data
data(myBaseData)
}
\keyword{datasets}
