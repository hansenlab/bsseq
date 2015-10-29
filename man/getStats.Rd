\name{getStats}
\alias{getStats}
\title{
  Obtain statistics from a BSseqTstat object
}
\description{
  Essentially an accessor function for the statistics of a \code{BSseqTstat}
  object.
}
\usage{
getStats(bstat, regions = NULL, ...)
}
\arguments{
  \item{bstat}{An object of class \code{BSseqStat} or \code{BSseqTstat}.}
  \item{regions}{An optional \code{data.frame} or
    \code{GenomicRanges} object specifying a number of genomic
    regions.}
  \item{...}{Additional arguments passed to the different backends based
    on the class of \code{bstat}; see Details.}
}
\value{
  An object of class \code{data.frame} possible restricted to the
  regions specified.
}
\details{
  Additional argument when the \code{bstat} object is of class \code{BSseqTstat}:
  \describe{
    \item{stat}{Which statistics column should be obtained.}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  \code{\linkS4class{BSseqTstat}} for the \code{BSseqTstat} class, and
  \code{\link{getCoverage}} and \code{\link{getMeth}} for similar
  functions, operating on objects of class \code{BSseq}.
}
\examples{
if(require(bsseqData)) {
  data(BS.cancer.ex.tstat)
  head(getStats(BS.cancer.ex.tstat))
  reg <- GRanges(seqnames = c("chr22", "chr22"),
     ranges = IRanges(start = c(1, 2*10^7), end = c(2*10^7 +1, 4*10^7)))
  head(getStats(BS.cancer.ex.tstat, regions = reg))
}
}
