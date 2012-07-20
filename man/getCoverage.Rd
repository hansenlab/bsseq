\name{getCoverage}
\alias{getCoverage}
\title{
  Obtain coverage for BSseq objects.
}
\description{
  Obtain coverage for BSseq objects.
}
\usage{
getCoverage(BSseq, regions = NULL, type = c("Cov", "M"),
  what = c("perBase", "perRegionAverage", "perRegionTotal"))
}
\arguments{
  \item{BSseq}{An object of class \code{"BSseq"}.}
  \item{regions}{An optional \code{"data.frame"} or 
    \code{"GenomicRanges"} object specifying a number of genomic regions.}
  \item{type}{This returns either coverage or the total
    evidence for methylation at the loci.}
  \item{what}{The type of return object, see details.}
}
\value{
  If \code{regions} are not specified (\code{regions = NULL}) a matrix
  (\code{what = "perBase"}) or a vector (otherwise) is returned.  This
  will either contain the per-base coverage or the genome total or
  average coverage.

  If \code{what = "perBase"} and \code{regions} are specified, a list is
  returned.  Each element of the list is a matrix corresponding to the
  genomic loci inside the region.  It is conceptually the same as
  splitting the coverage by region.

  If \code{what = "perRegionAverage"} or \code{what = "perRegionTotal"}
  and \code{regions} are specified the return value is a matrix.  Each
  row of the matrix corresponds to a region and contains either the
  total coverage of the average coverage in the region.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{BSseq}} for the \code{"BSseq"} class.
}
\examples{
  head(getCoverage(BS.chr22, type = "M"))
  reg <- GRanges(seqnames = c("chr22", "chr22"),
    ranges = IRanges(start = c(1, 2*10^7), end = c(2*10^7 +1, 4*10^7)))
  getCoverage(BS.chr22, regions = reg, what = "perRegionAverage")
  cList <- getCoverage(BS.chr22, regions = reg)
  length(cList)
  head(cList[[1]])
}
