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
  what = c("perBase", "perRegionAverage", "perRegionTotal"),
  withDimnames = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{regions}{An optional \code{data.frame} or
    \code{GenomicRanges} object specifying a number of genomic regions.}
  \item{type}{This returns either coverage or the total
    evidence for methylation at the loci.}
  \item{what}{The type of return object, see details.}
  \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames should
    be applied to extracted coverage elements. Setting
    \code{withDimnames = FALSE} increases the speed and memory efficiency with
    which coverage is extracted.}
}
\value{
  \strong{NOTE:} The return type of \code{getCoverage} varies depending on its
  arguments.

  If \code{regions} are not specified (\code{regions = NULL}) a
  \linkS4class{DelayedMatrix} object (\code{what = "perBase"}) is returned.
  This will either contain the per-base coverage, the average coverage, or the
  genome total coverage (depending on value of \code{what}).

  If \code{what = "perBase"} and \code{regions} are specified, a list is
  returned.  Each element of the list is a \linkS4class{DelayedMatrix} object
  corresponding to the genomic loci inside the region.  It is conceptually the
  same as splitting the coverage by region.

  If \code{what = "perRegionAverage"} or \code{what = "perRegionTotal"}
  and \code{regions} are specified the return value is a
  \linkS4class{DelayedMatrix} object. Each row of the
  \linkS4class{DelayedMatrix} corresponds to a region and contains either the
  average coverage or the total coverage in the region.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{BSseq}} for the \code{BSseq} class.
}
\examples{
data(BS.chr22)
head(getCoverage(BS.chr22, type = "M"))
reg <- GRanges(seqnames = c("chr22", "chr22"),
  ranges = IRanges(start = c(1, 2*10^7), end = c(2*10^7 +1, 4*10^7)))
getCoverage(BS.chr22, regions = reg, what = "perRegionAverage")
  cList <- getCoverage(BS.chr22, regions = reg)
length(cList)
head(cList[[1]])

#-------------------------------------------------------------------------------
# An example using a HDF5Array-backed BSseq object
#

library(HDF5Array)
# See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
hdf5_BS.chr22 <- saveHDF5SummarizedExperiment(x = BS.chr22,
                                              dir = tempfile())
head(getCoverage(hdf5_BS.chr22, type = "M"))
reg <- GRanges(seqnames = c("chr22", "chr22"),
               ranges = IRanges(start = c(1, 2 * 10 ^ 7),
               end = c(2 * 10 ^ 7 + 1, 4 * 10 ^ 7)))
getCoverage(hdf5_BS.chr22, regions = reg, what = "perRegionAverage")
hdf5_cList <- getCoverage(hdf5_BS.chr22, regions = reg)
length(hdf5_cList)
head(hdf5_cList[[1]])
}
