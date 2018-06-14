\name{getMeth}
\alias{getMeth}
\title{
  Obtain methylation estimates for BSseq objects.
}
\description{
  Obtain methylation estimates for BSseq objects, both smoothed and raw.
}
\usage{
getMeth(BSseq, regions = NULL, type = c("smooth", "raw"),
  what = c("perBase", "perRegion"), confint = FALSE, alpha = 0.95,
  withDimnames = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{regions}{An optional \code{data.frame} or
    \code{GenomicRanges} object specifying a number of genomic regions.}
  \item{type}{This returns either smoothed or raw estimates of the methylation level.}
  \item{what}{The type of return object, see details.}
  \item{confint}{Should a confidence interval be return for the
    methylation estimates (see below).  This is only supported if
    \code{what} is equal to \code{perBase}.}
  \item{alpha}{alpha value for the confidence interval.}
  \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames should
    be applied to extracted coverage elements. Setting
    \code{withDimnames = FALSE} increases the speed and memory efficiency with
    which coverage is extracted.}
}
\note{
  A \code{BSseq} object needs to be smoothed by the function
  \code{BSmooth} in order to support \code{type = "smooth"}.
}
\value{
  \strong{NOTE:} The return type of \code{getMeth} varies depending on its
  arguments.

  If \code{region = NULL} the \code{what} argument is ignored.  This is
  also the only situation in which \code{confint = TRUE} is supported.
  The return value is either a \linkS4class{DelayedMatrix}
  (\code{confint = FALSE} or a list with three \linkS4class{DelayedMatrix}
  components \code{confint = TRUE} (\code{meth}, \code{upper} and
  \code{lower}), giving the methylation estimates and (optionally) confidence
  intervals.

  Confidence intervals for \code{type = "smooth"} is based on standard
  errors from the smoothing algorithm (if present).  Otherwise it is
  based on pointwise confidence intervals for binomial distributions
  described in Agresti (see below), specifically the score confidence
  interval.

  If \code{regions} are specified, \code{what = "perBase"} will make the
  function return a list, each element of the list being a
  \linkS4class{DelayedMatrix} corresponding to a genomic region (and each row
  of the \linkS4class{DelayedMatrix} being a loci inside the region). If
  \code{what = "perRegion"} the function returns a \linkS4class{DelayedMatrix},
  with each row corresponding to a region and containing the average
  methylation level in that region.
}
\references{
  A Agresti and B Coull.
  \emph{Approximate Is Better than "Exact" for Interval Estimation of
    Binomial Proportions}.
  The American Statistician (1998) 52:119-126.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{BSseq}} for the \code{BSseq} class and
  \code{\link{BSmooth}} for smoothing such an object.
}
\examples{
data(BS.chr22)
head(getMeth(BS.chr22, type = "raw"))
reg <- GRanges(seqnames = c("chr22", "chr22"),
  ranges = IRanges(start = c(1, 2*10^7), end = c(2*10^7 +1, 4*10^7)))
head(getMeth(BS.chr22, regions = reg, type = "raw", what = "perBase"))

#-------------------------------------------------------------------------------
# An example using a HDF5Array-backed BSseq object
#

library(HDF5Array)
# See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
hdf5_BS.chr22 <- saveHDF5SummarizedExperiment(x = BS.chr22,
                                              dir = tempfile())
head(getMeth(hdf5_BS.chr22, type = "raw"))
head(getMeth(hdf5_BS.chr22, regions = reg, type = "raw", what = "perBase"))
}
