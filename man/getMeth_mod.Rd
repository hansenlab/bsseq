\name{getMeth_mod}
\alias{getMeth_mod}
\title{
  Obtain methylation estimates for modified BSseq objects.
}
\description{
  Obtain methylation estimates for modified BSseq objects, both smoothed and raw.
}
\usage{
getMeth_mod(BSseq_mod, regions = NULL, type = c("smooth", "raw"),
  what = c("perBase", "perRegion"), confint = FALSE, alpha = 0.95,
  withDimnames = TRUE)
}
\arguments{
  \item{BSseq_mod}{An object of class \code{BSseq_mod}.}
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
  A \code{BSseq_mod} object needs to be smoothed by the function
  \code{BSmooth} in order to support \code{type = "smooth"}.
}
\value{
  \strong{NOTE:} The return type of \code{getMeth_mod} varies depending on its
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
  \code{\linkS4class{BSseq_mod}} for the \code{BSseq_mod} class and
  \code{\link{BSmooth}} for smoothing such an object.
}
\examples{
M <- matrix(0:8, 3, 3)
U <- matrix(1:9, 3, 3)
Other_mod <- matrix(1:9, 3, 3)
BS <- BSseq_mod(chr = c("chr1", "chr2", "chr1"), pos = c(1,2,3),
        M = M, U = U, Other_mod = Other_mod, sampleNames = c("A","B", "C"))
head(getMeth_mod(BS, type = "raw"))
}
