\name{getCoverage_mod}
\alias{getCoverage_mod}
\title{
  Obtain coverage for modified BSseq objects.
}
\description{
  Obtain coverage for modified BSseq objects.
}
\usage{
getCoverage_mod(BSseq_mod, regions = NULL, type = c("Cov", "M"),
  what = c("perBase", "perRegionAverage", "perRegionTotal"),
  withDimnames = TRUE)
}
\arguments{
  \item{BSseq_mod}{An object of class \code{BSseq_mod}.}
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
  \strong{NOTE:} The return type of \code{getCoverage_mod} varies depending on its
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
  \code{\linkS4class{BSseq_mod}} for the \code{BSseq_mod} class.
}
\examples{
M <- U <- Other_mod <- matrix(1:9, 3,3)
colnames(M) <- c("A1", "A2", "A3")
BStest <- BSseq_mod(pos = 1:3, chr = c("chr1", "chr2", "chr1"), M = M, U = U,
                    Other_mod = Other_mod)
head(getCoverage_mod(BStest, type = "M"))
head(getCoverage_mod(BStest, type = "Cov", what = "perRegionAverage"))
}
