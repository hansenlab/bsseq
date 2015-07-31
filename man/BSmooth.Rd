\name{BSmooth}
\alias{BSmooth}
\title{
  BSmooth, smoothing bisulfite sequence data
}
\description{
  This implements the BSsmooth smoothing algorithm for bisulfite
  sequencing data.}
\usage{
BSmooth(BSseq, ns = 70, h = 1000, maxGap = 10^8,
  parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
  mc.cores = 1, keep.se = FALSE, verbose = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{ns}{The minimum number of methylation loci in a smoothing window.}
  \item{h}{The minimum smoothing window, in bases.}
  \item{maxGap}{The maximum gap between two methylation loci, before the
    smoothing is broken across the gap.  The default smoothes each
    chromosome separately.}
  \item{parallelBy}{Should the computation be parallel by chromosome or
    sample, see details.}
  \item{mc.preschedule}{Passed to \code{mclapply} (should the tasks be prescheduled).}
  \item{mc.cores}{Passed to \code{mclapply} (the number of cores used).
  Note that setting \code{mc.cores} to a value greater than 1 is not
  supported on MS Windows, see the help page for \code{mclapply}.}
  \item{keep.se}{Should the estimated standard errors from the smoothing
    algorithm be kept.  This will make the return object roughly 30
    percent bigger and may not be used for anything.}
  \item{verbose}{Should the function be verbose.}
}
\details{
  \code{ns} and \code{h} are passed to the \code{locfit} function.  The
  bandwidth used is the maximum (in genomic distance) of the \code{h}
  and a width big enough to contain \code{ns} number of methylation
  loci.

  The function uses the parallel package to do parallel computations.
  In order to use this, make sure your system have enough RAM, these are
  typically big objects.  The computation can either be split by
  chromosome or by sample, which is better depends on the number of
  samples and how many concurrent smoothings may be done.

}
\value{
  An object of class \code{BSseq}, containing smoothed values and
  optionally standard errors for these.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\references{
  KD Hansen, B Langmead, and RA Irizarry (2012). \emph{BSmooth: from
  whole genome bisulfite sequencing reads to differentially methylated
  regions}. Genome Biology 13:R83.
  \href{http://www.dx.doi.org/10.1186/gb-2012-13-10-r83}{doi:10.1186/gb-2012-13-10-r83}.
}
\seealso{
  \code{\link[locfit]{locfit}} in the locfit package, as well as \code{\linkS4class{BSseq}}.
}
\examples{
\dontrun{
data(BS.chr22)
BS.fit <- BSmooth(BS.chr22, verbose = TRUE)
BS.fit
}
}
