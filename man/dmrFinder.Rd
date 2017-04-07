\name{dmrFinder}
\alias{dmrFinder}
\title{
  Finds differentially methylated regions for whole genome bisulfite
  sequencing data.
}
\description{
  Finds differentially methylated regions for whole genome bisulfite
  sequencing data.  Essentially identifies regions of the genome where
  all methylation loci have an associated t-statistic that is beyond a
  (low, high) cutoff.
}
\usage{
dmrFinder(bstat, cutoff = NULL, qcutoff = c(0.025, 0.975),
  maxGap=300, stat = "tstat.corrected", verbose = TRUE)
}
\arguments{
  \item{bstat}{An object of class \code{BSseqStat} or \code{BSseqTstat}.}
  \item{cutoff}{The cutoff of the t-statistics.  This should be a vector
    of length two giving the (low, high) cutoff.  If \code{NULL}, see
    \code{qcutoff}.}
  \item{qcutoff}{In case \code{cutoff} is \code{NULL}, compute the
    cutoff using these quantiles of the t-statistic.}
  \item{maxGap}{If two methylation loci are separated by this distance,
    break a possible DMR.  This guarantees that the return DMRs have CpGs
    that are this distance from each other.}
  \item{stat}{Which statistic should be used?}
  \item{verbose}{Should the function be verbose?}
}
\details{
  The workhorse function is \code{BSmooth.tstat} which sets up a
  t-statistic for a comparison between two groups.

  Note that post-processing of these DMRs are likely to be necessary,
  filtering for example for length (or number of loci).
}
\value{
  A \code{data.frame} with columns
    \item{start,end,width,chr}{genomic locations and width.}
    \item{n}{The number of methylation loci.}
    \item{invdensity}{Average length per loci.}
    \item{group1.mean}{The mean methylation level across samples and
      loci in 'group1'.}
    \item{group2.mean}{The mean methylation level across samples and
      loci in 'group2'.}
    \item{meanDiff}{The mean difference in methylation level; the
      difference between \code{group1.mean} and \code{group2.mean}.}
    \item{idxStart, idxEnd, cluster}{Internal use.}
    \item{areaStat}{The 'area' of the t-statistic; equal to the sum of
      the t-statistics for the individual methylation loci.}
    \item{direction}{either \sQuote{hyper} or \sQuote{hypo}.}
    \item{areaStat.corrected}{Only present if \code{column =
	"tstat.corrected"}, contains the area of the corrected t-statistics.}
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\references{
  KD Hansen, B Langmead, and RA Irizarry.
  \emph{BSmooth: from whole genome bisulfite sequencing reads to
    differentially methylated regions}.
  Genome Biology (2012) 13:R83.
  doi:\href{http://www.dx.doi.org/10.1186/gb-2012-13-10-r83}{10.1186/gb-2012-13-10-r83}.
}
\seealso{
  \code{\link{BSmooth.tstat}} for the function constructing the input
    object, and \code{\linkS4class{BSseqTstat}} for its class.  In the
    example below, we use \code{\link[bsseqData]{BS.cancer.ex.tstat}} as
    the actual input object.  Also see the package vignette(s) for a
    detailed example.
}
\examples{
if(require(bsseqData)) {
  dmrs0 <- dmrFinder(BS.cancer.ex.tstat, cutoff = c(-4.6, 4.6), verbose = TRUE)
  dmrs <- subset(dmrs0, abs(meanDiff) > 0.1 & n >= 3)
}
}
