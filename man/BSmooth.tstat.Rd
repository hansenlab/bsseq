\name{BSmooth.tstat}
\alias{BSmooth.tstat}
\title{
  Compute t-statistics based on smoothed whole-genome bisulfite
  sequencing data.
}
\description{
  Compute t-statistics based on smoothed whole-genome bisulfite
  sequencing data.
}
\usage{
BSmooth.tstat(BSseq, group1, group2,
  estimate.var = c("same", "paired", "group2"), local.correct = TRUE,
  maxGap = NULL, qSd = 0.75, k = 101, mc.cores = 1, verbose = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{group1}{A vector of sample names or indexes for the \sQuote{treatment}
    group.}
  \item{group2}{A vector of sample names or indexes for the \sQuote{control}
    group.}
  \item{estimate.var}{How is the variance estimated, see details.}
  \item{local.correct}{A logical; should local correction be used, see
    details.}
  \item{maxGap}{A scalar greater than 0, see details.}
  \item{qSd}{A scalar between 0 and 1, see details.}
  \item{k}{A positive scalar, see details.}
  \item{mc.cores}{The number of cores used.  Note that setting
    \code{mc.cores} to a value greater than 1 is not supported on MS
    Windows, see the help page for \code{mclapply}.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  T-statistics are formed as the difference in means between group 1 and
  group 2 divided by an estimate of the standard deviation, assuming
  that the variance in the two groups are the same (\code{same}), that
  we have paired samples (\code{paired}) or only estimate the variance
  based on group 2 (\code{group2}).  The standard deviation estimates
  are then smoothed (using a running mean with a width of \code{k}) and
  thresholded (using \code{qSd} which sets the minimum standard
  deviation to be the \code{qSd}-quantile).  Optionally, the
  t-statistics are corrected for low-frequency patterns.

  It is sometimes useful to use \code{local.correct} even if no large
  scale changes in methylation have been found; it makes the marginal
  distribution of the t-statistics more symmetric.

  Additional details in the reference.
}
\value{
  An object of class \code{BSseqTstat}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\references{
 KD Hansen, B Langmead, and RA Irizarry.
  \emph{BSmooth: from whole genome bisulfite sequencing reads to
    differentially methylated regions}.
  Genome Biology (2012) 13:R83.
  doi:\href{http://www.dx.doi.org/10.1186/gb-2012-13-10-r83}{10.1186/gb-2012-13-10-r83}.
}
\seealso{
  \code{\link{BSmooth}} for the input object and
  \code{\linkS4class{BSseq}} for its class.
  \code{\linkS4class{BSseqTstat}} describes the return class.  This
  function is likely to be followed by the use of
  \code{\link{dmrFinder}}.  And finally, see the package vignette(s) for
  more information on how to use it.
}
\examples{
if(require(bsseqData)) {
  data(keepLoci.ex)
  data(BS.cancer.ex.fit)
  BS.cancer.ex.fit <- updateObject(BS.cancer.ex.fit)
  ## Remember to subset the BSseq object, see vignette for explanation
  BS.tstat <- BSmooth.tstat(BS.cancer.ex.fit[keepLoci.ex,],
                            group1 = c("C1", "C2", "C3"),
                            group2 = c("N1", "N2", "N3"),
                            estimate.var = "group2")
  BS.tstat
  ## This object is also stored as BS.cancer.ex.tstat in the
  ## bsseqData package

  #---------------------------------------------------------------------------
  # An example using a HDF5Array-backed BSseq object
  #
  \donttest{
    library(HDF5Array)
    # See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
    hdf5_BS.cancer.ex.fit <- saveHDF5SummarizedExperiment(
        x = BS.cancer.ex.fit[keepLoci.ex, ],
        dir = tempfile())
    hdf5_BS.tstat <- BSmooth.tstat(hdf5_BS.cancer.ex.fit,
                                   group1 = c("C1", "C2", "C3"),
                                   group2 = c("N1", "N2", "N3"),
                                   estimate.var = "group2")
  }
}
}
