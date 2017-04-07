\name{GoodnessOfFit}
\alias{poissonGoodnessOfFit}
\alias{poissonGoodnessOfFit}
\alias{chisqGoodnessOfFit}
\alias{binomialGoodnessOfFit}
\alias{print.chisqGoodnessOfFit}
\alias{plot.chisqGoodnessOfFit}
\title{
  Binomial and poisson goodness of fit statistics for BSSeq objects
}
\description{
  Binomial and poisson goodness of fit statistics for BSSeq objects,
  including plotting capability.
}
\usage{
poissonGoodnessOfFit(BSseq, nQuantiles = 10^5)
binomialGoodnessOfFit(BSseq, method = c("MLE"), nQuantiles = 10^5)
\S3method{print}{chisqGoodnessOfFit}(x, \ldots)
\S3method{plot}{chisqGoodnessOfFit}(x, type = c("chisq", "pvalue"), plotCol = TRUE, qqline = TRUE,
  pch = 16, cex = 0.75, \ldots)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{x}{A chisqGoodnessOfFit object (as produced by
    \code{poissonGoodnessOfFit} or \code{binomialGoodnessOfFit}).}
  \item{nQuantiles}{The number of (evenly-spaced) quantiles stored in
    the return object.}
  \item{method}{How is the parameter estimated.}
  \item{type}{Are the chisq or the p-values being plotted.}
  \item{plotCol}{Should the extreme quantiles be colored.}
  \item{qqline}{Add a \code{qqline}.}
  \item{pch, cex}{Plotting symbols and size.}
  \item{...}{Additional arguments being passed to \code{qqplot} (for
    \code{plot}) or ignored (for \code{print}).}
}
\details{

  These functions compute and plot goodness of fit statistics for
  \code{BSseq} objects.  For each methylation loci, the Poisson
  goodness of fit statistic tests whether the coverage (at that loci) is
  independent and identically Poisson distributed across the samples.
  In a similar fashion, the Binomial goodness of fit statistic tests
  whether the number of reads supporting methylation are independent and
  identically binomial distributed across samples (with different size
  parameters given by the coverage vector).

  These functions do not handle \code{NA} values.
}
\value{
  The plotting method is invoked for its side effect.  Both
  \code{poissonGoodnessOfFit} and \code{binomialGoodnessOfFit} returns
  an object of class \code{chisqGoodnessOfFit} which is a list with components
 \item{chisq}{a vector of Chisq values.}
 \item{quantiles}{a vector of quantiles (of the chisq values).}
 \item{df}{degress of freedom}
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  For the plotting method, see \code{qqplot}.
}
\examples{
if(require(bsseqData)) {
  data(BS.cancer.ex)
  BS.cancer.ex <- updateObject(BS.cancer.ex)
  gof <- poissonGoodnessOfFit(BS.cancer.ex)
  plot(gof)

#-------------------------------------------------------------------------------
# An example using a HDF5Array-backed BSseq object
#
\donttest{
    library(HDF5Array)

    # See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
    hdf5_BS.cancer.ex <- saveHDF5SummarizedExperiment(x = BS.cancer.ex,
                                                      dir = tempfile())
    hdf5_gof <- poissonGoodnessOfFit(hdf5_BS.cancer.ex)
    plot(hdf5_gof)
}
}
}
