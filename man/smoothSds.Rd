\name{smoothSds}
\alias{smoothSds}
\title{
    Smooth the standard deviations using a thresholded running mean based on
    smoothed whole-genome bisulfite sequencing data.
}
\description{
    Smooth the standard deviations using a thresholded running mean based on
    smoothed whole-genome bisulfite sequencing data.
}
\usage{
    smoothSds(BSseqStat, k = 101, qSd = 0.75, mc.cores = 1, maxGap = 10^8,
              verbose = TRUE)
}

\arguments{
    \item{BSseqStat}{An object of class \code{BSseqStat}, typically an object
        returned by \code{\link{BSmooth.fstat}(...)} and not constructed by
        the user.}
    \item{k}{A positive scalar, see details.}
    \item{qSd}{A scalar between 0 and 1, see details.}
    \item{mc.cores}{The number of cores used.  Note that setting
        \code{mc.cores} to a value greater than 1 is not supported on MS
        Windows, see the help page for \code{mclapply}.}
    \item{maxGap}{A scalar greater than 0, see details.}
    \item{verbose}{Should the function be verbose?}
}
\details{
    The standard deviation estimates are smoothed using a running mean with a
    width of \code{k} and thresholded using \code{qSd} which sets the minimum
    standard deviation to be the \code{qSd}-quantile.
}
\value{
    An object of class \linkS4class{BSseqStat}. More speciically, the input
    \linkS4class{BSseqStat} object with the computed statistics added to the
    \code{stats} slot (accessible with \code{\link{getStats}}).
}
\author{
    Kasper Daniel Hansen \email{khansen@jhsph.edu}
}

\seealso{
    \code{\link{BSmooth.fstat}} for the function to create the appropriate
    \code{\linkS4class{BSseqStat}} input object.
    \code{\linkS4class{BSseqStat}} also describes the return class.  This
    function is likely to be followed by the use of \code{\link{computeStat}}.}
\examples{
  \donttest{
    if(require(bsseqData)) {
        # library(limma) required for makeContrasts()
        library(limma)
        data(keepLoci.ex)
        data(BS.cancer.ex.fit)
        BS.cancer.ex.fit <- updateObject(BS.cancer.ex.fit)
        ## Remember to subset the BSseq object, see vignette for explanation
        ## TODO: Kind of a forced example
        design <- model.matrix(~0 + BS.cancer.ex.fit$Type)
        colnames(design) <- gsub("BS\\\.cancer\\\.ex\\\.fit\\\$Type", "",
                                 colnames(design))
        contrasts <- makeContrasts(
            cancer_vs_normal = cancer - normal,
            levels = design
        )
        BS.stat <- BSmooth.fstat(BS.cancer.ex.fit[keepLoci.ex,],
                                 design,
                                 contrasts)
        BS.stat <- smoothSds(BS.stat)
        ## Comparing the raw standard deviations to the smoothed standard
        ## deviations
        summary(getStats(BS.stat, what = "rawSds"))
        summary(getStats(BS.stat, what = "smoothSds"))
    }
  }
}
\keyword{internal}
