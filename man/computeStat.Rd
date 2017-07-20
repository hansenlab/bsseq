\name{computeStat}
\alias{computeStat}
\title{
    Compute a test statistic based on smoothed whole-genome bisulfite
    sequencing data.
 }
\description{
    Compute a test statistic based on smoothed whole-genome bisulfite
    sequencing data.
 }
\usage{
    computeStat(BSseqStat, coef = NULL)
}

\arguments{
    \item{BSseqStat}{An object of class \code{BSseqStat}, typically an object
        returned by \code{\link{smoothSds}(...)} and not constructed by
        the user.}
    \item{coef}{A vector indicating for which coefficients the statistic is to
        be computed (\code{coef = NULL} corresponds to testing all
        coefficients). If the length of the \code{coef} is 1 then the
        corresponding t-statistic is computed, otherwise the corresponding
        F-statistic is computed.}
}
\details{
    \strong{TODO}
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
  \code{\link{smoothSds}} for the function to create the appropriate
  \code{\linkS4class{BSseqStat}} input object.
  \code{\linkS4class{BSseqStat}} also describes the return class.  This
  function is likely to be followed by the use of \code{\link{dmrFinder}}.}
\examples{
  \donttest{
    if(require(bsseqData)) {
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
        BS.stat <- computeStat(BS.stat)
        BS.stat
    }
  }
}
\keyword{internal}
