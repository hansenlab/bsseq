\name{BSmooth.fstat}
\alias{BSmooth.fstat}
\title{
    Compute F-statistics based on smoothed whole-genome bisulfite
    sequencing data.
}
\description{
    Compute F-statistics based on smoothed whole-genome bisulfite
    sequencing data.
}
\usage{
BSmooth.fstat(BSseq, design, contrasts, verbose = TRUE)
}
\arguments{
    \item{BSseq}{An object of class \code{BSseq}.}
    \item{design}{The design matrix of the bisulfite-sequencing experiment,
        with rows corresponding to arrays and columns to coefficients to be
        estimated.}
    \item{contrasts}{Numeric matrix with rows corresponding to columns in
        \code{design} and columns containing contrasts. May be a vector if
        there is only one contrast.}
    \item{verbose}{Should the function be verbose?}
}
\details{
    \strong{TODO}
}
\value{
    An object of class \linkS4class{BSseqStat}.
}
\author{
    Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
    \code{\link{BSmooth}} for the input object and
    \linkS4class{BSseq} for its class.
    \linkS4class{BSseqStat} describes the return class. This
    function is likely to be followed by the use of
    \code{\link{smoothSds}}, \code{\link{computeStat}}, and
    \code{\link{dmrFinder}}.
}
\examples{
\donttest{
    if(require(bsseqData)) {
    # limma required for makeContrasts()
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
    BS.stat

    #---------------------------------------------------------------------------
    # An example using a HDF5Array-backed BSseq object
    #
    library(HDF5Array)
    # See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
    hdf5_BS.cancer.ex.fit <- saveHDF5SummarizedExperiment(
        x = BS.cancer.ex.fit[keepLoci.ex, ],
        dir = tempfile())
    hdf5_BS.stat <- BSmooth.fstat(hdf5_BS.cancer.ex.fit,
                                  design,
                                  contrasts)
    hdf5_BS.stat
    }
}
}
\keyword{internal}
