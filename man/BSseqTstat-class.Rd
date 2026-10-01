\name{BSseqTstat-class}
\Rdversion{1.1}
\docType{class}
\alias{BSseqTstat}
\alias{BSseqTstat-class}
\alias{[,BSseqTstat-method}
\alias{[,BSseqTstat,ANY,ANY,ANY-method}
\alias{show,BSseqTstat-method}
\alias{updateObject,BSseqTstat-method}
\title{Class BSseqTstat}
\description{
  A class for representing t-statistics for smoothed whole-genome
  bisulfite sequencing data.
}
\usage{
  BSseqTstat(gr = NULL, stats = NULL, parameters = NULL)
}
\arguments{
  \item{gr}{The genomic locations as an object of class \code{GRanges}.}
  \item{stats}{A matrix-like object (see 'Details' below) with columns
    representing various statistics for methylation loci along the
    genome.}
  \item{parameters}{A list of parameters, as an object of class \code{list},
    representing how the t-statistics were computed.}
}

\details{
    The 'stats' matrix-like object will be coerced to a
    \code{\link[DelayedArray]{DelayedMatrix}} objects; see
    \code{\link[DelayedArray]{DelayedMatrix}} in the \pkg{DelayedArray}
    package for the full list of supported matrix-like objects. We
    recommend using \code{\link[base]{matrix}} objects for in-memory
    storage of data and \code{\link[HDF5Array]{HDF5Matrix}} for on-disk
    storage of data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{BSseqTstat(...)}.
However, usually objects are returned by \code{BSmooth.tstat(...)} and
not constructed by the user..
}
\section{Extends}{
Class \code{\link[GenomicRanges]{GRanges}}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{The subsetting operator; one may only subset in one
      dimension, corresponding to methylation loci.}
    \item{show}{The show method.}
  }
}
\section{Utilities}{
  This class extends \code{\link[GenomicRanges]{GRanges}} and therefore
  inherits the full \code{GRanges} API. This means that \code{BSseqTstat}
  objects support any operation supported by \code{GRanges} objects, like
  \code{seqnames}, \code{start}, \code{end}, \code{width}, \code{strand},
  \code{[}, \code{findOverlaps}, \code{subsetByOverlaps}, etc...
  }
\section{Coercion}{
  Package version 1.11.1 introduced a new version of representing
  \sQuote{BSseqTstat} objects. You can update old serialized (saved)
  objects by invoking \code{x <- updateObject(x)}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  The package vignette(s).  \code{\link[GenomicRanges]{GRanges}} in the
  \pkg{GenomicRanges} package for accessing the genomic locations.
  \code{\link{BSmooth.tstat}} for a function
  that returns objects of class \code{BSseqTstat}, and \code{\link{dmrFinder}}
  for a function that computes DMRs based on the t-statistics.  Also see
  \code{\link[bsseqData]{BS.cancer.ex.tstat}} for an example of the
  class in the \pkg{bsseqData} package.
}
\examples{
if(require(bsseqData)) {
    data(BS.cancer.ex.tstat)
    BS.cancer.ex.tstat <- updateObject(BS.cancer.ex.tstat)
    dmrFinder(BS.cancer.ex.tstat)
}
}
\keyword{classes}
