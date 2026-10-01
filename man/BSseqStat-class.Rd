\name{BSseqStat-class}
\Rdversion{1.1}
\docType{class}
\alias{BSseqStat}
\alias{BSseqStat-class}
\alias{[,BSseqStat-method}
\alias{[,BSseqStat,ANY,ANY,ANY-method}
\alias{show,BSseqStat-method}
\alias{updateObject,BSseqStat-method}
\title{Class BSseqStat}
\description{
    A class for representing statistics for smoothed whole-genome
    bisulfite sequencing data.
}
\usage{
    BSseqStat(gr = NULL, stats = NULL, parameters = NULL)
}
\arguments{
    \item{gr}{The genomic locations as an object of class \code{GRanges}.}
    \item{stats}{A list of matrix-like objects (see 'Details' below) with
      list elements representing various statistics for methylation loci
      along the genome.}
    \item{parameters}{A list of parameters, as an object of class \code{list},
      representing how the statistics were computed.}
}
\details{
    The matrix-like elements of the list in the 'stats' slot will be
    coerced to \code{\link[DelayedArray]{DelayedMatrix}} objects; see
    \code{\link[DelayedArray]{DelayedMatrix}} in the \pkg{DelayedArray}
    package for the full list of supported matrix-like objects. We
    recommend using \code{\link[base]{matrix}} objects for in-memory
    storage of data and \code{\link[HDF5Array]{HDF5Matrix}} for on-disk
    storage of data.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{BSseqStat(...)}.
    However, usually objects are returned by \code{BSmooth.fstat(...)} and
    not constructed by the user.
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
    inherits the full \code{GRanges} API. This means that \code{BSseqStat}
    objects support any operation supported by \code{GRanges} objects, like
    \code{seqnames}, \code{start}, \code{end}, \code{width}, \code{strand},
    \code{[}, \code{findOverlaps}, \code{subsetByOverlaps}, etc...
}
\section{Coercion}{
  Package version 1.11.1 introduced a new version of representing
  \sQuote{BSseqStat} objects. You can update old serialized (saved)
  objects by invoking \code{x <- updateObject(x)}.
}
\author{
    Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
    \code{\link[GenomicRanges]{GRanges}} in the \pkg{GenomicRanges} package
    for accessing the genomic locations.
    \code{\link{BSmooth.fstat}} for a function
    that returns objects of class \code{BSseqStat}, and \code{\link{smoothSds}},
    \code{\link{computeStat}} and \code{\link{dmrFinder}}
    for functions that operate based on these statistics.  Also see
    the more specialised \code{\linkS4class{BSseqTstat}}.
}
\keyword{classes}
