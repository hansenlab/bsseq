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
    \item{stats}{The statistics, as a list of matrix-like objects (see
        'Details' below).}
    \item{parameters}{A list of parameters.}
}
\details{
  The matrix-like elements of the list in the 'stats' slot will be coerced to
  \linkS4class{DelayedMatrix} objects; see
  \code{?DelayedArray::\linkS4class{DelayedMatrix}} for the full list of
  supported matrix-like objects. We recommend using \link[base]{matrix} objects
  for in-memory storage of data and \linkS4class{HDF5Matrix} for on-disk
  storage of data.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{BSseqStat(...)}.
    However, usually objects are returned by \code{BSmooth.fstat(...)} and
    not constructed by the user.
}
\section{Slots}{
    \describe{
        \item{\code{stats}:}{This is a list of \linkS4class{DelayedMatrix}
            objects with list elements representing various statistics for
            methylation loci along the genome.}
        \item{\code{parameters}:}{Object of class \code{list}.  A list of
            parameters representing how the statistics were computed.}
        \item{\code{gr}:}{Object of class \code{GRanges} giving genomic
            locations.}
    }
}
\section{Extends}{
    Class \code{\linkS4class{hasGRanges}}, directly.
}
\section{Methods}{
    \describe{
        \item{[}{The subsetting operator; one may only subset in one
            dimension, corresponding to methylation loci.}
        \item{show}{The show method.}
    }
}
\section{Utilities}{
    This class extends \code{hasGRanges} and therefore inherits a number
    of useful \code{GRanges} methods that operate on the \code{gr} slot,
    used for accessing and setting the genomic locations and also do
    \code{subsetByOverlaps}.
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
    \code{\linkS4class{hasGRanges}} for accessing the genomic locations.
    \code{\link{BSmooth.fstat}} for a function
    that returns objects of class \code{BSseqStat}, and \code{\link{smoothSds}},
    \code{\link{computeStat}} and \code{\link{dmrFinder}}
    for functions that operate based on these statistics.  Also see
    the more specialised \code{\linkS4class{BSseqTstat}}.
}
\keyword{classes}
