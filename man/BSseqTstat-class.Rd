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
  \item{stats}{The statistics, as a matrix-like object (see 'Details' below).}
  \item{parameters}{A list of parameters.}
}

\details{
  The 'stats' matrix-like object will be coerced to a
  \linkS4class{DelayedMatrix} object; see
  \code{?DelayedArray::\linkS4class{DelayedMatrix}} for the full list of
  supported matrix-like objects. We recommend using \link[base]{matrix} objects
  for in-memory storage of data and \linkS4class{HDF5Matrix} for on-disk
  storage of data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{BSseqTstat(...)}.
However, usually objects are returned by \code{BSmooth.tstat(...)} and
not constructed by the user..
}
\section{Slots}{
  \describe{
    \item{\code{stats}:}{This is a \linkS4class{DelayedMatrix} object with
      columns representing various statistics for methylation loci along the
      genome.}
    \item{\code{parameters}:}{Object of class \code{list}.  A list of
      parameters representing how the t-statistics were computed.}
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
  \sQuote{BSseqTstat} objects. You can update old serialized (saved)
  objects by invoking \code{x <- updateObject(x)}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  The package vignette(s).  \code{\linkS4class{hasGRanges}} for accessing
  the genomic locations.  \code{\link{BSmooth.tstat}} for a function
  that returns objects of class \code{BSseqTstat}, and \code{\link{dmrFinder}}
  for a function that computes DMRs based on the t-statistics.  Also see
  \code{\link[bsseqData]{BS.cancer.ex.tstat}} for an example of the
  class in the \pkg{bsseqData} package.
}
\keyword{classes}
