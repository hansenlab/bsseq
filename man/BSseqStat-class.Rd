\name{BSseqStat-class}
\Rdversion{1.1}
\docType{class}
\alias{BSseqStat}
\alias{BSseqStat-class}
\alias{[,BSseqStat-method}
\alias{[,BSseqStat,ANY-method}
\alias{show,BSseqStat-method}
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
    \item{stats}{The statistics, as a matrix.}
    \item{parameters}{A list of parameters.}
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{BSseqStat(...)}.
    However, usually objects are returned by \code{BSmooth.fstat(...)} and
    not constructed by the user.
}
\section{Slots}{
    \describe{
        \item{\code{stats}:}{This is a matrix with columns representing
            various statistics for methylation loci along the genome.}
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
\author{
    Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
    The package vignette(s).  \code{\linkS4class{hasGRanges}} for accessing
    the genomic locations.  \code{\link{BSmooth.fstat}} for a function
    that objects of class \code{BSseqStat}, and \code{\link{dmrFinder}}
    for a function that computes DMRs based on the statistics.  Also see
    the mpre specialised \code{\linkS4class{BSseqTstat}}.
}
\keyword{classes}
