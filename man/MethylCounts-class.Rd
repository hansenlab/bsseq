\name{MethylCounts-class}
\Rdversion{1.1}
\docType{class}
\alias{MethylCounts-class}
\alias{show,MethylCounts-method}
\alias{getMethylCounts}
\alias{MethylCounts}
\title{Class MethylCounts}
\description{
  A class for representing whole-genome or capture sequencing
  methylation data.
}
\section{Objects from the Class}{
  An object from the class links together several pieces of information.
  (1) genomic locations stored as a \code{GRanges} object, a location by
  samples matrix of M values, a location by samples matrix of Cov
  (coverage) values, a location by samples matrix of Filtered (ambiguous
  modification status) values, and phenodata information.In addition,
  there are slots for representing smoothed data. This class is an
  extension of \link[SummarizedExperiment]{RangedSummarizedExperiment}
  from the \pkg{SummarizedExperiment} package.
}
\section{Slots}{
  In addition to the slots from \code{RangedSummarizedExperiment} we have
  \describe{
    \item{\code{parameters}:}{Object of class \code{list}.  A list of
      parameters representing for example how the data was smoothed.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{The show method.}
  }
}
\usage{
  MethylCounts(M = NULL, U = NULL, H = NULL, D = NULL, 
               gr = NULL, pos = NULL, chr = NULL,
               colData = NULL, sampleNames = NULL,
               parameters = NULL, rmZeroCov = FALSE)

  getMethylCounts(MethylCounts,
                  type = c("M", "U", "H", "D", "Cov", "gr", "parameters"),
                  withDimnames = TRUE)
}
\arguments{
  \item{M}{A matrix-like object of M values (see 'Details' below).}
  \item{U}{A matrix-like object of U values (see 'Details' below).}
  \item{H}{An optional matrix-like object of H values (see 'Details'
    below).}
  \item{D}{An optional matrix-like object of D values (see 'Details' below).}
  \item{gr}{An object of type \code{\link[GenomicRanges]{GRanges}}.}
  \item{pos}{A vector of locations.}
  \item{chr}{A vector of chromosomes.}
  \item{colData}{Typically a \code{DataFrame} with sample
    coveriates. See the documentation for \code{SummarizedExperiment}.}
  \item{sampleNames}{A character vector of sample names.}
  \item{parameters}{A list of parameters, for record-keeping.}
  \item{rmZeroCov}{\code{logical(1)}. Remove loci with zero coverage}
  \item{MethylCounts}{An object of class \code{MethylCounts}.}
  \item{type}{What should be extracted from the object?}
  \item{withDimnames}{\code{logical(1)}. Should dimnames be added to the
                  output when relevent?} 
}
\section{Details}{
  The \code{M}, \code{U}, \code{H}, and \code{D} matrix-like objects will be
  coerced to \code{\link[DelayedArray]{DelayedMatrix}} objects; see
  \code{\link[DelayedArray]{DelayedMatrix}} in the \pkg{DelayedArray}
  package for the full list of supported matrix-like objects. We
  recommend using \link[base]{matrix} objects for in-memory storage of
  data and \code{\link[HDF5Array]{HDF5Matrix}} for on-disk storage of
  data. These objects should have the same dimension and be integer-values.

  \code{M}, \code{U} are required methylation and 'unmethylation'
  values. \code{H} is optional hydroxymethylation values. \code{D}
  represents addition coverage for a loci, but where the data supports
  neither M, U nor H.
  
  The argument \code{rmZeroCov} may be useful in order to reduce the
  size of the return object be removing methylation loci with zero
  coverage.

  \code{pos} and \code{chr} can be supplied instead of \code{gr}.
}
\section{Utilities}{
  This class extends
  \link[SummarizedExperiment]{RangedSummarizedExperiment}
  from the \pkg{SummarizedExperiment} package and therefore
  inherits a number of useful \code{GRanges} methods that operate on the
  \code{rowRanges} slot, used for accessing and setting the genomic locations
  and also do \code{subsetByOverlaps}.

  There are a number of almost methods-like functions for operating on
  objects of class \code{MethylCounts}, including
  \code{getMethylCounts}. They are detailed below.

  \describe{
    \item{\code{getMethylCounts(MethylCounts, type, withDimnames)}}{
      A general accessor: is used to obtain a specific slot of an object of
      class \code{MethylCounts}.}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  The package vignette.
  \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
  (in the \pkg{SummarizedExperiment} package) for the underlying class.
}
\examples{
M <- matrix(1:9, 3,3)
colnames(M) <- c("A1", "A2", "A3")
MCtest <- MethylCounts(pos = 1:3,
                       chr = c("chr1", "chr2", "chr1"),
                       M = M, U = M + 2)

#-------------------------------------------------------------------------------
# An example using a HDF5-backed MethylCounts object
#
hdf5_NCtest <- realize(MCtest, "HDF5Array")
}
\keyword{classes}
