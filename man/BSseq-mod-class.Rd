\name{BSseq-mod-class}
\Rdversion{1.1}
\docType{class}
\alias{BSseq_mod-class}
\alias{[,BSseq_mod-method}
\alias{combine,BSseq_mod,BSseq_mod-method}
\alias{length,BSseq_mod-method}
\alias{pData,BSseq_mod-method}
\alias{pData<-,BSseq_mod,data.frame-method}
\alias{pData<-,BSseq_mod,DataFrame-method}
\alias{sampleNames,BSseq_mod-method}
\alias{sampleNames<-,BSseq_mod,ANY-method}
\alias{updateObject,BSseq_mod-method}
\alias{assays,BSseq_mod-method}
\alias{assayNames,BSseq_mod-method}
\alias{show,BSseq_mod-method}
\alias{getBSseq_mod}
\alias{strandCollapse_mod}
\alias{orderBSseq_mod}
\alias{chrSelectBSseq_mod}
\alias{hasBeenSmoothed_mod}

\title{Class BSseq_mod}
\description{
  A modified class for representing Nanopore data created by modkit to
  store other modification calls as a separate matrix based on users' choice.
}
\section{Objects from the Class}{
An object from the class links together several pieces of information.
(1) genomic locations stored as a \code{GRanges} object, a location by
samples matrix of M values, a location by samples matrix of U
(unmethylated) values, a location by samples matrix of Other_mod (other
modification calls), a location by samples matrix of Filtered (ambiguous
modification status) values, and phenodata information. In addition, there are
slots for representing smoothed data. This class is an extension of
\linkS4class{RangedSummarizedExperiment}.
}
\section{Slots}{
  \describe{
    \item{\code{trans}:}{Object of class \code{function}.  This
      function transforms the \code{coef} slot from the scale the
      smoothing was done to the 0-1 methylation scale.}
    \item{\code{parameters}:}{Object of class \code{list}.  A list of
      parameters representing for example how the data was smoothed.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Modified BSseq")}: Subsetting by location
      (using integer indices) or sample (using integers or sample
      names).}
    \item{length}{Unlike for \code{RangedSummarizedExperiment},
      \code{length()} is the number of methylation loci (equal to
      \code{length(granges(x))}).}
    \item{sampleNames,sampleNames<-}{Sample names and its replacement
      function for the object.  This is an alias for \code{colnames}.}
    \item{pData,pData<-}{Obtain and replace the \code{pData} slot of the
      \code{phenoData} slot. This is an alias for \code{colData}.}
    \item{show}{The show method.}
    \item{combine}{This function combines two \code{BSseq_mod} objects.
      The genomic locations of the new object is the union of the
      genomic locations of the individual objects.  In addition, the
      methylation data matrices are placed next to each other (as
      appropriate wrt. the new genomic locations) and zeros are entered
      into the matrices as needed.}
}}
\section{Utilities}{
  This class extends \linkS4class{RangedSummarizedExperiment} and therefore
  inherits a number of useful \code{GRanges} methods that operate on the
  \code{rowRanges} slot, used for accessing and setting the genomic locations
  and also do \code{subsetByOverlaps}.

  There are a number of almost methods-like functions for operating on
  objects of class \code{BSseq_mod}, including \code{getBSseq_mod},
  \code{collapseBSseq_mod}, and \code{orderBSseq_mod}.  They are detailed below.

  \describe{
  \item{\code{collapseBSseq_mod(BSseq_mod, columns)}}{

      is used to collapse an object of class \code{BSseq_mod}.  By
      collapsing we simply mean that certain columns (samples) are merge
      together by summing up the methylation, unmethylation, and other
      modification evidence. This is a useful function if you start by
      reading in a dataset based on say flowcells and you (after QC) want
      to simply add a number of flowcells into one sample. The argument
      \code{columns} specify which samples are to be merged, in the
      following way: it is a character vector of new sample names, and the
      names of the column vector indicates which samples in the
      \code{BSseq_mod} object are to be collapsed.  If \code{columns} have
      the same length as the number of rows of \code{BSseq_mod} (and has no
      names) it is assumed that the ordering corresponds to the sample
      ordering in \code{BSseq_mod}.}

    \item{\code{orderBSseq_mod(BSseq_mod, seqOrder = NULL)}}{

      simply orders an object of class \code{BSseq_mod} according to
      (increasing) genomic locations. The \code{seqOrder} vector is a
      character vector of \code{seqnames(BSseq_mod)} describing the order of
      the chromosomes. This is useful for ordering \code{chr1} before
      \code{chr10}.}

    \item{\code{chrSelectBSseq_mod(BSseq_mod, seqnames = NULL, order = FALSE)}}{

      subsets and optionally reorders an object of class \code{BSseq_mod}.
      The \code{seqnames} vector is a character vector of
      \code{seqnames(BSseq_mod)} describing which chromosomes should be
      retained.  If \code{order} is \code{TRUE}, the chromosomes are
      also re-ordered using \code{orderBSseq_mod}.}

    \item{\code{getBSseq_mod(BSseq_mod, type = c("U", "M", "Other_mod", "Filtered",
    "gr", "coef", "se.coef", "trans",  "parameters"))}}{

      is a general accessor: is used to obtain a specific slot of an
      object of class \code{BSseq_mod}.  It is primarily intended for
      internal use in the package, for users we recommend \code{granges}
      to get the genomic locations, \code{getCoverage_mod} to get the
      coverage slots and \code{getMeth_mod} to get the smoothed values (if
      they exist).  }

    \item{\code{hasBeenSmoothed_mod(BSseq_mod)}}{

      This function returns a logical depending on whether or not the
      \code{BSseq_mod} object has been smoothed using \code{BSmooth}.  }

    \item{\code{strandCollapse_mod(BSseq_mod, shift = TRUE)}}{

      This function operates on a \code{BSseq_mod} objects which has
      stranded loci (i.e. loci where the strand is one of \sQuote{+} or
      \sQuote{-}). It will collapse the methylation and coverage
      information across the two strands, unstranding the loci in the process
      and potentially re-ordering them.

      The argument \code{shift} indicates whether the positions for the loci
      on the reverse strand should be shifted one (i.e. the positions for
      these loci are the positions of the \sQuote{G} in the
      \sQuote{CpG}; this is the case for Bismark output for example).}
  }
}
\section{Coercion}{
  Package versions 1.5.2 and 1.11.1 introduced a new version of representing
  \sQuote{BSseq} objects.  You can update old serialized (saved)
  objects by invoking \code{x <- updateObject(x)}.
}
\section{Assays}{
  This class overrides the default implementation of \code{assays} to
  make it faster.  Per default, no names are added to the returned data
  matrices.

  Assay names can conveniently be obtained by the function
  \code{assayNames(x)}
}

\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{

  The package vignette.  \code{\link{BSseq_mod}} for the constructor
  function.  \code{\linkS4class{RangedSummarizedExperiment}}
  for the underlying class.  \code{\link{getBSseq_mod}},
  \code{\link{getCoverage_mod}}, and \code{\link{getMeth_mod}} for accessing the
  data stored in the object and finally \code{\link{BSmooth}} for
  smoothing the bisulfite sequence data.

}
\examples{
M <- U <- Other_mod <- matrix(1:9, 3,3)
colnames(M) <- c("A1", "A2", "A3")
BStest <- BSseq_mod(pos = 1:3, chr = c("chr1", "chr2", "chr1"), M = M, U = U,
                    Other_mod = Other_mod)
}
\keyword{classes}
