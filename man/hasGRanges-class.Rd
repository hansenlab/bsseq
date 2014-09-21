\name{hasGRanges-class}
\Rdversion{1.1}
\docType{class}
\alias{hasGRanges-class}
\alias{[,hasGRanges-method}
\alias{[,hasGRanges,ANY-method}
\alias{end,hasGRanges-method}
\alias{end<-,hasGRanges-method}
\alias{granges,hasGRanges-method}
\alias{length,hasGRanges-method}
\alias{seqlengths,hasGRanges-method}
\alias{seqlengths<-,hasGRanges-method}
\alias{seqlevels,hasGRanges-method}
\alias{seqlevels<-,hasGRanges-method}
\alias{seqnames,hasGRanges-method}
\alias{seqnames<-,hasGRanges-method}
\alias{start,hasGRanges-method}
\alias{start<-,hasGRanges-method}
\alias{strand,hasGRanges-method}
\alias{strand<-,hasGRanges-method}
\alias{strand<-,hasGRanges,ANY-method}
\alias{subsetByOverlaps,GenomicRanges,hasGRanges-method}
\alias{subsetByOverlaps,hasGRanges,GenomicRanges-method}
\alias{subsetByOverlaps,hasGRanges,hasGRanges-method}
\alias{findOverlaps,GenomicRanges,hasGRanges-method}
\alias{findOverlaps,hasGRanges,GenomicRanges-method}
\alias{findOverlaps,hasGRanges,hasGRanges-method}
\alias{width,hasGRanges-method}
\alias{width<-,hasGRanges-method}

\title{Class hasGRanges}
\description{
  A class with a GRanges slot, used as a building block for other
  classes.  Provides basic accessor functions etc.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hasGRanges", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{gr}:}{Object of class \code{GRanges}.}
  }
}
\section{Methods}{
  \describe{
    \item{"["}{Subsets a single dimension.}
    \item{granges}{Get the \code{GRanges} object representing 
      genomic locations.}
    \item{start,start<-,end,end<-,width,width<-}{Start, end and width
      for the genomic locations of the object, also replacement
      functions.  This accessor functions operate directly on the
      \code{gr} slot.} 
    \item{strand,strand<-}{Getting and setting the
      \code{strand} of the genomic locations (the \code{gr} slot).} 
    \item{seqlengths,seqlengths<-}{Getting and setting the
      \code{seqlengths} of the genomic locations (the \code{gr} slot).} 
    \item{seqlevels,seqlevels<-}{Getting and setting the
      \code{seqlevels} of the genomic locations (the \code{gr} slot).} 
    \item{seqnames,seqnames<-}{Getting and setting the
      \code{seqnames} of the genomic locations (the \code{gr} slot).} 
    \item{show}{The show method.}
    \item{findOverlaps}{\code{(query = "hasGRanges", subject =
	"hasGRanges")}: finds overlaps between the \code{granges()} of
      the two objects.}
    \item{findOverlaps}{\code{(query = "GenomicRanges", subject =
	"hasGRanges")}: finds overlaps between query and the \code{granges()} of
      the subject.}
    \item{findOverlaps}{\code{(query = "hasGRanges", subject =
	"GenomicRanges")}: finds overlaps between the \code{granges()}
      of the query and the subject.}
    \item{subsetByOverlaps}{\code{(query = "hasGRanges", subject =
	"hasGRanges")}: Subset the query, keeping the genomic locations that 
      overlaps the subject.}
    \item{subsetByOverlaps}{\code{(query = "hasGRanges", subject =
	"GenomicRanges")}: Subset the query, keeping the genomic
      locations that overlaps the subject.}
    \item{subsetByOverlaps}{\code{(query = "GenomicRanges",
	subject = "hasGRanges")}: Subset the query, keeping the genomic
      locations that overlaps the subject.}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\note{
  If you extend the \code{hasGRanges} class, you should consider writing
  a subset method (\code{[}), and a \code{show} method.  If the new
  class supports single index subsetting, the \code{subsetByOverlaps}
  methods show extend without problems.
}
\examples{
showClass("hasGRanges")
}
\keyword{classes}
