\name{FWGRanges-class}
\Rdversion{1.1}
\docType{class}
\alias{FWIRanges-class}
\alias{end,FWIRanges-method}
\alias{end<-,FWIRanges-method}
\alias{names,FWIRanges-method}
\alias{names<-,FWIRanges-method}
\alias{start,FWIRanges-method}
\alias{start<-,FWIRanges-method}
\alias{width,FWIRanges-method}
\alias{width<-,FWIRanges-method}
\alias{FWGRanges-class}
\alias{findOverlaps,FWGRanges,FWGRanges-method}
\alias{end,FWGRanges-method}
\alias{start,FWGRanges-method}
\alias{strand,FWGRanges-method}
\alias{seqnames,FWGRanges-method}
\alias{width,FWGRanges-method}
\alias{dimnames,arrayRealizationSink-method}
\title{Classes FWIRanges and FWGRanges}
\description{
  Classes for fixed-width IRanges and GRanges, ie. objects where all
  ranges have the same width. The intention is for these classes to be
  added to GenomicRanges. Documented here temporarily.
}
\details{
  See description. Otherwise works like IRanges and GRanges, except
  there are many unimplemented methods.

  This is really a private class, with private methods, but R's
  NAMESPACE handling means they get unintentionally exported. Hence this
  documentation. 
}
\examples{
showClass("FWIRanges")
}
\keyword{classes}
