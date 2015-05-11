\name{data.frame2GRanges}
\alias{data.frame2GRanges}
\title{
  Converts a data frame to a GRanges.
}
\description{
  Converting a data.frame to a GRanges object.  The data.frame needs
  columns like chr, start and end (strand is optional).  Additional
  columns may be kept in the GRanges object.
}
\usage{
data.frame2GRanges(df, keepColumns = FALSE, ignoreStrand = FALSE)
}
\arguments{
  \item{df}{A \code{data.frame} with columns \code{chr} or
    \code{seqnames}, \code{start}, \code{end} and optionally a
    \code{strand column}.}    
  \item{keepColumns}{In case \code{df} has additional columns, should
    these columns be stored as metadata columns on the return GRanges
    or should they be discarded.}
  \item{ignoreStrand}{In case \code{df} has a \code{strand} column,
    should this column be ignored.}
}
\note{
  In case \code{df} has \code{rownames}, they will be used as
  \code{names} for the return object.
}
\value{
  An object of class \code{GRanges}
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\examples{
df <- data.frame(chr = "chr1", start = 1:3, end = 2:4,
                 strand = c("+","-","+"))
data.frame2GRanges(df, ignoreStrand = TRUE)
}
