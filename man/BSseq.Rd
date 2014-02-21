\name{BSseq}
\alias{BSseq}
\title{
  The constructor function for BSseq objects.
}
\description{
  The constructor function for BSseq objects.
}
\usage{
BSseq(M = NULL, Cov = NULL, coef = NULL, se.coef = NULL,
  trans = NULL, parameters = NULL, pData = NULL, gr = NULL,
  pos = NULL, chr = NULL, sampleNames = NULL, rmZeroCov = FALSE)
}
\arguments{
  \item{M}{A matrix of methylation evidence.}
  \item{Cov}{A matrix of coverage.}
  \item{coef}{Smoothing estimates.}
  \item{se.coef}{Smoothing standard errors.}
  \item{trans}{A smoothing transformation.}
  \item{parameters}{A list of smoothing parameters.}
  \item{pData}{An \code{data.frame} or \code{DataFrame}.}
  \item{sampleNames}{A vector of sample names.}
  \item{gr}{An object of type \code{GRanges}.}
  \item{pos}{A vector of locations.}
  \item{chr}{A vector of chromosomes.}
  \item{rmZeroCov}{Should genomic locations with zero coverage in all
    samples be removed.}
}
\details{
  Genomic locations are specified either through \code{gr} or through
  \code{chr} and \code{pos} but not both.  There should be the same
  number of genomic locations as there are rows in the \code{M} and
  \code{Cov} matrix.

  The argument \code{rmZeroCov} may be useful in order to reduce the
  size of the return object be removing methylation loci with zero
  coverage.

  In case one or more methylation loci appears multiple times, the
  \code{M} and \code{Cov} matrices are summed over rows linked to the
  same methylation loci.  See the example below.

  Users should never have to specify \code{coef}, \code{se.coef},
  \code{trans}, and \code{parameters}, this is for internal use (they
  are added by \code{BSmooth}).

  \code{phenoData} is a way to specify pheno data (as known from the
  \code{ExpressionSet} and \code{eSet} classes), at a minimum
  \code{sampleNames} should be given (if they are not present, the
  function uses \code{col.names(M)}).  
}
\value{
  An object of class \code{BSseq}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  \code{\linkS4class{BSseq}}
}
\examples{
  M <- matrix(0:8, 3, 3)
  Cov <- matrix(1:9, 3, 3)
  BS1 <- BSseq(chr = c("chr1", "chr2", "chr1"), pos = c(1,2,3),
               M = M, Cov = Cov, sampleNames = c("A","B", "C"))
  BS1
  BS2 <- BSseq(chr = c("chr1", "chr1", "chr1"), pos = c(1,1,1),
               M = M, Cov = Cov, sampleNames = c("A","B", "C"))
  BS2
}
