\name{BSseq_mod}
\alias{BSseq_mod}
\title{
  The constructor function for modified BSseq objects.
}
\description{
  The constructor function for modified BSseq objects.
}
\usage{
BSseq_mod(M = NULL, U = NULL, Other_mod = NULL, Filtered = NULL, coef = NULL,
  se.coef = NULL, trans = NULL, parameters = NULL, pData = NULL, gr = NULL,
  pos = NULL, chr = NULL, sampleNames = NULL, rmZeroCov = FALSE)
}
\arguments{
  \item{M}{A matrix-like object of methylation evidence (see 'Details' below).}
  \item{U}{A matrix-like object of unmethylation evidence (see 'Details' below)).}
  \item{Other_mod}{A matrix-like object of other modification bases obtained
    from modkit when users choose to store it separately.}
  \item{Filtered}{A matrix-like object of ambiguous modification bases obtained
    from modkit.}
  \item{coef}{A matrix-like object of smoothing estimates (see 'Details'
    below).}
  \item{se.coef}{A matrix-like object of smoothing standard errors (see
    'Details' below).}
  \item{trans}{A smoothing transformation.}
  \item{parameters}{A list of smoothing parameters.}
  \item{pData}{An \code{data.frame} or \linkS4class{DataFrame}.}
  \item{sampleNames}{A vector of sample names.}
  \item{gr}{An object of type \linkS4class{GRanges}.}
  \item{pos}{A vector of locations.}
  \item{chr}{A vector of chromosomes.}
  \item{rmZeroCov}{Should genomic locations with zero coverage in all
    samples be removed.}
}
\details{
  The 'M', 'U', 'Other_mod', 'Filtered', 'coef', and 'se.coef' matrix-like
  objects will be coerced to \linkS4class{DelayedMatrix} objects; see
  \code{?DelayedArray::\linkS4class{DelayedMatrix}} for the full list of
  supported matrix-like objects. We recommend using \link[base]{matrix} objects
  for in-memory storage of data and \linkS4class{HDF5Matrix} for on-disk
  storage of data.

  Genomic locations are specified either through \code{gr} or through
  \code{chr} and \code{pos} but not both. There should be the same
  number of genomic locations as there are rows in the \code{M},
  \code{U}, and \code{Other_mod} matrix.

  The argument \code{rmZeroCov} may be useful in order to reduce the
  size of the return object be removing methylation loci with zero
  coverage.

  In case one or more methylation loci appears multiple times, the
  \code{M}, \code{U} and \code{Other_mod} matrices are summed over
  rows linked to the same methylation loci. See the example below.

  Users should never have to specify \code{coef}, \code{se.coef},
  \code{trans}, and \code{parameters}, this is for internal use (they
  are added by \code{BSmooth}).

  \code{phenoData} is a way to specify pheno data (as known from the
  \code{ExpressionSet} and \code{eSet} classes), at a minimum
  \code{sampleNames} should be given (if they are not present, the
  function uses \code{col.names(M)}).
}
\value{
  An object of class \code{BSseq_mod}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  \code{\linkS4class{BSseq_mod}}
}
\examples{
  M <- matrix(0:8, 3, 3)
  U <- matrix(1:9, 3, 3)
  Other_mod <- matrix(1:9, 3, 3)
  BS1 <- BSseq_mod(chr = c("chr1", "chr2", "chr1"), pos = c(1,2,3),
                    M = M, U = U, Other_mod = Other_mod,
                    sampleNames = c("A","B", "C"))
  BS1
  BS2 <- BSseq_mod(chr = c("chr1", "chr1", "chr1"), pos = c(1,1,1),
                    M = M, U = U, Other_mod = Other_mod,
                    sampleNames = c("A","B", "C"))
  BS2

  #-------------------------------------------------------------------------------
  # An example using a HDF5Array-backed BSseq object
  #

  library(HDF5Array)
  hdf5_M <- realize(M, "HDF5Array")
  hdf5_U <- realize(U, "HDF5Array")
  hdf5_Other_mod <- realize(Other_mod, "HDF5Array")
  hdf5_BS1 <- BSseq_mod(chr = c("chr1", "chr2", "chr1"),
                        pos = c(1, 2, 3),
                        M = hdf5_M,
                        U = hdf5_U,
                        Other_mod = hdf5_Other_mod,
                        sampleNames = c("A", "B", "C"))
  hdf5_BS1
  hdf5_BS2 <- BSseq_mod(chr = c("chr1", "chr1", "chr1"),
                        pos = c(1, 1, 1),
                        M = hdf5_M,
                        U = hdf5_U,
                        Other_mod = hdf5_Other_mod,
                        sampleNames = c("A", "B", "C"))
  hdf5_BS2
}
