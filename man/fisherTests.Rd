\name{fisherTests}
\alias{fisherTests}
\title{
  Compute Fisher-tests for a BSseq object
}
\description{
  A function to compute Fisher-tests for an object of class \code{BSseq}.
}
\usage{
fisherTests(BSseq, group1, group2, lookup = NULL,
  returnLookup = TRUE, mc.cores = 1, verbose = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{group1}{A vector of sample names or indexes for the
    \sQuote{treatment} group.}
  \item{group2}{A vector of sample names or indexes for the
    \sQuote{control} group.}
  \item{lookup}{A \sQuote{lookup} object, see details.}
  \item{returnLookup}{Should a \sQuote{lookup} object be returned, see
    details.}
  \item{mc.cores}{The number of cores used.  Note that setting
    \code{mc.cores} to a value greater than 1 is not supported on MS
    Windows, see the help page for \code{mclapply}.}
  \item{verbose}{Should the function be verbose.}
}
\details{
  This function computes row-wise Fisher's exact tests.  It uses an
  internal lookup table so rows which forms equivalent 2x2 tables are
  group together and only a single test is computed.  If
  \code{returnLookup} is \code{TRUE} the return object contains the
  lookup table which may be feed to another call to the function using
  the \code{lookup} argument.

  If \code{group1}, \code{group2} designates more than 1 sample, the
  samples are added together before testing.

  This function can use multiple cores on the same computer.

  This test cannot model biological variability.
}
\value{
  if \code{returnLookup} is \code{TRUE}, a list with components
  \code{results} and \code{lookup}, otherwise just the \code{results}
  component.  The \code{results} (component) is a matrix with the same
  number of rows as the \code{BSseq} argument and 2 columns
  \code{p.value} (the unadjusted p-values) and \code{log2OR} (log2
  transformation of the odds ratio).
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  \link[stats]{fisher.test} for information about Fisher's test.
  \code{\link{mclapply}} for the \code{mc.cores} argument.
}
\examples{
M <- matrix(1:9, 3,3)
colnames(M) <- c("A1", "A2", "A3")
BStest <- BSseq(pos = 1:3, chr = c("chr1", "chr2", "chr1"),
                M = M, Cov = M + 2)
results <- fisherTests(BStest, group1 = "A1", group2 = "A2",
                       returnLookup = TRUE)
results

#-------------------------------------------------------------------------------
# An example using a HDF5Array-backed BSseq object
#
library(HDF5Array)
# See ?SummarizedExperiment::saveHDF5SummarizedExperiment for details
hdf5_BStest <- saveHDF5SummarizedExperiment(x = BStest,
                                            dir = tempfile())
results <- fisherTests(hdf5_BStest,
                       group1 = "A1",
                       group2 = "A2",
                       returnLookup = TRUE)
results
}
