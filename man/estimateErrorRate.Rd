\name{estimateErrorRate}
\alias{estimateErrorRate}
\title{
  Estimate CpG-specific error rate from BSseq object.
}
\description{
  This function estimates the CpG-specific error rate from a single sample BSseq object generated using \code{read.bedMethyl}.
}
\usage{
estimateErrorRate(BSseq, minCov = 10, maxCov = 100, minRatio = 0.8, plotErrorProfile = FALSE)
}
\arguments{
  \item{BSseq}{A single sample object of class \code{BSseq}.}
  \item{minCov}{A non-negative integer specifying the minimum coverage required for CpG loci to be considered.}
  \item{maxCov}{A non-negative integer specifying the maximum coverage allowed for CpG loci to be considered.}
  \item{minRatio}{A numeric value between 0 and 1 specifying the minimum ratio of CpG sites to non-CpG sites for a loci to be considered.}
  \item{plotErrorProfile}{A logical value indicating whether to plot the CpG to non-CpG ratio distribution for the filtered sites.}
}
\value{
  A numeric value representing the estimated CpG-specific error rate the BSseq object.
}
\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
\seealso{
  \code{\linkS4class{BSseq}} for the \code{BSseq} class, 
  \code{\link{read.bedMethyl}} for details on reading data into a \code{BSseq} object.
}
\examples{
# Example input files
infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
                         package = "bsseq"),
             system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
                         package = "bsseq"))

# Run the function to import data
bsseq <- read.bedMethyl(files = infiles,
                        colData = DataFrame(row.names = c("test_nanopore", 
                                                          "test_pacbio")),
                        strandCollapse = TRUE,
                        verbose = TRUE)

# Estimate error rate
estimateErrorRate(bsseq[, 1], plotErrorProfile = FALSE)
}
