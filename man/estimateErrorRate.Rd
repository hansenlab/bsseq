\name{estimateErrorRate}
\alias{estimateErrorRate}
\title{
  Estimate sample specific error rate for CpG calls in a MethylCounts object.
}
\description{
  Estimates the CpG-specific error rate from a single-sample \code{MethylCounts} object. The input must be created using \code{read.bedMethyl} setting \code{output="MethylCounts"}.
}
\usage{
  estimateErrorRate(MethylCounts,
                    minCov = 10,
                    maxCov = 100,
                    minRatio = 0.8,
                    plotErrorProfile = FALSE)
}
\arguments{
  \item{MethylCounts}{A single sample object of class \code{MethylCounts}.}
  \item{minCov}{A non-negative integer specifying the minimum coverage required for CpG loci to be considered.}
  \item{maxCov}{A non-negative integer specifying the maximum coverage allowed for CpG loci to be considered.}
  \item{minRatio}{A numeric value between 0 and 1 specifying the minimum ratio of CpG sites to non-CpG sites for a locus to be considered.}
  \item{plotErrorProfile}{A logical value indicating whether to generate a plot showing the CpG-to-non-CpG ratio distribution for the filtered sites. This can help in visualizing filtering decisions and potential error sources.}
}
\value{
  A single numeric value representing the estimated CpG-specific error rate for the given sample. The error rate is calculated based on the filtered CpG loci in the \code{MethylCounts} object.
}
\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
\seealso{
  \code{\linkS4class{MethylCounts}} for the \code{MethylCounts} class,
  \code{\link{read.bedMethyl}} for details on reading data into a \code{MethylCounts} object.
}
\examples{
# Import example MethylCounts data
infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
                         package = "bsseq"),
             system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
                         package = "bsseq"))

# Read data into a MethylCounts object
mc <- read.bedMethyl(files = infiles,
                     colData = DataFrame(row.names = c("test_nanopore",
                                                       "test_pacbio")),
                     strandCollapse = TRUE,
                     output = "MethylCounts",
                     verbose = TRUE)

# Estimate error rate for the first sample
estimateErrorRate(mc[, 1], plotErrorProfile = FALSE)
}
