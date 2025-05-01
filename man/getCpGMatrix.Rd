\name{getCpGMatrix}
\alias{getCpGMatrix}
\title{
  Generate a matrix of the most likely CpG status for a multi-sample BSseq object.
}
\description{
  This function generates a matrix of the most likely CpG status for all loci and samples in a BSseq object. Each element of the matrix represents the most likely CpG status (0 for homozygous CpG, 1 for heterozygous CpG, and 2 for non-CpG or if allCpG = TRUE 0 for homozygous or heterozygous CpG, and 2 for non-CpG) for a specific locus and sample.
}
\usage{
getCpGMatrix(BSseq, e = NULL, allCpG = FALSE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{e}{An optional numeric vector representing error rates for each sample. If \code{NULL}, the error rate for each sample is estimated using \code{\link{estimateErrorRate}}.}
  \item{allCpG}{A logical value indicating whether to classify loci as allCpG (i.e. combine homozygous or heterozygous CpG) and non-CpG based on their likelihoods. Should be the same for \code{getCpGMatrix} and \code{getMaxLikelihoodMatrix} }
}
\value{
  A numeric matrix where each row represents a locus, and each column represents a sample, and the values correspond to the CpG status (same order as the BSseq object in input).
}
\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
\seealso{
  \code{\linkS4class{BSseq}} for the \code{BSseq} class, 
  \code{\link{read.bedMethyl}} for details on reading data into a \code{BSseq} object,
  \code{\link{estimateErrorRate}} for estimating the CpG-specific error rate.
  \code{\link{getCpGs}} for filtering a single-sample BSseg object.
  \code{\link{getMaxLikelihoodMatrix}} for generating a matrix with the maximum scaled likelihoods matching the CpGMatrix.  
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

# Single samples can be filtered using the getCpGs function
bsseq_nano <- bsseq[, 1]
bsseq_nano_99All_filtered <- bsseq[getCpGs(bsseq_nano, 
                                           type = "allCpG", threshold = 0.99)]

bsseq_pacbio <- bsseq[, 2]
bsseq_pacbio_99All_filtered <- bsseq[getCpGs(bsseq_pacbio, 
                                             type = "allCpG", threshold = 0.99)]

# For filtering multiple samples, we can use a CpGMatrix and a MaxLikelihoodMatrix
# Construct the CpGMatrix and getMaxLikelihoodMatrix for the bsseq object
CpGMatrix <- getCpGMatrix(bsseq, allCpG = TRUE)
MaxLikelihoodMatrix <- getMaxLikelihoodMatrix(bsseq, allCpG = TRUE)

# Filter for allCpG loci with a likelihood > 0.99 in both samples
bsseq_combined_99All_filtered <- bsseq[which(rowAlls(CpGMatrix == 0) 
                                          & rowMins(MaxLikelihoodMatrix) > 0.99)]
}
