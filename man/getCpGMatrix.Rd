\name{getCpGMatrix}
\alias{getCpGMatrix}
\title{
  Generate a matrix of the most likely CpG status for a multi-sample MethylCounts object.
}
\description{
  This function generates a matrix of the most likely CpG status for all loci and samples in a MethylCounts object. Each element of the matrix represents the most likely CpG status (0 for homozygous CpG, 1 for heterozygous CpG, and 2 for non-CpG or if allCpG = TRUE 0 for homozygous or heterozygous CpG, and 2 for non-CpG) for a specific locus and sample.
}
\usage{
  getCpGMatrix(MethylCounts,
               e = NULL,
               allCpG = FALSE)
}
\arguments{
  \item{MethylCounts}{An object of class \code{MethylCounts}.}
  \item{e}{An optional numeric vector representing error rates for each sample. If \code{NULL}, the error rate for each sample is estimated using \code{\link{estimateErrorRate}}.}
  \item{allCpG}{A logical value indicating whether to classify loci as allCpG (i.e. combine homozygous or heterozygous CpG) and non-CpG based on their likelihoods. Should be the same for \code{getCpGMatrix} and \code{getMaxLikelihoodMatrix} }
}
\value{
  A numeric matrix where each row represents a locus, and each column represents a sample, and the values correspond to the CpG status (same order as the MethylCounts object in input).
}
\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
\seealso{
  \code{\linkS4class{MethylCounts}} for the \code{MethylCounts} class,
  \code{\link{read.bedMethyl}} for details on reading data into a \code{MethylCounts} object,
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
mc <- read.bedMethyl(files = infiles,
                        colData = DataFrame(row.names = c("test_nanopore",
                                                          "test_pacbio")),
                        strandCollapse = TRUE,
                        output = "MethylCounts",
                        verbose = TRUE)

# Single samples can be filtered using the getCpGs function
mc_nano <- mc[, 1]
mc_nano_99All_filtered <- mc_nano[getCpGs(mc_nano,
                                             type = "allCpG",
                                             threshold = 0.99)]

mc_pacbio <- mc[, 2]
mc_pacbio_99All_filtered <- mc_pacbio[getCpGs(mc_pacbio,
                                              type = "allCpG",
                                              threshold = 0.99)]

# For filtering multiple samples, we can use a CpGMatrix and a MaxLikelihoodMatrix
# Construct the CpGMatrix and getMaxLikelihoodMatrix for the MethylCounts object

CpGMatrix <- getCpGMatrix(mc, allCpG = TRUE)
MaxLikelihoodMatrix <- getMaxLikelihoodMatrix(mc, allCpG = TRUE)

# Filter for allCpG loci with a likelihood > 0.99 in both samples
mc_combined_99All_filtered <- mc[which(rowAlls(CpGMatrix == 0) &
                                       rowMins(MaxLikelihoodMatrix) > 0.99)]
}
