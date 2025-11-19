\name{getCpGs}
\alias{getCpGs}
\title{
  Get CpG (and non-CpG) loci from a single sample MethylCounts object.
}
\description{
  This function identifies CpG (and non-CpG) loci from a single sample MethylCounts object using scaled likelihoods computed from the count of CpG sites and Non-CpG sites mapped to the loci, based on the specified type, and minimum scaled likelihood (threshold). The function only work for MethylCounts objects generated using read.bedMethyl.
}
\usage{
  getCpGs(MethylCounts,
          type = c("homozygous", "heterozygous", "allCpG", "nonCpG"),
          threshold = 0.99,
          e = NULL)
}
\arguments{
  \item{MethylCounts}{An single sample object of class \code{MethylCounts}.}
  \item{type}{A character string specifying the type of loci to extract. Must be one of \code{"homozygous"}, \code{"heterozygous"}, \code{"allCpG"}, or \code{"nonCpG"}.}
  \item{threshold}{A numeric value between 0 and 1 specifying the minimum likelihood threshold required for loci to be included.}
  \item{e}{An optional numeric value representing the error rate. If \code{NULL}, the error rate is estimated using \code{\link{estimateErrorRate}}.}
}
\value{
  An integer vector of indices representing the loci that match the criteria.
}
\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
\seealso{
  \code{\linkS4class{MethylCounts}} for the \code{MethylCounts} class,
  \code{\link{read.bedMethyl}} for details on reading data into a \code{MethylCounts} object,
  \code{\link{estimateErrorRate}} for estimating the CpG-specific error rate.
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

# Filter CpG sites for the Nanopore dataset
mc_nano <- mc[, 1]
mc_nano_99All_filtered <- mc_nano[getCpGs(mc_nano,
                                           type = "allCpG", threshold = 0.99)]

# Filter CpG sites for the PacBio dataset
mc_pacbio <- mc[, 2]
mc_pacbio_99All_filtered <- mc_pacbio[getCpGs(mc_pacbio,
                                             type = "allCpG", threshold = 0.99)]
}
