\name{read.modkit}
\alias{read.modkit}
\title{Construct BSseq objects from nanopore BED files}
\description{
Construct BSseq objects from nanopore BED files
}
\usage{
read.modkit(
    files,
    colData = NULL,
    rmZeroCov = FALSE,
    strandCollapse = TRUE
)
}
\arguments{
    \item{files}{vector, BED files}
    \item{colData}{data frame, phenotypic data with samples as rows and variables as columns}
    \item{rmZeroCov}{A logical (1) indicating whether methylation loci that have zero coverage in all samples be removed}
    \item{strandCollapse}{A logical (1) indicating whether stand-symmetric methylation loci (i.e. CpGs) should be collapsed across strands}
}

\value{
BSseq objects
}
\details{
This function reads in nanopore sequencing modified BED files
to Bsseq objects. Nanopore sequencing data (i.e. aggregated modified base
counts) is stored in modified-base BAM files. These modified-base BAM files
are converted to bedMethyl (BED) files using \href{https://github.com/nanoporetech/modkit}{modkit}.

\subsection{Details for modkit}{
    Modkit outputs modified reads, unmodified reads, ambiguous modification reads (reads where the probability was below the threshold and usually failing the lowest 10th percentile), and other modified reads.}

\subsection{modkit to Bsseq object}{
    After creating BED files using modkit, the BED files are read in and the Bsseq object is constructed via \code{read.modkit()} function. The function reads in BED files, extract genomic regions, methylation, coverage, ambiguous modification status data and sample information and then construct Bsseq object using \code{BSseq} function within the package. Other modification bases such as hydroxymethylation are extracted and added to the methylation matrix when present.}
}

\examples{
# No other modification present
files <- c(system.file("extdata/modkit/chr21.chr22.HG002.top1000.bed.gz", package = "bsseq"))
bsseq_nano <- read.modkit(files, rmZeroCov = FALSE, strandCollapse=FALSE)

# Other modification present
files <- c(system.file("extdata/modkit/Hypo1.first50Bed.txt",package = "bsseq"))
bsseq_nano <- read.modkit(files, rmZeroCov = FALSE, strandCollapse=FALSE)
}
