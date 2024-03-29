\name{read.modbam2bed}
\alias{read.modbam2bed}
\title{Construct BSseq objects from nanopore BED files}
\description{
Construct BSseq objects from nanopore BED files
}
\usage{
read.modbam2bed(
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
BSseq object
}
\details{
This function reads in nanopore sequencing modified BED files
to Bsseq objects. Nanopore sequencing data (i.e. aggregated modified base
counts) is stored in modified-base BAM files. These modified-base BAM files
are converted to bedMethyl (BED) files using \href{https://github.com/epi2me-labs/modbam2bed}{modbam2bed}.

\subsection{Details for using modbam2bed}{
    After installing modbam2bed, a conda environment is activated. Index files
for BAM files are created using \code{samtools index}. The code requires aligned reads with the Mm and Ml tags (MM and ML also supported), and the reference sequence used for alignment (<reference.fasta>).

\itemize{
    \item \code{-e, -- extended} to output canonical, modified, and filtered bases;
    \item \code{-m, -- mod_base=BASE} to output modified base of interest, one of: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao. (Or modA, modC, modG, modT, modU, modN for generic modified base);
    \item \code{-r, --region=chr:start-end} to output chromosome or genomic region of interest;
    \item \code{-f, --threshold=THRESHOLD} to output filtered bases for
probability lower than threshold (default = 0.66)
}
}

  \subsection{modbam2bed to Bsseq object}{
    After creating BED files using modbam2bed, the BED files are read in and the Bsseq object is constructed via \code{read.modbam2bed()} function. The function reads in BED files, extract genomic regions, methylation, coverage, ambiguous modification status data and sample information and then construct Bsseq object using \code{BSseq} function within the package.}
}
\examples{
files <- c(system.file("extdata/modbam2bed/ctr1.chr10.chr11.bed.gz", package = "bsseq"),
    system.file("extdata/modbam2bed/ctr2.chr10.chr11.bed.gz", package = "bsseq"),
    system.file("extdata/modbam2bed/ctr3.chr10.chr11.bed.gz", package = "bsseq"),
    system.file("extdata/modbam2bed/tret1.chr10.chr11.bed.gz", package = "bsseq"),
    system.file("extdata/modbam2bed/tret2.chr10.chr11.bed.gz", package = "bsseq"),
    system.file("extdata/modbam2bed/tret3.chr10.chr11.bed.gz", package = "bsseq"))
pd <- data.frame(condition = rep(c("control", "treatment"), each = 3),
    replicate = rep(c("rep1", "rep2", "rep3"), times = 2))
bsseq_nano <- bsseq::read.modbam2bed(files,colData=pd,rmZeroCov = FALSE,
    strandCollapse=TRUE)
}
