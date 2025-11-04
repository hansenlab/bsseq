\name{read.bedMethyl}
\alias{read.bedMethyl}
\title{
  Parsing bedMethyl output from modkit pileup.
}
\description{
  Parsing bedMethyl output from modkit pileup.
}
\usage{
read.bedMethyl(files,
               loci = NULL,
               colData = NULL,
               rmZeroCov = TRUE,
               strandCollapse = TRUE,
               BPPARAM = bpparam(),
               BACKEND = NULL,
               dir = tempfile("BSseq"),
               replace = FALSE,
               chunkdim = NULL,
               level = NULL,
               nThread = 1L,
               output = c("BSseq", "MethylCounts"),
               verbose = getOption("verbose"))
}
\arguments{
  \item{files}{The path to the files created by running modkit pileup, one sample per file. See the methods section of [link to preprint] for validated output.}
  \item{loci}{\code{NULL} (default) or a \code{\link[GenomicRanges]{GenomicRanges}} instance containing methylation loci (all with width equal to 1). If \code{loci = NULL}, then \code{read.bedMethyl()} will perform a first pass over the bedMethyl files to identify candidate loci. If \code{loci} is a \code{\link[GenomicRanges]{GenomicRanges}} instance, then these form the candidate loci. The candidate loci will be collapsed if \code{strandCollapse = TRUE}.}
  \item{colData}{An optional \code{\link[S4Vectors]{DataFrame}} describing the samples. Row names, if present, become the column names of the \linkS4class{BSseq} object. If \code{NULL}, then a \code{\link[S4Vectors]{DataFrame}} will be created with \code{files} used as the row names.}
  \item{rmZeroCov}{A \code{logical(1)} indicating whether methylation loci that have zero coverage in all samples should be removed. Default setting is \code{rmZeroCov = TRUE}}
  \item{strandCollapse}{A \code{logical(1)} indicating whether strand-symmetric methylation loci (i.e. CpGs) should be collapsed across strands.}
  \item{BPPARAM}{An optional \code{\link[BiocParallel]{BiocParallelParam}} instance determining the parallel back-end to be used during evaluation.}
  \item{BACKEND}{\code{NULL} or a single string specifying the name of the realization backend. Currently, the backend is not supported for downstream applications.}
  \item{dir}{\strong{Only applicable if \code{BACKEND == "HDF5Array"}.} The path (as a single string) to the directory where to save the HDF5-based \linkS4class{BSseq} object.}
  \item{replace}{\strong{Only applicable if \code{BACKEND == "HDF5Array"}.} If the directory \code{dir} already exists, should it be replaced with a new one?}
  \item{chunkdim}{\strong{Only applicable if \code{BACKEND == "HDF5Array"}.} The dimensions of the chunks to use for writing the data to disk.}
  \item{level}{The compression level to use for writing the data to disk.}
  \item{nThread}{The number of threads used by \code{\link[data.table]{fread}} when reading the \code{files}.}
  \item{output}{The type of output returned when reading the \code{files}. Allowed values are "BSseq" and "MethylCounts")}
  \item{verbose}{A \code{logical(1)} indicating whether progress messages should be printed (default \code{TRUE}).}
}

\section{File formats}{
  The format of each file should be similar to the examples in [link to preprint]. Files ending in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip} will be automatically decompressed to \code{\link{tempdir}()}.
  \subsection{Supported file formats}{
    Modkit bedMethyl files from modkit pileup. For downstream likelihood functions we recommend running modkit pileup on output from bam files modification/basecalled using a CG context model and not using a reference genome for pileup.
  }
  \subsection{Unsupported file formats}{
    Other types of output.
  }
  \subsection{One-based vs. zero-based genomic co-ordinates}{
    The genomic co-ordinates of bedMethyl files are zero-based. Since Bioconductor packages typically use one-based co-ordinates, the co-ordinates from the bedMethyl files are converted to one-based in the BSseq object.
  }
}

\examples{
# Example: Reading bedMethyl files included in the bsseq package
# Paths to example bedMethyl files in the package's extdata directory
infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
                         package = "bsseq"),
             system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
                         package = "bsseq"))

# Run the function to import data
bsseq <- read.bedMethyl(files = infiles,
                        colData = DataFrame(row.names = c("test_nanopore",
                                                          "test_pacbio")),
                        rmZeroCov = FALSE,
                        strandCollapse = TRUE,
                        verbose = TRUE)

# View the resulting BSseq object
bsseq
}

\author{
  Søren Blikdal Hansen (soren.blikdal.hansen@sund.ku.dk)
}
