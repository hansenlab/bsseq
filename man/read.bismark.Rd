\name{read.bismark}
\alias{read.bismark}
\title{
  Parsing output from the Bismark alignment suite.
}
\description{
  Parsing output from the Bismark alignment suite.
}
\usage{read.bismark(files,
             sampleNames,
             rmZeroCov = FALSE,
             strandCollapse = TRUE,
             fileType = c("cov", "oldBedGraph", "cytosineReport"),
             mc.cores = 1,
             verbose = TRUE,
             BACKEND = NULL)
}
\arguments{
  \item{files}{Input files. Each sample is in a different file. Input
    files are created by running Bismark's methylation extractor; see Note for
    details.}
  \item{sampleNames}{sample names, based on the order of \code{files}.}
  \item{rmZeroCov}{Should methylation loci that have zero coverage in
    all samples be removed. This will result in a much smaller object if
    the data originates from (targeted) capture bisulfite sequencing.}
  \item{strandCollapse}{Should strand-symmetric methylation loci, e.g., CpGs,
    be collapsed across strands. This option is only available if
    \code{fileType = "cytosineReport"} since the other file types do not
    contain the necessary strand information.}
  \item{fileType}{The format of the input file; see Note for details.}
  \item{mc.cores}{The number of cores used.  Note that setting
    \code{mc.cores} to a value greater than 1 is not  supported on MS
    Windows, see the help page for \code{mclapply}.}
  \item{verbose}{Make the function verbose.}
  \item{BACKEND}{Which backend should be used for the 'M' and 'Cov' matrices?
    The default, \code{NULL}, corresponds to using \link[base]{matrix}
    objects. See \code{?DelayedArray::\link[DelayedArray]{setRealizeBackend}}
    for alternative backends.}
}
\note{
  Input files can either be gzipped or not.

  The user must specify the relevant file format via the \code{fileType}
  argument. The format of the output of the Bismark alignment suite will depend
  on the version of Bismark and on various user-specified options. Please
  consult the Bismark documentation and the Bismark RELEASE NOTES
  (\url{http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/RELEASE_NOTES.txt})
  for the definitive list of changes between versions. When possible, it is
  strongly recommended that you use the most recent version of Bismark.

  The "\code{cov}" and "\code{oldBedGraph}" formats both have six columns
  ("\code{chromosome}", "\code{position}", "\code{strand}",
  "\code{methylation percentage}", "\code{count methylated}",
  "\code{count unmethylated}"). If you are using a recent version of Bismark
  (\code{v>=0.10.0}) then the standard file extension for this file is
  "\code{.cov}". If, however, you are using an older version of Bismark
  (\code{v<0.10.0}) then the file extension will be "\code{.bedGraph}". Please
  note that the "\code{.bedGraph}" file created in recent versions of Bismark
  (\code{v>=0.10.0}) is \strong{not} suitable for analysis with bsseq because
  it only contains the "\code{methylation percentage}" and not
  "\code{count methylated}" nor "\code{count unmethylated}".

  The "\code{cytosineReport}" format has seven columns
  ("\code{chromosome}", "\code{position}", "\code{strand}",
  "\code{count methylated}", "\code{count unmethylated}", "\code{C-context}",
  "\code{trinucleotide context}").
  There is no standard file extension for this file. The "\code{C-context}" and
  "\code{trinucleotide context}" columns are not currently used by bsseq.

  The following is a list of some issues to be aware of when using
  output from Bismark's methylation extractor:

  \itemize{
    \item The program to extract methylation counts was named
    \code{methylation_extractor} in older versions of Bismark (\code{v<0.8.0})
    and re-named \code{bismark_methylation_extractor} in recent versions of
    Bismark (\code{v>=0.8.0}). Furthermore, very old versions of Bismark
    (\code{v<0.7.7}) required that user run a separate script (called
    something like \code{genome_methylation_bismark2bedGraph}) to create the
    six-column "\code{cov}"/"\code{oldBedGraph}" file.
    \item The \code{--counts} and \code{--bedGraph} arguments must be supplied
    to \code{methylation_extractor}/\code{bismark_methylation_extractor} in
    order to use the output with \code{bsseq::read.bismark()}.
    \item The genomic co-ordinates of the Bismark output file may be zero-based
    or one-based depending on whether the \code{--zero_based} argument is used.
    Furthermore, the default co-ordinate system varies by version of Bismark.
    bsseq makes no assumptions about the basis of the genomic co-ordinates and
    it is left to the user to ensure that the appropriate basis is used in the
    analysis of their data. Since Bioconductor packages and
    \code{\linkS4class{GRanges}} use one-based co-ordinates, it
    is recommended that your Bismark files are also one-based.
  }
}
\value{
  An object of class \code{\link{BSseq}}.
}

\seealso{
  \code{\link{read.bsmooth}} for parsing output from the BSmooth
  alignment suite. \code{\link{read.umtab}} for parsing legacy (old)
  formats from the  BSmooth alignment suite.
  \code{\link{collapseBSseq}} for collapse (merging or summing) the
  data in two different directories.
}

\examples{
  infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
                        package = 'bsseq')
  bismarkBSseq <- read.bismark(files = infile,
                               sampleNames = "test_data",
                               rmZeroCov = FALSE,
                               strandCollapse = FALSE,
                               fileType = "cov",
                               verbose = TRUE)
  bismarkBSseq

  #-----------------------------------------------------------------------------
  # An example constructing a HDF5Array-backed BSseq object
  #
  library(HDF5Array)
  # See ?DelayedArray::setRealizeBackend for details
  hdf5_bismarkBSseq <- read.bismark(files = infile,
                                    sampleNames = "test_data",
                                    rmZeroCov = FALSE,
                                    strandCollapse = FALSE,
                                    fileType = "cov",
                                    verbose = TRUE,
                                    BACKEND = "HDF5Array")
}

\author{
  Peter Hickey \email{peter.hickey@gmail.com}
}
