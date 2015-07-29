\name{read.bismark}
\alias{read.bismark}
\title{
  Parsing output from the Bismark alignment suite
}
\description{
  Parsing output from the Bismark alignment suite.
}
\usage{
  read.bismark(files,
               sampleNames,
               rmZeroCov = FALSE,
               strandCollapse = TRUE,
               fileType = c("coverage", "cytosineReport"),
               verbose = TRUE)
}
\arguments{
  \item{files}{Input files. Each sample is in a different file. Input
    files are created by running Bismark's
    \code{bismark_methylation_extractor}; see Note for details.}
  \item{sampleNames}{sample names, based on the order of \code{files}.}
  \item{rmZeroCov}{Should methylation loci that have zero coverage in
    all samples be removed. This will result in a much smaller object if
    the data originates from (targeted) capture bisulfite sequencing.}
  \item{strandCollapse}{Should strand-symmetric methylation loci, e.g., CpGs,
    be collapsed across strands.}
  \item{fileType}{The format of the input file; see Note for details.}
  \item{verbose}{Make the function verbose.}
}
\note{
  The format of the output of the Bismark alignment suite will depend on
  various user-specified options. Most relevant to \code{bsseq::read.bismkar()}
  is the options provided to \code{bismark_methylation_extractor}, part of the
  Bismark alignment suite.

  When calling \code{bsseq::read.bismark()}, the user specifies the relevant
  file format by the \code{fileType} argument. The information provided herein
  is accurate for Bismark \code{v0.14.3}; much of it is true for earlier
  versions of the software but may require that the user consults the relevant
  Bismark documentation.

  The "\code{coverage}" format has 6 columns (\code{chromosome},
  \code{position}, \code{strand}, \code{methylation percentage},
  \code{count methylated}, \code{count unmethylated}) and uses 1-based
  co-ordinates by default. The standard file extension for this file is
  \code{.cov}.

  % TODO: Give example of creating this output/input

  Please note that the bedGraph output created in the above call (file
  extension \code{.bedGraph}) is \strong{not} suitable for analysis with bsseq because
  it only contains the "\code{methylation percentage}" and not
  "\code{count methylated}" nor "\code{count unmethylated}".

  The "\code{cytosineReport}" format has 7 columns
  (\code{chromosome}, \code{position}, \code{strand}, \code{count methylated},
  \code{count unmethylated}, \code{C-context}, \code{trinucleotide context})
  and uses 1-based co-ordinates by default. There is no standard file extension
  for this file.

  % TODO: Give example of creating this output/input

  Input files can either be gzipped or not.
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
  infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov",
                        package = 'bsseq')
  bismarkBSseq <- read.bismark(files = infile,
                               sampleNames = "test_data",
                               rmZeroCov = FALSE,
                               strandCollapse = TRUE,
                               fileType = "coverage",
                               verbose = TRUE)
  bismarkBSseq
}

\author{
  Peter Hickey \email{peter.hickey@gmail.com}
}
