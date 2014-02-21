\name{read.umtab}
\alias{read.umtab}
\alias{read.umtab2}
\title{
  Parsing UM tab files (legacy output) containing output from the
  BSmooth aligner. 
}
\description{
  Parsing UM tab files containing output from the bisulfite aligner
  \code{Merman}.  This is two different legacy formats, which we keep
  around.  These functions are likely to be deprecated in the future.
}
\usage{
read.umtab(dirs, sampleNames = NULL, rmZeroCov = FALSE,
  pattern = NULL, keepU = c("U10", "U20", "U30", "U40"),
  keepM = c("M10", "M20", "M30", "M40"), verbose = TRUE)

read.umtab2(dirs, sampleNames = NULL, rmZeroCov = FALSE,
  readCycle = FALSE, keepFilt = FALSE,
  pattern = NULL, keepU, keepM, verbose = TRUE)
}
\arguments{
  \item{dirs}{Input directories.  Usually each sample is in a different
    directory.} 
  \item{pattern}{An optional pattern, see \code{list.files} in the base
    package.} 
  \item{sampleNames}{sample names, based on the order of \code{dirs}.}
  \item{rmZeroCov}{Should methylation loci that have zero coverage in
    all samples be removed.  This will result in a much smaller object
    if the data originates from (targeted) capture bisulfite
    sequencing.} 
  \item{keepU}{A vector of U columns which are kept.}
  \item{keepM}{A vector of M columns which are kept.}
  \item{readCycle}{Should the cycle columns be returned?}
  \item{keepFilt}{Should the filter columns be returned?}
  \item{verbose}{Make the function verbose.}
}
\value{
  Both functions returns lists, the components are
  \item{BSdata}{An object of class \code{BSseq} containing the
    methylation evidence.}
  \item{GC}{A vector of local GC content values.}
  \item{Map}{A vector of local mapability values.}
  \item{Mcy}{A matrix of the number of unique M cycles.}
  \item{Ucy}{A matrix of the number of unique U cycles.}
  \item{chr}{A vector of chromosome names.}
  \item{pos}{A vector of genomic positions.}
  \item{M}{A matrix representing methylation evidence.}
  \item{U}{A matrix representing un-methylation evidence.}
  \item{csums}{Description of 'comp2'}
}
\details{
  \code{read.umtab2} is newer than \code{read.umtab} and both process
  output from older versions of the BSmooth alignment suite (versions
  older than 0.6.1).  These functions are likely to be deprecated in the
  future.  Newer output from the BSmooth alignment suite can be parsed
  using \code{read.bsmooth}.

  A script using this function can be found in the \code{bsseqData}
  package, in the file \file{scripts/create_BS.cancer.R}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
  \code{\link{read.bsmooth}}.
}
\examples{
\dontrun{
require(bsseqData)
umDir <- system.file("umtab", package = "bsseqData")
sampleNames <- list.files(umDir)
dirs <- file.path(umDir, sampleNames, "umtab")
umData <- read.umtab(dirs, sampleNames)
}
}
