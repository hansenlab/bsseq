\name{read.bsmooth}
\alias{read.bsmooth}
\title{
  Parsing output from the BSmooth alignment suite
}
\description{
  Parsing output from the BSmooth alignment suite.
}
\usage{
read.bsmooth(dirs, sampleNames = NULL, seqnames = NULL,
  returnRaw = FALSE, qualityCutoff = 20, rmZeroCov = FALSE,
  verbose = TRUE)
}
\arguments{
  \item{dirs}{Input directories.  Usually each sample is in a different
    directory, or perhaps each (sample, lane) is a different directory.} 
  \item{sampleNames}{sample names, based on the order of \code{dirs}.
    If \code{NULL} either set to \code{basename(dirs)} (if unique) or
    \code{dirs}.}
  \item{seqnames}{The default is to read all BSmooth output files in
    \code{dirs}.  Using this argument, it is possible to restrict this
    to only files with names in \code{seqnames} (apart from
    \code{.cpg.tsv} and optionally \code{.gz}).}
  \item{returnRaw}{Should the function return the complete information
    in the output files?}
  \item{qualityCutoff}{Only use evidence (methylated and unmethylated
    evidence) for a given methylation loci, if the base in the read has
    a quality greater than this cutoff.}
  \item{rmZeroCov}{Should methylation loci that have zero coverage in
    all samples be removed.  This will result in a much smaller object
    if the data originates from (targeted) capture bisulfite
    sequencing.} 
  \item{verbose}{Make the function verbose.}
}
\note{
  Input files can either be gzipped or not.  Gzipping the input files
  results in much greater speed of reading (and saves space), so it is
  recommended.
  
  We are working on making this function faster and less memory hungry.
}
\value{
  Either an object of class \code{BSseq} (if \code{returnRaw = FALSE})
  or a list of \code{GRanges} which each component coming from a
  directory. 
}
\seealso{
  \code{\link{read.umtab}} for parsing legacy (old) formats from the
  BSmooth alignment suite.  \code{\link{collapseBSseq}} for collapse
  (merging or summing) the data in two different directories.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
