\name{read.bismark}
\alias{read.bismark}
\title{
  Parsing output from the Bismark alignment suite
}
\description{
  Parsing output from the Bismark alignment suite.
}
\usage{
  read.bismark(files, sampleNames, returnRaw = FALSE, rmZeroCov = FALSE, verbose = TRUE)
}
\arguments{
  \item{files}{Input files. Each sample is in a different file. Input files are created by running Bismark's \code{methylation_extractor}; see Note for details.} 
  \item{sampleNames}{sample names, based on the order of \code{files}.}
  \item{returnRaw}{Should the function return the complete information in the output files?}
  \item{rmZeroCov}{Should methylation loci that have zero coverage in all samples be removed. This will result in a much smaller object if the data originates from (targeted) capture bisulfite sequencing.} 
  \item{verbose}{Make the function verbose.}
}
\note{
  Input files can either be gzipped or not.
  
  Input files are created by running Bismark's \code{methylation_extractor} and \code{genome_methylation_bismark2bedGraph_v4.pl} scripts over the Bismark alignment file. For example (run from the command line):
  
  \code{methylation_extractor -s --comprehensive test_data.fastq_bismark.sam}
  
  \code{genome_methylation_bismark2bedGraph_v4.pl --counts CpG_context_test_data.fastq_bismark.txt > CpG_context_test_data.fastq_bismark.bedGraph}
  
  The \code{--comprehensive} argument to \code{methylation_extractor} and the \code{--counts} argument to \code{genome_methylation_bismark2bedGraph_v4.pl} are required.
  
  In this example, the file \code{CpG_context_test_data.fastq_bismark.bedGraph} is then the input file to \code{read.bismark}. 
  
  See \url{http://rpubs.com/PeteHaitch/readBismark} for a worked example using Bismark and \code{read.bismark}. Please consult the Bismark website for full details of these scripts and the latest versions (\url{http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark})
}
\value{
  Either an object of class \code{"BSseq"} (if \code{returnRaw = FALSE})
  or a list of \code{"GRanges"} which each component coming from a
  file. 
}

\seealso{
  \code{\link{read.bsmooth}} for parsing output from the BSmooth alignment suite. \code{\link{read.umtab}} for parsing legacy (old) formats from the
  BSmooth alignment suite.  \code{\link{collapseBSseq}} for collapse
  (merging or summing) the data in two different directories.
}

\examples{
  \dontrun{
  bismarkBedGraph <- system.file('data/CpG_context_test_data.fastq_bismark.bedGraph', package = 'bsseq')
  bismarkBSseq <- read.bismark(files = bismarkBedGraph, sampleNames = 'test_data', returnRaw = FALSE, rmZeroCov = FALSE, verbose = TRUE)
  bismarkBSseq
  }
}

\author{
  Peter Hickey \email{peter.hickey@gmail.com}

}
