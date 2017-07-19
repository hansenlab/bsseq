\name{BS.chr22}
\alias{BS.chr22}
\docType{data}
\title{
  Whole-genome bisulfite sequencing for chromosome 22 from Lister et al.
}
\description{
  This dataset represents chromosome 22 from the IMR90 cell line
  sequenced in Lister et al.  Only CpG methylation are included (there
  were very few non-CpG loci).  The two samples are two different
  extractions from the same cell line (ie. technical replicates), and
  are pooled in the analysis in the original paper.
}
\usage{data(BS.chr22)}
\format{
  An object of class \code{BSseq}.
}
\details{
  All coordinates are in hg18.
}
\source{
  Obtained from \url{http://neomorph.salk.edu/human_methylome/data.html}
  specifically the files \url{mc_imr90_r1.tar.gz} and
  \url{mc_imr90_r2.tar.gz}.  A script which downloads these files and
  constructs the \code{BS.chr22} object may be found in
  \file{inst/scripts/get_BS.chr22.R}, see the example.
}
\references{
  R Lister et al.
  \emph{Human DNA methylomes at base resolution show widespread
    epigenomic differences}.
  Nature (2009) 462, 315-322.
}
\examples{
data(BS.chr22)
BS.chr22

script <- system.file("scripts", "get_BS.chr22.R", package = "bsseq")
script
readLines(script)
}
\keyword{datasets}
