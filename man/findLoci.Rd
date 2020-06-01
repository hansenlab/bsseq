\name{findLoci}
\alias{findLoci}
\title{
  Find methylation loci in a genome
}
\description{
    This is a convenience function to find methylation loci, such as CpGs, in a reference genome.
    The result is useful as the value of the \code{loci} argument for \code{\link{read.bismark}()}.
}
\usage{
findLoci(pattern,
         subject,
         include = seqlevels(subject),
         strand = c("*", "+", "-"),
         fixed = "subject",
         resize = TRUE)
}
\arguments{
  \item{pattern}{
    A string specifying the pattern to search for, e.g. \code{"CG"}.
    Can contain IUPAC ambiguity codes, e.g., \code{"CH"}.
  }
  \item{subject}{
    A string containing a file path to the genome sequence, in FASTA or 2bit format, to be searched.
    Alternatively, a \linkS4class{BSgenome} or \linkS4class{DNAStringSet} object storing the genome sequence to be searched.
  }
  \item{include}{
    A character vector indicating the seqlevels of \code{subject} to be used.
  }
  \item{strand}{
    A character scaler specifying the strand of \code{subject} to be used.
    If \code{strand = "*"}, then both the positive (\code{strand = "+"}) and negative (\code{strand = "-"} strands will be searched.)
    It is assumed that \code{subject} contains the sequence with respect to the \code{+} strand.
  }
  \item{fixed}{
    If \code{"subject"} (the default), IUPAC ambiguity codes in the pattern only are interpreted as wildcards, e.g., a \code{pattern} containing \code{CH} matches a \code{subject} containing \code{CA} but not vice versa.
    See \code{?Biostrings::`\link[Biostrings]{lowlevel-matching}`} for more information
  }
  \item{resize}{
    A logical scalar specifying whether the ranges should be shifted to have width 1 and anchored by the start of the locus, e.g., resize a CpG dinucleotide to give the co-ordinates of the cytosine.
  }
}
\details{
  This function provides a convenience wrapper for finding methylation loci in
  a genome, based on running
  \code{\link[Biostrings:matchPattern]{vmatchPattern}()}.
  Users requiring finer-grained control should directly use the
  \code{\link[Biostrings:matchPattern]{vmatchPattern}()} function and coerce
  the result to a \linkS4class{GRanges} object.
}
\value{
  A \linkS4class{GRanges} object storing the found loci.
}
\author{
  Peter F. Hickey
}
\seealso{
  \itemize{
    \item \code{Biostrings::\link[Biostrings:matchPattern]{vmatchPattern}()}
    \item \code{?BSgenome::`\link[BSgenome]{BSgenome-utils}`}
    }
}
\examples{
  library(Biostrings)
  # Find CpG dinucleotides on the both strands of an artificial sequence
  my_seq <- DNAStringSet("ATCAGTCGC")
  names(my_seq) <- "test"
  findLoci(pattern = "CG", subject = my_seq)
  # Find CHG trinucleotides on the both strands of an artificial sequence
  findLoci(pattern = "CHG", subject = my_seq)

  # Find CpG dinucleotides on the + strand of chr17 from the human genome (hg38)
  if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) {
    findLoci(
        pattern = "CG",
        subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        include = "chr17",
        strand = "+")
  }
}
