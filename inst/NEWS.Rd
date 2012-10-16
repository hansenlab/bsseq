\name{bsseqnews}
\title{bsseq news}
\encoding{UTF-8}

\section{Version 0.7.x}{
  \itemize{
    \item Exposed combineList as a faster alternative to Reduce(combine,
    list). 
  }
}

\section{Version 0.6.x}{
  \itemize{
    \item Fixed a bug in getMeth, where type="raw" resulted in an error
    for non-smoothed data objects.
    \item Updated CITATION and citations in the vignettes.
    \item Now read.bsmooth supports both gzipped and non-gzipped files,
    whereas previously it assumed the output files to be gzipped.
    Thanks to Andreas Schoenegger for reporting this problem.
    \item Fixed a bug in combine() that also resulted in a bug in
    read.bsmooth when multiple input directories were specified.  Bug
    reported by Andreas Schoenegger.
  }
}

\section{Version 0.4}{
  \itemize{
    \item Improved combine and fixed a bug.  Also added a non-exported
    combineList for testing.
    \item Bug fix to read.bsmooth; it now works correctly for the
    default settings (= returning a single object of class BSseq and not a
    list). 
    \item Getting ready for initial release on Bioconductor.
    \item Updated the citations in the vignette(s) and the CITATION file.
  }
}

