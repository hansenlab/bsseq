\name{bsseqnews}
\title{bsseq news}
\encoding{UTF-8}

\section{Version 0.9.x}{
  \itemize{
    \item Fixed a problem with "width" in the title of bsseq plots.
    \item plot.BSseqTstat now allows for BSseqTstat objects computed
    without correction.
    \item validObject(BSseq) has been extended to also check for
    sampleNames consistency.
    \item Fixed a bug related to validity checking.
    \item Increased maxk from 10,000 to 50,000 in calls to locfit, to
    allowing fitting the model on genomes with unusally long chromosomes
    (Thanks to Brian Herb for reporting).
    \item The class representation for class 'BSseq' has changed
    completely. The new class is build on 'SummarizedExperiment' from
    GenomicRanges instead of 'hasGRanges'.  Use 'x <- updateObject(x)' to
    update serialized (saved) objects.
    \item Fixing a problem in orderBSseq related to chromosome names.
    \item Allowed user specification of maxk, with a default of 10,000
    in BSmooth.
    \item Many bugfixes made necessary by the new class representation.
  }
}

\section{Version 0.7.x}{
  \itemize{
    \item Removed the returnRaw argument to read.bismark() as it was
    unnecessary (Bismark output files does not have additional
    information beyond M and Cov and genomic positions, unlike BSmooth).
    \item Moved the Bismark example data from data to inst/extdata.
    \item combineList() now deals with the case where the list of BSseq
    objects have different genomic locations.  This speeds up
    read.bismark() substantially.
    \item Exposed combineList() as a faster alternative to
    Reduce(combine, list).
    \item Updated the code for the plotting routines (plotRegion).  This
    should not have an impact on user-visible code.
    \item Added read.bismark() function to parse output from the Bismark
    alignment suit [thanks to Pete Hickey].
    \item Refactorized plotting code.
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

