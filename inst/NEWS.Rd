\name{bsseqnews}
\title{bsseq news}
\encoding{UTF-8}

\section{Version 0.11.x}{
  \itemize{
    \item Converted to using Authors@R in the DESCRIPTIOn file.
    \item plotRegion did not respect the col/lty/lwd argument if given
    explicitely as opposed to through pData().  Reported by Karen
    Conneely <kconnee@emory.edu>.
    \item Fixed an issue introduced by the previous change (to
    plotRegion).  Reported (with a fix) by Tim Triche Jr
    <tim.triche@gmail.com>.
    \item Fixed a serious bug in collapseBSseq reported by Julien Roux
    <jroux@uchicago.edu>: it would use the Meth values instead of Cov
    for the final object's Cov.  However, this will result in the return
    object having a methylation of 100 percent across all loci, so
    hopefully users will have seen this for themselves.
    \item Fixed a bug in combineList which made combineList use "slow"
    code even in a special case where a faster shortcut was possible.
    This bug fix does not change the output of the function under any
    input, it only makes it faster.  Reported by Julien Roux
    <jroux@uchicago.edu>.
    \item validity now checks for the presence of colnames (sampleNames)
    which was assumed to be set.  Reported by Kevin Rue-Albrecht
    <kevin.rue@ucdconnect.ie>.
    \item Fixed a man page issue.
    \item Slightly changed CITATION.
    \item Added GitHub URL to DESCRIPTION.
  }
}

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
    \item Better argument checking in BSmooth.tstat.
    \item A few undocumented functions are now documented.
    \item Rewrote orderBSseq
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

