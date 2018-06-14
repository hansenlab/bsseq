\name{read.bismark}
% TODO: alias required?
\alias{read.bismark}
\title{
  Parsing output from the Bismark alignment suite.
}
\description{
  Parsing output from the Bismark alignment suite.
}
\usage{read.bismark(files,
             sampleNames = basename(files),
             loci = NULL,
             rmZeroCov = FALSE,
             strandCollapse = TRUE,
             verbose = TRUE,
             BPPARAM = bpparam(),
             BACKEND = getRealizationBackend(),
             ...,
             fileType = c("cov", "oldBedGraph", "cytosineReport"),
             mc.cores = 1)
}
\arguments{
  \item{files}{The path to the files created by running Bismark's methylation
    extractor, one sample per file.
    Files ending in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip} will be automatically uncompressed.
    We strongly recommend you use the 'genome wide cytosine report' output files.
    See section 'File formats' for further details.}
  \item{sampleNames}{A character vector of sample names, based on the order of \code{files}.}
  \item{loci}{\code{NULL} (default) or a \linkS4class{GenomicRanges} instance
    containing methylation loci (all with width equal to 1).
    If \code{loci = NULL}, then \code{read.bismark()} will perform a first pass over the Bismark file to identify candidate loci.
    If \code{loci} is a \linkS4class{GenomicRanges} instance, then these form the candidate loci.
    In either case, the candidate loci will be filtered if
    \code{rmZeroCov = TRUE} and collapsed if \code{strandCollapse = TRUE} to form the final set of methylation loci that form the \code{\link[SummarizedExperiment]{rowRanges}} of the returned \linkS4class{BSseq} object.
    See section 'Efficient use of \code{read.bismark()}' for further details.}
  \item{rmZeroCov}{A \code{logical(1)} indicating whether methylation loci that
    have zero coverage in all samples be removed.
    For several reasons, the default \code{rmZeroCov = FALSE} is recommended even in cases where you ultimately want to remove such loci.
    See section 'Efficient use of \code{read.bismark()}' for further details.}
  \item{strandCollapse}{A \code{logical(1)} indicating whether strand-symmetric
    methylation loci (i.e. CpGs) should be collapsed across strands.
    This is only applicable for stranded methylation loci, e.g., loci extracted
    from 'genome wide cytosine reports' (see section 'File formats' for further details).}
  \item{verbose}{A \code{logical(1)} indicating whether progress messages
    should be printed (default \code{TRUE}).}
\item{BPPARAM}{An optional \linkS4class{BiocParallelParam} instance
    determining the parallel back-end to be used during evaluation.
    Currently supported are \linkS4class{SerialParam} (Unix, Mac, Windows),
    \linkS4class{MulticoreParam} (Unix and Mac), \linkS4class{SnowParam}
    (Unix, Mac, and Windows, limited to single-machine clusters), and
    \linkS4class{BatchJobsParam} (Unix, Mac, Windows, only with the in-memory
    realization backend).
    See sections 'Parallelization and progress  monitoring' and 'Realization backends' for further details.}
  \item{BACKEND}{\code{NULL} or a single string specifying the name of the
    realization backend.
    When the backend is set to \code{NULL}, the \code{M} and \code{Cov} assays
    are realized in memory as ordinary matrices, otherwise these are realized with the given \code{BACKEND}.
    See section 'Realization backends' for further details.}
  \item{...}{Additional arguments passed to the \linkS4class{RealizationSink}
    constructor.
    See section 'Realization backends' for further details.}
  \item{fileType}{\strong{Deprecated}.
    The format of each file is now automatically detected using the internal function \code{bsseq:::.guessBismarkFileType()}.}
  \item{mc.cores}{\strong{Deprecated}.
    See section 'Parallelization and progress monitoring' for further details.}
}

\section{Efficient use of \code{read.bismark()}}{
  We recommend the following to achieve fast and efficient importing of Bismark files:

  \itemize{
    \item Specify the set of methylation loci via the \code{loci} argument.
    \item Leave \code{rmZeroCov = FALSE}.
    \item Use a \code{BPPARAM} with a large number of workers (cores).
    \item Use \code{BACKEND = "HDF5Array"}.
  }

  Each point is discussed below.
  \subsection{Specifying \code{loci}}{
    Specifying the set of methylation loci via the \code{loci} argument means that \code{read.bismark()} does not need to first parse all \code{files} to identify the set of candidate loci.
    Provided that \code{rmZeroCov = FALSE}, this means that each file is only read once.
    This may be a considerable saving when there are a large number of \code{files}.

    \strong{If you are unsure whether the below-described shortcuts apply to your data, leave \code{loci = NULL} and let \code{read.bismark()} identify the set of candidate loci from \code{files}.}

    If all \code{files} are 'genome wide cytosine reports' for samples aligned to the same reference genome, then all \code{files} contain the exact same set of methylation loci.
    In this case, you may wish to first construct \code{loci} using the internal function \code{bsseq:::.readBismarkAsFWGRanges()} applied to a single file, e.g., \code{loci = bsseq:::.readBismarkAsFWGRanges(files[1], rmZeroCov, strandCollapse)}.

    % TODO: Uncomment once findCytosines is exported and documented.
    % Alternatively, you may wish to use the \code{\link{findCytosines}()} function to construct all methylation loci of interest (e.g., all CpGs) in your reference genome and then pass this via the \code{loci} argument.
  }

  \subsection{Leaving \code{rmZeroCov = FALSE}}{
    If you set \code{rmZeroCov = TRUE}, then \code{read.bismark()} must first parse all the \code{files} to identify which loci have zero coverage in all samples and then filter these out from the set of candidate loci.
    \strong{This will happen regardless of whether \code{loci} is \code{NULL} or is supplied as a \linkS4class{GenomicRanges} of candidate loci.}

    Furthermore, any coverage-based filtering of methylation loci is best left until you have constructed your final \linkS4class{BSseq} object.
    In our experience, the final \linkS4class{BSseq} object is often the product of combining multiple \linkS4class{BSseq} objects, each constructed with a separate call to \code{read.bismark()}.
    In such cases, it is premature to use \code{rmZeroCov = TRUE} when running each \code{read.bismark()} and combining these objects will regretably often lead to an inefficiently constructed final \linkS4class{BSseq} object.
  }

  \subsection{Using a \code{BPPARAM} with a large number of workers (cores)}{
    Each file can be processed on its own, so you can process in parallel as many files as you have workers. We have good experience using as many as 40 cores via \code{BPPARAM = MulticoreParam(workers = 40)}.}

  \subsection{Using \code{BACKEND = "HDF5Array"}}{
    By using the HDF5Array realization backend from \pkg{HDF5Array}, we reduce the amount of data that is kept in-memory at any one time.
    Once each file is parsed, the data are written to the HDF5 file and are no longer needed in-memory.
    When combined with multiple workers (cores), this means that each file will only need to read and retain in-memory 1 sample's worth of data at a time.

    Conversely, if you opt for all data to be in-memory (via \code{BACKEND = NULL}), then each worker will pass each file's data back to the main process and memory usage will steadily accumulate to often unreasonable levels.
    }
}

\section{File formats}{
  The format of each file is automatically detected using the internal function \code{bsseq:::.guessBismarkFileType()}.
  Files ending in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip} will be automatically uncompressed.
  \subsection{Supported file formats}{
     We strongly recommend that you use the 'genome wide cytosine report' output format
     (\url{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-genome-wide-cytosine-report-optional-is-tab-delimited-in-the-following-format-1-based-coords}) as this includes trand information for each locus.
     The 'coverage' output format (\url{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-coverage-output-looks-like-this-tab-delimited-1-based-genomic-coords}) is also supported, although all loci from these files will have their \code{\link[BiocGenerics]{strand}} set to \code{*} in the returned \linkS4class{BSseq} object because this file format does not include strand information.
  }

  \subsection{Unsupported file formats}{
     Neither the 'bedGraph' output format (\url{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-bedgraph-output-optional-looks-like-this-tab-delimited-0-based-start-1-based-end-coords}) nor the 'bismark_methylation_extractor' output format (\url{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-bismark_methylation_extractor-output-is-in-the-form-tab-delimited-1-based-coords}) are supported.
     The former does not include the required counts of methylated and unmethylated reads hile the is an intermediate file containing read-level, rather than locus-level, data on methylation.
  }
  \subsection{One-based vs. zero-based genomic co-ordinates}{
    The genomic co-ordinates of the Bismark output files may be zero-based or one-based depending on whether the \code{--zero_based} argument was used when running Bismark's methylation extractor.
    Furthermore, the default co-ordinate counting system varies by version of Bismark.
    \pkg{bsseq} makes no assumptions about the basis of the genomic co-ordinates and it is left to the user to ensure that the appropriate basis is used in the analysis of their data.

    Since Bioconductor packages typically use one-based co-ordinates, we strongly recommend that your Bismark files are also one-based.
  }
}

\section{Realization backends}{
  The \code{read.bismark()} function creates a \linkS4class{BSseq} object with two assays, \code{M} and \code{Cov}.
  The choice of \emph{realization backend} controls whether these assays are stored in-memory as an ordinary \link[base]{matrix} or on-disk as a \linkS4class{HDF5Array}, for example.
  The choice of realization backend is controlled by the \code{BACKEND} argument, which defaults to the current value of \code{DelayedArray::\link[DelayedArray]{getRealizationBackend}()}.

  \code{read.bismark()} supports the following realization backends:

  \itemize{
    \item \code{NULL} (in-memory): This stores each new assay in-memory using an ordinary \link[base]{matrix}.
    \item \code{HDF5Array} (on-disk): This stores each new assay on-disk in a HDF5 file using an \linkS4class{HDF5Matrix} from \pkg{HDF5Array}.
  }

  Please note that certain combinations of realization backend and parallelization backend are currently not supported.
  For example, the \linkS4class{HDF5Array} realization backend is currently only compatible when used with a single-machine parallelization backend (i.e. it is not compatible
  with a \linkS4class{SnowParam} that specifies an \emph{ad hoc} cluster of
  \strong{multiple} machines).
  % TODO: Test this
  \code{BSmooth()} will issue an error when given such incompatible realization and parallelization backends.

  Additional arguments related to the realization backend can be passed via the \code{...} argument.
  These arguments must be named and are passed to the relevant \linkS4class{RealizationSink} constructor.
  For example, the \code{...} argument can be used to specify the path to the HDF5 file to be
  used by \code{BSmooth()}.
  Please see the examples at the bottom of the page.
}

\section{Parallelization and progress monitoring}{
  \code{read.bismark()} now uses the \pkg{BiocParallel} package to implement
  parallelization. This brings some notable improvements:

  \itemize{
    \item Imported files can now be written directly to an on-disk
      realization backend by the worker. This dramatically reduces memory
      usage compared to previous versions of \pkg{bsseq} that required all
      results be retained in-memory.
    \item Parallelization is now supported on Windows through the use of a
    \linkS4class{SnowParam} object as the value of \code{BPPARAM}.
    % TODO: Uncomment when/if implemented.
    %\item Improved error handling makes it possible to gracefully resume
    %\code{read.bismark()} jobs that encountered errors by re-doing only the necessary
    %tasks.
    \item Detailed and extensive job logging facilities.
  }

  All parallelization options are controlled via the \code{BPPARAM} argument.
  In general, we recommend that users combine multicore (single-machine)
  parallelization with an on-disk realization backend (see section,
  'Realization backend'). For Unix and Mac users, this means using
  a \linkS4class{MulticoreParam}. For Windows users, this means using a
  single-machine \linkS4class{SnowParam}. Please consult the \pkg{BiocParallel}
  documentation to take full advantage of the more advanced features.

  \subsection{Deprecated arguments}{
    \code{mc.cores} is deprecated and will be removed in subsequent releases of \pkg{bsseq}.
    This argument was necessary when \code{read.bismark()} used the \pkg{parallel} package to implement parallelization, but this functionality is superseded by the aforementioned use of \pkg{BiocParallel}.
    We recommend that users who previously relied on these arguments switch to \code{BPPARAM = MulticoreParam(workers = mc.cores, progressbar = TRUE)}.
  }

  \subsection{Progress monitoring}{
    A useful feature of \pkg{BiocParallel} are progress bars to monitor the
    status of long-running jobs, such as \code{BSmooth()}. Progress bars are
    controlled via the \code{progressbar} argument in the
    \linkS4class{BiocParallelParam} constructor. Progress bars replace the
    use of the deprecated \code{verbose} argument to print out information on
    the status of \code{BSmooth()}.

    \pkg{BiocParallel} also supports extensive and detailed logging facilities.
    Please consult the \pkg{BiocParallel} documentation to take full advantage
    these advanced features.
  }
}

\value{
  A \linkS4class{BSseq} object.
}

\seealso{
  \code{\link{read.bsmooth}()} for parsing output from the BSmooth
  alignment suite.
  \code{\link{read.umtab}()} for parsing legacy (old) formats from the BSmooth alignment suite.
  \code{\link{collapseBSseq}()} for collapsing (aggregating) data from sample's with multiple Bismark methylation extractor files (e.g., technical replicates).
}

\examples{
  # Run read.bismark() on a single sample to construct a matrix-backed BSseq
  # object.
  infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
                        package = 'bsseq')
  bsseq <- read.bismark(files = infile,
                        sampleNames = "test_data",
                        rmZeroCov = FALSE,
                        strandCollapse = FALSE,
                        verbose = TRUE)
  # This is a matrix-backed BSseq object.
  sapply(assays(bsseq, withDimnames = FALSE), class)
  bsseq

  \dontrun{
  # Run read.bismark() on a single sample to construct a HDF5Array-backed BSseq
  # object (with data written to the file 'read.bismark_example.h5')
  bsseq <- read.bismark(files = infile,
                        sampleNames = "test_data",
                        rmZeroCov = FALSE,
                        strandCollapse = FALSE,
                        verbose = TRUE,
                        BACKEND = "HDF5Array",
                        filepath = "read.bismark_example.h5")
  # This is a HDF5Array-backed BSseq object.
  sapply(assays(bsseq, withDimnames = FALSE), class)
  # The 'M' and 'Cov' assays are in the HDF5 file 'read.bismark_example.h5' (in
  # the current working directory).
  sapply(assays(bsseq, withDimnames = FALSE), path)
  }
}

\author{
  Peter Hickey \email{peter.hickey@gmail.com}
}
