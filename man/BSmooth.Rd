\name{BSmooth}
\alias{BSmooth}
\title{
  BSmooth, smoothing bisulfite sequence data
}
\description{
  This implements the BSmooth algorithm for estimating methylation levels from
  bisulfite sequencing data.}
\usage{
BSmooth(BSseq, ns = 70, h = 1000, maxGap = 10^8,
  parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
  mc.cores = 1, keep.se = FALSE, verbose = TRUE, BPREDO = list(),
  BPPARAM = bpparam())
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{ns}{The minimum number of methylation loci in a smoothing window.}
  \item{h}{The minimum smoothing window, in bases.}
  \item{maxGap}{The maximum gap between two methylation loci, before the
    smoothing is broken across the gap. The default smoothes each
    chromosome separately.}
  \item{parallelBy}{\strong{Deprecated}. See section, 'Parallelization and
    progress monitoring'.}
  \item{mc.preschedule}{\strong{Deprecated}. See section, 'Parallelization and
    progress monitoring'.}
  \item{mc.cores}{\strong{Deprecated}. See section, 'Parallelization and
    progress monitoring'.}
  \item{keep.se}{Should the estimated standard errors from the smoothing
    algorithm be kept. This will make the return object roughly 30
    percent bigger and is currently not be used for anything in \pkg{bsseq}.}
  \item{verbose}{\strong{Deprecated}. See section, 'Parallelization and progress
    monitoring'.}
  \item{BPREDO}{Typically not used. See section, 'Parallelization and progress
    monitoring'.}
  \item{BPPARAM}{An optional \linkS4class{BiocParallelParam} instance
    determining the parallel back-end to be used during evaluation. Currently
    supported are \linkS4class{SerialParam} (Unix, Mac, Windows),
    \linkS4class{MulticoreParam} (Unix and Mac), \linkS4class{SnowParam}
    (Unix, Mac, and Windows, limited to single-machine clusters), and
    \linkS4class{BatchJobsParam} (Unix, Mac, Windows, only with the in-memory
    realization backend). See sections 'Parallelization and progress
    monitoring' and 'Realization backends' for further details.}
}
\details{
  \code{ns} and \code{h} are passed to the \code{locfit} function. The
  bandwidth used is the maximum (in genomic distance) of the \code{h}
  and a width big enough to contain \code{ns} number of methylation
  loci.
}
\section{Realization backend}{
  The \code{BSmooth()} function creates a new assay to store the coefficients
  used to construct the smoothed methylation estimates ((\code{coef}). An
  additional assay is also created if \code{keep.se == TRUE} (\code{se.coef}).

  The choice of \emph{realization backend} controls whether these assay(s) are
  stored in-memory as an ordinary \link[base]{matrix} or on-disk as a
  \linkS4class{HDF5Array}, for example.

  The choice of realization backend is controlled by
  \code{\link[DelayedArray]{setRealizationBackend()}} from
  the \pkg{DelayedArray} package. The current realization backend can be
  queried by calling
  \code{DelayedArray::\link[DelayedArray]{getRealizationBackend()}}. Please see
  the examples at the bottom of this page.

  \code{BSmooth} supports the following realization backends:

  \itemize{
    \item \code{NULL} (in-memory): This stores each new assay in-memory using
      an ordinary \link[base]{matrix}.
    \item \code{HDF5Array} (on-disk): This stores each new assay on-disk in a
    HDF5 file using an \linkS4class{HDF5Matrix} from \pkg{HDF5Array}.
  }

  Please note that certain combinations of realization backend and
  parallelization backend are currently not supported. For example, the
  \linkS4class{HDF5Array} realization backend is currently only compatible when
  used with a single-machine parallelization backend (i.e. it is not compatible
  with a \linkS4class{SnowParam} that specifies an \emph{ad hoc} cluster of
  \strong{multiple} machines). \code{BSmooth()} will issue an error when given
  such incompatible realization and parallelization backends. Furthermore, to
  avoid memory usage blow-ups, \code{BSmooth()} will issue an error if an
  in-memory realization backend is used when smoothing a disk-backed
  \linkS4class{BSseq} object.
}

\section{Parallelization and progress monitoring}{
  \code{BSmooth} now uses the \pkg{BiocParallel} package to implement
  parallelization. This brings some notable improvements:

  \itemize{
    \item Smoothed results can now be written directly to an on-disk
      realization backend by the worker. This dramatically reduces memory
      usage compared to previous versions of \pkg{bsseq} that required all
      results be retained in-memory.
    \item Parallelization is now supported on Windows through the use of a
    \linkS4class{SnowParam} object as the value of \code{BPPARAM}.
    \item Improved error handling makes it possible to gracefully resume
    \code{BSmooth} jobs that encountered errors by re-doing only the necessary
    tasks.
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
    \code{parallelBy}, \code{mc.cores}, and \code{mc.preschedule} are
    deprecated and will be removed in subsequent releases of \pkg{bsseq}. These
    arguments were necessary when \code{BSmooth()} used the \pkg{parallel}
    package to implement parallelization, but this functionality is superseded
    by the aforementioned use of \pkg{BiocParallel}. We recommend that users
    who previously relied on these arguments switch to
    \code{BPPARAM = MulticoreParam(workers = mc.cores, progressbar = TRUE)}.
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

  \subsection{Re-running failed tasks}{
    By using \pkg{BiocParallel} for parallelization, \code{BSmooth()} may be
    able to resume smoothing runs that failed due to hardware issues, for
    example. To take advantage of this feature requires that
    \code{stop.on.error = FALSE} was used to construct the
    \linkS4class{BiocParallelParam} object passed to \code{BPPARAM}. When
    \code{BSmooth()} encounters a potentially recoverable error, it will return
    the intermediate results (including errors) and print instructions for how
    the job may be resumed. This typically means re-running the exact same
    command but with the addition of the the intermediate results passed via
    the \code{BPREDO} argument.

    Please note that not all errors can be recovered from in this way. We
    recommend inspecting the returned object of intermediate results to try to
    identify the source of the error (e.g., bad input data, failed node on a
    cluster job, memory-usage error on a shared-memory job). It may be helpful
    to use the advanced logging capabilities of \pkg{BiocParallel} when
    re-running the job.
  }
}

\value{
  An object of class \code{BSseq}, containing coefficients used to fit smoothed
  methylation values and optionally standard errors for these.
}
\author{
  Method and original implementation by Kasper Daniel Hansen
  \email{khansen@jhsph.edu}. Updated implementation to support disk-backed
  \linkS4class{BSseq} objects and more general parallelization by Peter Francis
  Hickey.
}
\references{
  KD Hansen, B Langmead, and RA Irizarry.
  \emph{BSmooth: from whole genome bisulfite sequencing reads to
    differentially methylated regions}.
  Genome Biology (2012) 13:R83.
  doi:\href{http://www.dx.doi.org/10.1186/gb-2012-13-10-r83}{10.1186/gb-2012-13-10-r83}.
}
\seealso{
  \code{\link[locfit]{locfit}} in the locfit package, as well as \code{\linkS4class{BSseq}}.
}
\examples{
\dontrun{
  # Run BSmooth() on a matrix-backed BSseq object using an in-memory realization
  # backend and use serial evaluation.
  data(BS.chr22)
  # This is a matrix-backed BSseq object
  sapply(assays(BS.chr22, withDimnames = FALSE), class)
  setRealizationBackend(NULL)
  BS.fit <- BSmooth(BS.chr22, BPPARAM = SerialParam(progressbar = TRUE))
  # The new 'coef' assay is an ordinary matrix
  sapply(assays(BS.fit, withDimnames = FALSE), class)
  BS.fit

  # Run BSmooth() on a disk-backed BSseq object using a disk-backed realization
  # backend and use multi-core parallel evaluation.
  BS.chr22 <- realize(BS.chr22, "HDF5Array")
  # This is a disk-backed BSseq object
  sapply(assays(BS.chr22, withDimnames = FALSE), class)
  setRealizationBackend("HDF5Array")
  BS.fit <- BSmooth(BS.chr22,
              BPPARAM = MulticoreParam(workers = 2, progressbar = TRUE))
  # The new 'coef' assay is an HDF5Matrix
  sapply(assays(BS.fit, withDimnames = FALSE), class)
  BS.fit
}
}
