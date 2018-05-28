# Internal functions -----------------------------------------------------------

.rowTickmarks <- function(hasGRanges, maxGap) {
    gr <- granges(hasGRanges)
    # NOTE: This relies on 'gr' being sorted, so really want to be sure of this!
    stopifnot(!is.unsorted(gr))
    clusters <- reduce(gr, min.gapwidth = maxGap, with.revmap = TRUE)
    cumsum(lengths(mcols(clusters)[["revmap"]]))
}

# TODO: Remove other uses of this function, if possible.
# TODO: Rename to .makeClusters() since it's an internal function
makeClusters <- function(hasGRanges, maxGap = 10^8) {
    chrOrder <- as.character(runValue(seqnames(hasGRanges)))
    if(anyDuplicated(chrOrder))
        stop("argument 'hasGRanges' is not properly order")
    grBase <- granges(hasGRanges)
    clusters <- reduce(resize(grBase, width = 2*maxGap + 1, fix = "center"))
    start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
    clusters.sp <- split(clusters, seqnames(clusters))
    stopifnot(all(sapply(clusters.sp, function(cluster.gr) {
        if(length(cluster.gr) <= 1) return(TRUE)
        all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
    }))) # are the clusters ordered within the chromosome? This is probably guranteed
    clusters <- Reduce(c, clusters.sp[chrOrder])
    stopifnot(all(chrOrder == runValue(seqnames(clusters))))
    ov <- findOverlaps(grBase, clusters)
    clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
    names(clusterIdx) <- NULL
    clusterIdx
}

.BSmooth <- function(b, M, Cov, pos, coef_sink, se.coef_sink,
                     coef_sink_lock, se.coef_sink_lock, grid, pos_grid,
                     ns, h, keep.se) {
    # Load packages on worker (required for SnowParam) -------------------------

    suppressPackageStartupMessages({
        library(BiocParallel)
    })
    suppressPackageStartupMessages({
        library(locfit)
    })
    suppressPackageStartupMessages({
        library(DelayedArray)
    })

    # Construct inputs for smoothing -------------------------------------------

    # NOTE: 'bb' indexes over elements of pos_grid by cycling `ncol(M)` times
    #       over 1, ..., nrow(pos_grid).
    bb <- (b - 1L) %% nrow(pos_grid) + 1L
    # NOTE: unname() is necessary because M and Cov may carry colnames
    sdata <- data.frame(
        pos = DelayedArray:::extract_block(pos, pos_grid[[bb]]),
        M = unname(DelayedArray:::extract_block(M, grid[[b]])),
        Cov = unname(DelayedArray:::extract_block(Cov, grid[[b]])))
    # Ensure 0 < M < Cov to avoid boundary issues (only relevant at loci with
    # non-zero coverage, so doesn't matter what M is for these loci).
    sdata[["M"]] <- pmin(pmax(sdata[["M"]], 0.01), pmax(sdata[["Cov"]] - 0.01))
    n_loci <- nrow(sdata)
    n_loci_with_non_zero_coverage <- sum(sdata[["Cov"]] > 0)

    # Fit local binomial regression model --------------------------------------

    # NOTE: Require (ns + 1) loci with non-zero coverage to run smooth.
    if (n_loci_with_non_zero_coverage <= ns) {
        coef <- rep(NA_real_, n_loci)
        if (keep.se) {
            se.coef <- rep(NA_real_, n_loci)
        } else {
            se.coef <- NULL
        }
    } else {
        # Set 'nearest neighbour' smoothing parameter.
        nn <- ns / n_loci_with_non_zero_coverage
        # Fit model only using loci with non-zero coverage.
        fit <- locfit(
            M ~ lp(pos, nn = nn, h = h),
            data = sdata,
            weights = Cov,
            family = "binomial",
            subset = Cov > 0,
            maxk = 10000)
        # Evaluate smooth at all loci (regardless of coverage).
        pp <- preplot(
            object = fit,
            newdata = sdata["pos"],
            band = "local")
        coef <- pp$fit
        if (keep.se) {
            se.coef <- pp$se.fit
        } else {
            se.coef <- NULL
        }
    }

    # Return the results of the smooth or write them to the RealizationSink ----

    if (is.null(coef_sink)) {
        return(list(coef = coef, se.coef = se.coef))
    }
    # Write to coef_sink and se.coef_sink while respecting the IPC lock(s).
    ipclock(coef_sink_lock)
    write_block_to_sink(as.matrix(coef), coef_sink, grid[[b]])
    ipcunlock(coef_sink_lock)
    if (keep.se) {
        ipclock(se.coef_sink_lock)
        write_block_to_sink(as.matrix(se.coef), se.coef_sink, grid[[b]])
        ipcunlock(se.coef_sink_lock)
    }
}

# Exported functions -----------------------------------------------------------

# TODO: Make 'mc.cores', 'mc.preschedule', and 'verbose' defunct one release
#       cycle with them deprecated.
# TODO: Consider having BSmooth() create a 'smoothed' assay in addition to or
#       instead of the 'coef' and 'se.coef' assays.
# TODO: To benefit from error recovery requires that bpStopOnError(BPPARAM) is
#       TRUE but the default is FALSE. How to help the user? Probably don't
#       want to override the user-specified value.
BSmooth <- function(BSseq, ns = 70, h = 1000, maxGap = 10^8,
                    parallelBy = c("sample", "chromosome"),
                    mc.preschedule = FALSE, mc.cores = 1, keep.se = FALSE,
                    verbose = TRUE, BPREDO = list(), BPPARAM = bpparam()) {
    # Argument checks-----------------------------------------------------------

    # Check if this is a re-do.
    # NOTE: Under the assumptions of a re-do (i.e. BSmooth() is being re-run
    #       with the same arguments), we skip straight ahead to re-running
    #       failed smoothing tasks with no further argument checks.
    if (length(BPREDO)) {
        if (!is.list(BPREDO) ||
            identical(names(BPREDO), c("smooth", "coef_sink", "se.coef_sink",
                                       "realization_backend"))) {
            stop("'BPREDO' should be a list with elements 'smooth', ",
                 "'coef_sink', 'se.coef_sink', and 'realization_backend'.")
        }
        is_redo <- TRUE
        coef_sink <- BPREDO[["coef_sink"]]
        se.coef_sink <- BPREDO[["se.coef_sink"]]
        realization_backend <- BPREDO[["realization_sink"]]
        BPREDO <- BPREDO[["smooth"]]
    } else {
        is_redo <- FALSE
        # Check validity of 'BSseq' argument
        if (!is(BSseq, "BSseq")) {
            stop("'BSseq' must be a BSseq object.")
        }
        if (!isTRUE(all(width(BSseq) == 1L))) {
            stop("All loci in 'BSseq' must have width == 1.")
        }
        if (is.unsorted(BSseq)) {
            stop("'BSseq' must be sorted before smoothing. Use 'sort(BSseq)'.")
        }
        # Check for deprecated arguments and issue warning(s) if found.
        if (!missing(parallelBy)) {
            warning(
                "'parallelBy' is deprecated.\n",
                "See help(\"BSmooth\") for details.",
                call. = FALSE,
                immediate. = TRUE)
        }
        if (!missing(mc.preschedule)) {
            warning(
                "'mc.preschedule' is deprecated.\n",
                "See help(\"BSmooth\") for details.",
                call. = FALSE,
                immediate. = TRUE)
        }
        if (!missing(mc.cores)) {
            warning(
                "'mc.cores' is deprecated.\n",
                "See help(\"BSmooth\").",
                call. = FALSE,
                immediate. = TRUE)
            BPPARAM <- MulticoreParam(workers = mc.cores)
        }
        if (!missing(verbose)) {
            warning(
                "'verbose' is deprecated.\n",
                "See help(\"BSmooth\") for details.",
                call. = FALSE,
                immediate. = TRUE)
            if (verbose) bpprogressbar(BPPARAM) <- TRUE
        }
        # Check compatability of realization backend with backend(s) of BSseq
        # object.
        realization_backend <- getRealizationBackend()
        BSseq_backends <- .getBSseqBackends(BSseq)
        if (.areBackendsInMemory(realization_backend) &&
            !.areBackendsInMemory(BSseq_backends)) {
            stop("Using an in-memory backend for a disk-backed BSseq object ",
                 "is not supported.\n",
                 "See help(\"BSmooth\") for details.",
                 call. = FALSE)
        }
        # Check compatability of 'BPPARAM' with the realization backend.
        if (!.areBackendsInMemory(realization_backend)) {
            if (!.isSingleMachineBackend(BPPARAM)) {
                stop("The parallelisation strategy must use a single machine ",
                     "when using an on-disk realization backend.\n",
                     "See help(\"BSmooth\") for details.",
                     call. = FALSE)
            }
        } else {
            if (!is.null(realization_backend)) {
                # NOTE: Currently do not support any in-memory realization
                #       backends. If the realization backend is NULL then an
                #       ordinary matrix is returned rather than a matrix-backed
                #       DelayedMatrix.
                stop("The '", realization_backend, "' realization backend is ",
                     "not supported.\n",
                     "See help(\"BSmooth\") for details.",
                     call. = FALSE)
            }
        }
    }

    # Smoothing ----------------------------------------------------------------

    # Extract pieces of BSseq object required for smoothing
    M <- getCoverage(BSseq, type = "M", withDimnames = FALSE)
    Cov <- getCoverage(BSseq, type = "Cov", withDimnames = FALSE)
    pos <- matrix(start(BSseq), ncol = 1)
    # Set up ArrayGrid so that each block contains data for a single sample and
    # single cluster of loci.
    row_tickmarks <- .rowTickmarks(BSseq, maxGap)
    col_tickmarks <- seq_len(ncol(M))
    grid <- ArbitraryArrayGrid(list(row_tickmarks, col_tickmarks))
    # Set up "parallel" ArrayGrid over pos
    pos_grid <- ArbitraryArrayGrid(list(row_tickmarks, 1L))
    # Construct RealizationSink objects (as required)
    if (!is_redo) {
        if (is.null(realization_backend)) {
            coef_sink <- NULL
            coef_sink_lock <- NULL
            se.coef_sink <- NULL
            se.coef_sink_lock <- NULL
        } else {
            coef_sink <- DelayedArray:::RealizationSink(dim(M), type = "double")
            on.exit(close(coef_sink), add = TRUE)
            coef_sink_lock <- ipcid()
            on.exit(ipcremove(coef_sink_lock), add = TRUE)
            if (keep.se) {
                se.coef_sink <- DelayedArray:::RealizationSink(dim(M),
                                                               type = "double")
                on.exit(close(se.coef_sink), add = TRUE)
                se.coef_sink_lock <- ipcid()
                on.exit(ipcremove(se.coef_sink_lock), add = TRUE)
            } else {
                se.coef_sink <- NULL
                se.coef_sink_lock <- NULL
            }
        }
    }

    # Set number of tasks to ensure the progress bar gives frequent updates.
    # NOTE: The progress bar increments once per task
    #       (https://github.com/Bioconductor/BiocParallel/issues/54).
    #       Although it is somewhat of a bad idea to overrides a user-specified
    #       bptasks(BPPARAM), the value of bptasks(BPPARAM) doesn't affect
    #       performance in this instance, and so we opt for a useful progress
    #       bar. Only SnowParam (and MulticoreParam by inheritance) have a
    #       bptasks<-() method.
    if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
        bptasks(BPPARAM) <- length(grid)
    }
    smooth <- bptry(bplapply(
        X = seq_along(grid),
        FUN = .BSmooth,
        M = M,
        Cov = Cov,
        pos = pos,
        coef_sink = coef_sink,
        se.coef_sink = se.coef_sink,
        coef_sink_lock = coef_sink_lock,
        se.coef_sink_lock = se.coef_sink_lock,
        grid = grid,
        pos_grid = pos_grid,
        ns = ns,
        h = h,
        keep.se = keep.se,
        BPREDO = BPREDO,
        BPPARAM = BPPARAM))
    if (!all(bpok(smooth))) {
        # TODO: Feels like stop() rather than warning() should be used, but
        #       stop() doesn't allow for the return of partial results;
        #       see https://support.bioconductor.org/p/109374/
        warning("BSmooth() encountered errors: ",
                sum(!bpok(smooth)), " of ", length(smooth),
                " smoothing tasks failed.\n",
                "BSmooth() has returned partial results, including errors, ",
                "for debugging purposes.\n",
                "It may be possible to re-run just these failed smoothing ",
                "tasks.\nSee help(\"BSmooth\")",
                call. = FALSE)
        # NOTE: Return intermediate results as well as all derived variables.
        return(list(smooth = smooth,
                    coef_sink = coef_sink,
                    se.coef_sink = se.coef_sink,
                    realization_backend = realization_backend))
    }
    # Construct coef and se.coef from results of smooth().
    if (is.null(realization_backend)) {
        # Returning matrix objects.
        coef <- do.call(c, lapply(smooth, "[[", "coef"))
        attr(coef, "dim") <- dim(M)
        if (keep.se) {
            se.coef <- do.call(c, lapply(smooth, "[[", "se.coef"))
            attr(se.coef, "dim") <- dim(M)
        } else {
            se.coef <- NULL
        }
    } else {
        # Returning DelayedMatrix objects.
        coef <- as(coef_sink, "DelayedArray")
        if (keep.se) {
            se.coef <- as(se.coef_sink, "DelayedArray")
        } else {
            se.coef <- NULL
        }
    }

    # Construct BSseq object ---------------------------------------------------

    assays <- c(assays(BSseq, withDimnames = FALSE), SimpleList(coef = coef))
    if (keep.se) assays <- c(assays, SimpleList(se.coef = se.coef))
    BiocGenerics:::replaceSlots(
        object = BSseq,
        assays = Assays(assays),
        trans = plogis,
        parameters = list(
            smoothText = sprintf(
                "BSmooth (ns = %d, h = %d, maxGap = %d)", ns, h, maxGap),
            ns = ns,
            h = h,
            maxGap = maxGap),
        check = FALSE)
}

# TODOs ------------------------------------------------------------------------

# TODO: Use the logging facilities of BiocParallel. This is a longterm goal.
#       For example, we could set custom messages within .BSmooth() using the
#       futile.logger syntax; see the BiocParalell vignette 'Errors, Logs and
#       Debugging in BiocParallel'.
# TODO: Remove NOTEs that are really documentation issues to the docs
