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

.BSmooth <- function(b, M, Cov, pos, coef_sink, se.coef_sink, sink_lock, grid,
                     pos_grid, ns, h, keep.se) {
    # Load packages on worker (required for SnowParam) -------------------------

    suppressPackageStartupMessages({
        requireNamespace("BiocParallel")
    })
    suppressPackageStartupMessages({
        requireNamespace("locfit")
    })
    suppressPackageStartupMessages({
        requireNamespace("DelayedArray")
    })

    # Construct inputs for smoothing -------------------------------------------

    # NOTE: 'bb' indexes over elements of pos_grid by cycling `ncol(M)` times
    #       over 1, ..., nrow(pos_grid).
    bb <- (b - 1L) %% nrow(pos_grid) + 1L
    # NOTE: unname() is necessary because M and Cov may carry colnames
    sdata <- data.frame(
        pos = read_block(pos, pos_grid[[bb]]),
        M = unname(read_block(M, grid[[b]])),
        Cov = unname(read_block(Cov, grid[[b]])))
    # Ensure 0 < M < Cov to avoid boundary issues (only relevant at loci with
    # non-zero coverage, so doesn't matter what M is for loci with zero
    # coverage).
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
    # Write to coef_sink and se.coef_sink while respecting the IPC lock.
    ipclock(sink_lock)
    write_block(coef_sink, viewport = grid[[b]], block = as.matrix(coef))
    if (keep.se) {
        write_block(
            se.coef_sink,
            viewport = grid[[b]],
            block = as.matrix(se.coef))
    }
    ipcunlock(sink_lock)
    NULL
}

# Exported functions -----------------------------------------------------------

# TODO: verbose = TRUE should report timings.
# TODO: BSmooth() should warn if BSseq object contains mix of strands.
# TODO: Consider having BSmooth() create a 'smoothed' assay in addition to or
#       instead of the 'coef' and 'se.coef' assays.
BSmooth <- function(BSseq,
                    ns = 70,
                    h = 1000,
                    maxGap = 10^8,
                    keep.se = FALSE,
                    BPPARAM = bpparam(),
                    chunkdim = NULL,
                    level = NULL,
                    verbose = getOption("verbose")) {

    # Argument checks-----------------------------------------------------------

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
    stopifnot(isTRUEorFALSE(keep.se))
    # Set appropriate BACKEND and check compatability with BPPARAM.
    BACKEND <- .getBSseqBackends(BSseq)
    if (!.areBackendsInMemory(BACKEND)) {
        if (!.isSingleMachineBackend(BPPARAM)) {
            stop("The parallelisation strategy must use a single machine ",
                 "when using an on-disk realization backend.\n",
                 "See help(\"BSmooth\") for details.",
                 call. = FALSE)
        }
    } else {
        if (!is.null(BACKEND)) {
            # NOTE: Currently do not support any in-memory realization
            #       backends. If 'BACKEND' is NULL then an ordinary matrix
            #       is returned rather than a matrix-backed DelayedMatrix.
            stop("The '", BACKEND, "' realization backend is not ",
                 "supported.\n",
                 "See help(\"BSmooth\") for details.",
                 call. = FALSE)
        }
    }
    # If using HDF5Array, check if BSseq object is updateable.
    if (identical(BACKEND, "HDF5Array")) {
        is_BSseq_updateable <- .isHDF5BackedBSseqUpdatable(BSseq)
        if (is_BSseq_updateable) {
            h5_path <- path(assay(BSseq, withDimnames = FALSE))
            if (any(c("coef", "se.coef") %in% rhdf5::h5ls(h5_path)[, "name"])) {
                # TODO: Better error message; what should be done in this case?
                stop("The HDF5 file '", h5_path, "' already contains a ",
                          "dataset named 'coef' or 'se.coef'.")
            }
        } else {
            h5_path <- NULL
            warning(
                wmsg("'BSseq' was not created with either read.bismark() or ",
                     "HDF5Array::saveHDF5SummarizedExperiment(). BSmooth() is ",
                     "using automatically created HDF5 file(s) (see ",
                     "?HDF5Array::setHDF5DumpFile) to store smoothing result. ",
                     "You will need to run ",
                     "HDF5Array::saveHDF5SummarizedExperiment() on the ",
                     "returned object if you wish to save the returned ",
                     "object."))
        }
    }
    stopifnot(isTRUEorFALSE(verbose))

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
    if (is.null(BACKEND)) {
        coef_sink <- NULL
        se.coef_sink <- NULL
        sink_lock <- NULL
    } else if (BACKEND == "HDF5Array") {
        coef_sink <- HDF5RealizationSink(
            dim = dim(M),
            dimnames = NULL,
            type = "double",
            filepath = h5_path,
            name = "coef",
            chunkdim = chunkdim,
            level = level)
        on.exit(close(coef_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
        if (keep.se) {
            se.coef_sink <- HDF5RealizationSink(
                dim = dim(M),
                dimnames = NULL,
                type = "double",
                filepath = h5_path,
                name = "se.coef",
                chunkdim = chunkdim,
                level = level)
            on.exit(close(se.coef_sink), add = TRUE)
        } else {
            se.coef_sink <- NULL
        }
    } else {
        # TODO: This branch should probably never be entered because we
        #       (implicitly) only support in-memory or HDF5Array backends.
        #       However, we retain it for now (e.g., fstArray backend would
        #       use this until a dedicated branch was implemented).
        coef_sink <- DelayedArray::AutoRealizationSink(dim(M), type = "double")
        on.exit(close(coef_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
        if (keep.se) {
            se.coef_sink <- DelayedArray::AutoRealizationSink(
                dim(M),
                type = "double")
            on.exit(close(se.coef_sink), add = TRUE)
        } else {
            se.coef_sink <- NULL
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
        sink_lock = sink_lock,
        grid = grid,
        pos_grid = pos_grid,
        ns = ns,
        h = h,
        keep.se = keep.se,
        BPPARAM = BPPARAM))
    if (!all(bpok(smooth))) {
        stop("BSmooth() encountered errors: ",
             sum(!bpok(smooth)), " of ", length(smooth),
             " smoothing tasks failed.")
    }
    # Construct coef and se.coef from results of smooth().
    if (is.null(BACKEND)) {
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

    # Construct BSseq object, saving it if it is HDF5-backed -------------------

    # NOTE: Using BiocGenerics:::replaceSlots(..., check = FALSE) to avoid
    #       triggering the validity method.
    assays <- c(assays(BSseq, withDimnames = FALSE)[c("M", "Cov")],
                SimpleList(coef = coef))
    if (keep.se) {
        assays <- c(assays, SimpleList(se.coef = se.coef))
    }
    BSseq <- BiocGenerics:::replaceSlots(
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
    if (identical(BACKEND, "HDF5Array") && is_BSseq_updateable) {
        # NOTE: Save BSseq object; mimicing
        #       HDF5Array::saveHDF5SummarizedExperiment().
        dir <- dirname(h5_path)
        x <- BSseq
        x@assays <- HDF5Array::shorten_assay2h5_links(x@assays)
        saveRDS(x, file = file.path(dir, "se.rds"))
    }
    BSseq
}

# TODOs ------------------------------------------------------------------------

# TODO: Use the logging facilities of BiocParallel. This is a longterm goal.
#       For example, we could set custom messages within .BSmooth() using the
#       futile.logger syntax; see the BiocParalell vignette 'Errors, Logs and
#       Debugging in BiocParallel'.
