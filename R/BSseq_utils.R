chrSelectBSseq <- function(BSseq, seqnames = NULL, order = FALSE) {
    seqlevels(BSseq, pruning.mode = "coarse") <- seqnames
    if (order) BSseq <- orderBSseq(BSseq, seqOrder = seqnames)
    BSseq
}

orderBSseq <- function(BSseq, seqOrder = NULL) {
    if (!is.null(seqOrder)) {
        seqlevels(BSseq, pruning.mode = "coarse") <- seqOrder
    }
    BSseq[order(granges(BSseq))]
}

# TODO: getMeth() realises the result in memory iff regions is not NULL;
#       discuss with Kasper
# TODO: Whether or not colnames are added to returned value depends on whether
#       regions is non-NULL; discuss with Kasper
# TODO: Add parallel support
# TODO: Document withDimnames
# TODO: Should is be an explicit withDimnames arg or simply passed via `...`?
getMeth <- function(BSseq, regions = NULL, type = c("smooth", "raw"),
                    what = c("perBase", "perRegion"), confint = FALSE,
                    alpha = 0.95, withDimnames = TRUE) {
    p.conf <- function(p, n, alpha) {
        z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
        upper <- (p + z ^ 2 / (2 * n) +
                      z * sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
            (1 + z ^ 2 / n)
        lower <- (p + z ^ 2 / (2 * n) -
                      z * sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
            (1 + z ^ 2 / n)
        return(list(meth = p, lower = lower, upper = upper))
    }

    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    if (type == "smooth" & !hasBeenSmoothed(BSseq)) {
        stop("'type=smooth' requires the object to have been smoothed.")
    }
    what <- match.arg(what)
    if (what == "perRegion" & is.null(regions)) {
        stop("'what=perRegion' but no 'regions' supplied")
    }
    z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
    if (is.null(regions) && type == "smooth") {
        coef <- getBSseq(BSseq, "coef", withDimnames)
        meth <- getBSseq(BSseq, "trans", withDimnames)(coef)
        if (confint) {
            upper <- meth + z * getBSseq(BSseq, "se.coef", withDimnames)
            lower <- meth - z * getBSseq(BSseq, "se.coef", withDimnames)
            return(list(meth = meth, lower = lower, upper = upper))
        } else {
            return(meth)
        }
    }
    if (is.null(regions) && type == "raw") {
        meth <- getBSseq(BSseq, "M", withDimnames) /
            getBSseq(BSseq, "Cov", withDimnames)
        if (confint) {
            return(p.conf(meth, getBSseq(BSseq, "Cov", withDimnames), alpha))
        } else {
            return(meth)
        }
    }

    ## At this point, regions have been specified
    if (is(regions, "data.frame")) {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))
    if (confint) {
        stop("'confint = TRUE' is not supported by 'getMeth' when regions is given")
    }
    grBSseq <- granges(BSseq)
    ov <- findOverlaps(grBSseq, regions)
    # NOTE: This realises a large object in memory (`meth`) - could do it in
    #       chunks if what = perRegion
    if (type == "smooth") {
        meth <- as.matrix(
            getBSseq(BSseq, "trans", withDimnames)(
                getBSseq(BSseq, "coef", withDimnames))[
                    queryHits(ov), , drop = FALSE])
    } else if (type == "raw") {
        meth <- as.matrix(
            (getBSseq(BSseq, "M", withDimnames) /
                 getBSseq(BSseq, "Cov", withDimnames))[
                     queryHits(ov), , drop = FALSE])
    }
    out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
    if (what == "perBase") {
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the results end up in the wrong order wrt to regions?
        outList <- vector("list", length(regions))
        outList[as.integer(names(out))] <- out
        return(outList)
    } else if (what == "perRegion") {
        out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the rows end up in the wrong order?
        outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
        if (withDimnames) colnames(outMatrix) <- sampleNames(BSseq)
        outMatrix[as.integer(rownames(out)), ] <- out
        outMatrix
    }
}

# TODO: getCoverage() realises the result in memory iff regions is not NULL;
#       discuss with Kasper
# TODO: Whether or not colnames are added to returned value depends on whether
#       regions is non-NULL; discuss with Kasper
# TODO: Document withDimnames
# TODO: Should is be an explicit withDimnames arg or simply passed via `...`?
getCoverage <- function(BSseq, regions = NULL, type = c("Cov", "M"),
                        what = c("perBase", "perRegionAverage",
                                 "perRegionTotal"),
                        withDimnames = TRUE) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    what <- match.arg(what)
    if (is.null(regions)) {
        if (what == "perBase") {
            return(getBSseq(BSseq, type, withDimnames))
        }
        if (what == "perRegionTotal") {
            return(colSums2(getBSseq(BSseq, type, withDimnames)))
        }
        if (what == "perRegionAverage") {
            return(colMeans2(getBSseq(BSseq, type, withDimnames)))
        }
    }
    if (is(regions,  "data.frame")) {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))
    grBSseq <- granges(BSseq)
    ov <- findOverlaps(grBSseq, regions)
    coverage <- getBSseq(BSseq, type, withDimnames)[
        queryHits(ov), , drop = FALSE]
    out <- lapply(split(coverage, subjectHits(ov)), matrix,
                  ncol = ncol(coverage))

    if (what == "perBase") {
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the results end up in the wrong order wrt to regions?
        outList <- vector("list", length(regions))
        outList[as.integer(names(out))] <- out
        return(outList)
    } else if (what == "perRegionAverage") {
        out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
    } else if (what == "perRegionTotal") {
        out <- do.call(rbind, lapply(out, colSums2, na.rm = TRUE))
    }
    # TODO: Don't really understand the logic of the remaining code; how
    #       could the rows end up in the wrong order?
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    if (withDimnames) colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)), ] <- out
    outMatrix
}

strandCollapse <- function(x, shift = TRUE, BPPARAM = bpparam(),
                           BACKEND = getAutoRealizationBackend(),
                           dir = tempfile("BSseq"), replace = FALSE,
                           chunkdim = NULL, level = NULL,
                           type = c("double", "integer")) {

    # Argument checks ----------------------------------------------------------

    if (all(runValue(strand(x)) == "*")) {
        warning("All loci are unstranded, nothing to collapse.", call. = FALSE)
        return(BSseq)
    }
    if (!(all(runValue(strand(x)) %in% c("+", "-")))) {
        stop("'x' object has a mix of stranded and unstranded loci.")
    }
    # Register 'BACKEND' and return to current value on exit.
    # TODO: Is this strictly necessary?
    current_BACKEND <- getAutoRealizationBackend()
    on.exit(setAutoRealizationBackend(current_BACKEND), add = TRUE)
    setAutoRealizationBackend(BACKEND)
    # Check compatability of 'BPPARAM' with 'BACKEND'.
    if (!.areBackendsInMemory(BACKEND)) {
        if (!.isSingleMachineBackend(BPPARAM)) {
            stop("The parallelisation strategy must use a single machine ",
                 "when using an on-disk realization backend.\n",
                 "See help(\"read.bismark\") for details.",
                 call. = FALSE)
        }
    } else {
        if (!is.null(BACKEND)) {
            # NOTE: Currently do not support any in-memory realization
            #       backends. If the realization backend is NULL then an
            #       ordinary matrix is returned rather than a matrix-backed
            #       DelayedMatrix.
            stop("The '", BACKEND, "' realization backend is not supported.",
                 "\n  See help(\"read.bismark\") for details.",
                 call. = FALSE)
        }
    }
    # If using HDF5Array as BACKEND, check remaining options are sensible.
    if (identical(BACKEND, "HDF5Array")) {
        # NOTE: Most of this copied from
        #       HDF5Array::saveHDF5SummarizedExperiment().
        if (!isSingleString(dir)) {
            stop(wmsg("'dir' must be a single string specifying the path to ",
                      "the directory where to save the BSseq object (the ",
                      "directory will be created)."))
        }
        if (!isTRUEorFALSE(replace)) {
            stop("'replace' must be TRUE or FALSE")
        }
        if (!dir.exists(dir)) {
            HDF5Array::create_dir(dir)
        } else {
            HDF5Array::replace_dir(dir, replace)
        }
        h5_path <- file.path(dir, "assays.h5")
    } else if (identical(BACKEND, NULL)) {
        h5_path <- NULL
    }

    # Collapse loci ------------------------------------------------------------

    loci <- rowRanges(x)
    if (shift) {
        loci <- shift(loci, shift = as.integer(-1L * (strand(loci) == "-")))
    }
    collapsed_loci <- reduce(loci, min.gapwidth = 0L, ignore.strand = TRUE)

    # Collapse 'M' and 'Cov' matrices ------------------------------------------

    ol <- findOverlaps(loci, collapsed_loci, type = "equal")
    group <- subjectHits(ol)
    M <- .rowsum(
        x = assay(x, "M", withDimnames = FALSE),
        group = group,
        # NOTE: reorder = TRUE to ensure same row-order as collapsed_loci.
        reorder = TRUE,
        BPPARAM = BPPARAM,
        filepath = h5_path,
        name = "Cov",
        chunkdim = chunkdim,
        level = level,
        type = type)
    if(is(x, "BSseq")) {
        Cov <- .rowsum(
            x = assay(x, "Cov", withDimnames = FALSE),
            group = group,
            ## NOTE: reorder = TRUE to ensure same row-order as collapsed_loci.
            reorder = TRUE,
            BPPARAM = BPPARAM,
            filepath = h5_path,
            name = "Cov",
            chunkdim = chunkdim,
            level = level,
            type = type)
        ## Construct BSseq object, saving it if it is HDF5-backed -------------------

        se <- SummarizedExperiment(
            assays = SimpleList(M = unname(M), Cov = unname(Cov)),
            rowRanges = collapsed_loci,
            colData = colData(x))
        ## TODO: Is there a way to use the internal constructor with `check = FALSE`?
        ##       Assuming input was valid, the output is valid, too.
        ## .BSseq(se, trans = function(x) NULL, parameters = list())
        out <- new2("BSseq", se, check = FALSE)
    }
    if(is(x, "MethylCounts")) {
        U <- .rowsum(
            x = assay(x, "U", withDimnames = FALSE),
            group = group,
            ## NOTE: reorder = TRUE to ensure same row-order as collapsed_loci.
            reorder = TRUE,
            BPPARAM = BPPARAM,
            filepath = h5_path,
            name = "U",
            chunkdim = chunkdim,
            level = level,
            type = type)
        if("H" %in% assayNames(x)) {
            H <- .rowsum(
                x = assay(x, "H", withDimnames = FALSE),
                group = group,
                ## NOTE: reorder = TRUE to ensure same row-order as collapsed_loci.
                reorder = TRUE,
                BPPARAM = BPPARAM,
                filepath = h5_path,
                name = "H",
                chunkdim = chunkdim,
                level = level,
                type = type)
        } else {
            H <- NULL
        }
        if("D" %in% assayNames(x)) {
            D <- .rowsum(
                x = assay(x, "D", withDimnames = FALSE),
                group = group,
                ## NOTE: reorder = TRUE to ensure same row-order as collapsed_loci.
                reorder = TRUE,
                BPPARAM = BPPARAM,
                filepath = h5_path,
                name = "D",
                chunkdim = chunkdim,
                level = level,
                type = type)
        } else {
            D <- NULL
        }
        se <- SummarizedExperiment(
            assays = SimpleListExcludeNULL(M = unname(M), U = unname(U), H = unname(H), D = unname(D)),
            rowRanges = collapsed_loci,
            colData = colData(x))
        ## TODO: Is there a way to use the internal constructor with `check = FALSE`?
        ##       Assuming input was valid, the output is valid, too.
        ## .BSseq(se, trans = function(x) NULL, parameters = list())
        out <- new2("MethylCounts", se, check = FALSE)
    }

    if (!is.null(BACKEND) && BACKEND == "HDF5Array") {
        # NOTE: Save BSseq object; mimicing
        #       HDF5Array::saveHDF5SummarizedExperiment().
        xtmp <- out
        xtmp@assays <- HDF5Array::shorten_assay2h5_links(xtmp@assays)
        base::saveRDS(xtmp, file = file.path(dir, "se.rds"))
    }
    out
}

