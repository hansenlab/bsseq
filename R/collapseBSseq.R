# Internal functions -----------------------------------------------------------

.collapseColData <- function(x, group, reorder = TRUE) {
    if (length(group) != NROW(x)) {
        stop("incorrect length for 'group'")
    }
    if (anyNA(group)) {
        warning("missing values for 'group'")
    }
    ugroup <- unique(group)
    if (reorder) {
        ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
    }
    idx <- split(seq_along(group), group)[ugroup]
    if (ncol(x)) {
        ans <- endoapply(x, function(xx) {
            sapply(X = idx, function(i) unique(xx[i]))
        })
        rownames(ans) <- names(idx)
        return(ans)
    }
    DataFrame(row.names = names(idx))
}

# Exported functions -----------------------------------------------------------

# TODO: Remove `replace`? It'd be a bad idea to use replace = TRUE if `dir` is
#       the same as the current dir for the BSseq object.
# NOTE: This is similar to edgeR::sumTechReps().
# TODO: Make optional the collapsing of colData?
# TODO: Document (and warn) that coef and se.coef aren't collapsed?
collapseBSseq <- function(BSseq, group, BPPARAM = bpparam(),
                          dir = tempfile("BSseq"), replace = FALSE,
                          chunkdim = NULL, level = NULL,
                          type = c("double", "integer")) {

    # Argument checks ----------------------------------------------------------

    if (!anyDuplicated(group)) {
        return(BSseq)
    }
    if (hasBeenSmoothed(BSseq)) {
        warning("Collapsing a smoothed BSseq object. You will need to ",
                "re-smooth using 'BSmooth()' on the returned object.")
    }
    if (length(group) != NCOL(BSseq)) {
        stop("incorrect length for 'group'")
    }
    if (anyNA(group)) {
        warning("missing values for 'group'")
    }
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
    # TODO: Additional argument checks (e.g., `replace`,
    #       HDF5Array:::.create_dir)

    # Collapse 'M' and 'Cov' matrices ------------------------------------------

    h5_path <- file.path(dir, "assays.h5")
    M <- colsum(
        x = getCoverage(BSseq, type = "M", withDimnames = FALSE),
        group = group,
        reorder = FALSE,
        BPPARAM = BPPARAM,
        filepath = h5_path,
        name = "M",
        chunkdim = chunkdim,
        level = level,
        type = type)
    Cov <- colsum(
        x = getCoverage(BSseq, type = "Cov", withDimnames = FALSE),
        group = group,
        reorder = FALSE,
        BPPARAM = BPPARAM,
        filepath = h5_path,
        name = "Cov",
        chunkdim = chunkdim,
        level = level,
        type = type)

    # Collapse 'colData' -------------------------------------------------------

    colData <- .collapseColData(
        x = colData(BSseq),
        group = group,
        reorder = FALSE)

    # Construct BSseq object, saving it if it is HDF5-backed -------------------

    se <- SummarizedExperiment(
        assays = SimpleList(M = unname(M), Cov = unname(Cov)),
        rowRanges = rowRanges(BSseq),
        colData = colData)
    # TODO: Is there a way to use the internal constructor with `check = FALSE`?
    #       Don't need to check M and Cov because this has already happened
    #       when files were parsed.
    # .BSseq(se, trans = function(x) NULL, parameters = list())
    bsseq <- new2("BSseq", se, check = FALSE)
    if (identical(BACKEND, "HDF5Array")) {
        # NOTE: Save BSseq object; mimicing
        #       HDF5Array::saveHDF5SummarizedExperiment().
        x <- bsseq
        x@assays <- HDF5Array:::.shorten_h5_paths(x@assays)
        saveRDS(x, file = file.path(dir, "se.rds"))
    }
    bsseq
}
