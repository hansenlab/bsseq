# Functions/methods that would be good to have in DelayedArray -----------------

.rowVars <- function(x, rows = NULL, cols = NULL, ...) {
    if (is(x, "DelayedArray")) {
        if (!is.null(rows)) {
            x <- x[rows, ]
        }
        if (!is.null(cols)) {
            x <- x[, cols]
        }
        row_vars <- rowVars(as.array(x), ...)
    } else {
        row_vars <- rowVars(x, rows = rows, cols = cols, ...)
    }
    row_vars
}

.rowSds <- function(x, rows = NULL, cols = NULL, ...) {
    row_vars <- .rowVars(x, rows = rows, cols = cols, ...)
    sqrt(row_vars)
}

.quantile <- function(x, ...) {
    if (is(x, "DelayedArray")) {
        x <- as.array(x)
    }
    quantile(x, ...)
}

.isSimpleDelayedMatrix <- function(x) {
    is(x@seed, "matrix")
}

.zero_type <- function(type) {
    if (identical(type, "integer")) {
        fill <- 0L
    } else if (identical(type, "double")) {
        fill <- 0
    } else {
        stop("'type' = ", type, " is not supported")
    }
}

# A temporary workaround to
# https://github.com/Bioconductor/DelayedArray/issues/41.
.colsum <- function(x, group, reorder = TRUE, na.rm = FALSE, filepath = NULL,
                    name = NULL, chunkdim = NULL, level = NULL,
                    type = c("double", "integer"), BPPARAM = bpparam()) {

    # NOTE: Special case for HDF5Matrix, otherwise defer to rowsum().
    if (is(x, "HDF5Matrix")) {
        # Check arguments ------------------------------------------------------

        type <- match.arg(type)
        if (any(!c(type(x), type) %in% c("integer", "double"))) {
            stop("'type(x)' must be 'integer' or 'double'.")
        }
        if (length(group) != NCOL(x)) {
            stop("incorrect length for 'group'")
        }
        if (anyNA(group)) {
            warning("missing values for 'group'")
        }
        ugroup <- unique(group)
        if (reorder) {
            ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
        }
        # TODO: Default is type = "double" because rowSums2() returns numeric,
        #       but it can be useful to manually override this when you know
        #       the result is integer.

        # Construct RealizationSink --------------------------------------------

        # NOTE: This is ultimately coerced to the output DelayedMatrix
        #       object
        ans_nrow <- nrow(x)
        ans_ncol <- length(ugroup)
        ans_dim <- c(ans_nrow, ans_ncol)
        sink <- HDF5RealizationSink(
            dim = ans_dim,
            dimnames = list(rownames(x), as.character(ugroup)),
            type = type,
            filepath = filepath,
            name = name,
            chunkdim = chunkdim,
            level = level)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)

        # Construct ArrayGrid --------------------------------------------------

        sink_grid <- colAutoGrid(x = sink, ncol = 1L)
        list_of_cols <- split(seq_along(group), group)[ugroup]

        # Compute colsum() -----------------------------------------------------

        bplapply(
            X = seq_along(sink_grid),
            FUN = function(b, x, sink, sink_lock, sink_grid, list_of_cols) {
                cols <- list_of_cols[[b]]
                if (length(cols) == 1L) {
                    ans <- as.matrix(x[, cols, drop = FALSE])
                    if (na.rm) {
                        ans[is.na(ans)] <- 0L
                    }
                } else {
                    ans <- matrix(
                        rowSums2(x, cols = cols, na.rm = na.rm),
                        ncol = 1)
                }
                ipclock(sink_lock)
                write_block(sink, viewport = sink_grid[[b]], block = ans)
                ipcunlock(sink_lock)
                NULL
            },
            x = x,
            sink = sink,
            sink_lock = sink_lock,
            sink_grid = sink_grid,
            list_of_cols = list_of_cols,
            BPPARAM = BPPARAM)
        return(as(sink, "DelayedArray"))
    }

    colsum(x, group, reorder)
}

# A temporary workaround to
# https://github.com/Bioconductor/DelayedArray/issues/41.
.rowsum <- function(x, group, reorder = TRUE, na.rm = FALSE, filepath = NULL,
                    name = NULL, chunkdim = NULL, level = NULL,
                    type = c("double", "integer"), BPPARAM = bpparam()) {

    # NOTE: Special case for HDF5Matrix, otherwise defer to rowsum().
    if (is(x, "HDF5Matrix")) {

        # Check arguments ------------------------------------------------------

        if (any(!c(type(x), type) %in% c("integer", "double"))) {
            stop("'type(x)' must be 'integer' or 'double'.")
        }
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
        # NOTE: Default is type = "double" because colSums2() returns numeric,
        #       but it can be useful to manually override this when you know the
        #       result is integer.
        type <- match.arg(type)

        # Construct RealizationSink --------------------------------------------

        # NOTE: This is ultimately coerced to the output DelayedMatrix
        #       object
        ans_nrow <- length(ugroup)
        ans_ncol <- ncol(x)
        ans_dim <- c(ans_nrow, ans_ncol)
        sink <- HDF5RealizationSink(
            dim = ans_dim,
            dimnames = list(as.character(ugroup), colnames(x)),
            type = type,
            filepath = filepath,
            name = name,
            chunkdim = chunkdim,
            level = level)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)

        # Construct ArrayGrid --------------------------------------------------

        sink_grid <- rowAutoGrid(x = sink, nrow = 1L)
        list_of_rows <- split(seq_along(group), group)[as.character(ugroup)]

        # Compute colsum() -----------------------------------------------------

        bplapply(
            X = seq_along(sink_grid),
            FUN = function(b, x, sink, sink_lock, sink_grid, list_of_rows) {
                rows <- list_of_rows[[b]]
                if (length(rows) == 1L) {
                    ans <- as.matrix(x[rows, , drop = FALSE])
                    if (na.rm) {
                        ans[is.na(ans)] <- 0L
                    }
                } else {
                    ans <- matrix(
                        colSums2(x, rows = rows, na.rm = na.rm),
                        nrow = 1)
                }
                ipclock(sink_lock)
                write_block(sink, viewport = sink_grid[[b]], block = ans)
                ipcunlock(sink_lock)
                NULL
            },
            x = x,
            sink = sink,
            sink_lock = sink_lock,
            sink_grid = sink_grid,
            list_of_rows = list_of_rows,
            BPPARAM = BPPARAM)
        return(as(sink, "DelayedArray"))
    }

    rowsum(x, group, reorder)
}


# Missing methods --------------------------------------------------------------

# NOTE: Copied from minfi
# TODO: Perhaps move this to DelayedMatrixStats?
# TODO: DelayedArray::type() for all RealizationSink subclasses
setMethod("type", "HDF5RealizationSink", function(x) {
    x@type
})
# NOTE: Copied from minfi
# TODO: Perhaps move this to DelayedMatrixStats?
setMethod("type", "arrayRealizationSink", function(x) {
    DelayedArray::type(x@result_envir$result)
})
# NOTE: Copied from minfi
# TODO: Perhaps move this to DelayedMatrixStats?
setMethod("type", "RleRealizationSink", function(x) {
    x@type
})
# NOTE: Copied from minfi
# TODO: Perhaps move this to DelayedMatrixStats?
# TODO: dimnames() for all RealizationSink subclasses
setMethod("dimnames", "arrayRealizationSink", function(x) {
    dimnames(x@result_envir$result)
})

# Advanced block processing routines -------------------------------------------

# NOTE: Copy of minfi:::blockApplyWithRealization()
# TODO: Perhaps move this to DelayedMatrixStats?
# NOTE: DelayedArray::blockApply() with the option to write the blocks to
#       'sink'. Useful, for example, to apply a function across column-blocks
#       of a DelayedMatrix, write these results to disk, and then wrap
#       these in a DelayedMatrix.
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/10
blockApplyWithRealization <- function(x, FUN, ..., sink = NULL, x_grid = NULL,
                                      sink_grid = NULL, BPREDO = list(),
                                      BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)

    # Check conformable dots_grids and sinks_grids
    x_grid <- DelayedArray:::normarg_grid(x_grid, x)
    sink_grid <- DelayedArray:::normarg_grid(sink_grid, sink)
    if (!identical(dim(x_grid), dim(sink_grid))) {
        stop("non-conformable 'x_grid' and 'sink_grid'")
    }

    # Loop over blocks of `x` and write to `sink`
    nblock <- length(x_grid)
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        x_viewport <- x_grid[[b]]
        sink_viewport <- sink_grid[[b]]
        block <- read_block(x, x_viewport)
        attr(block, "from_grid") <- x_grid
        attr(block, "block_id") <- b
        block_ans <- FUN(block, ...)
        # NOTE: This is the only part different from DelayedArray::blockApply()
        if (!is.null(sink)) {
            write_block(sink, viewport = sink_viewport, block = block_ans)
            block_ans <- NULL
        }
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}

# TODO: Needed?
.getSEDir <- function(x) {
    paths <- lapply(assays(x, withDimnames = FALSE), function(a) {
        try(path(a), silent = TRUE)
    })
    if (any(vapply(paths, is, logical(1L), "try-error"))) {
        stop("Cannot extract 'dir'.")
    }
    unique_paths <- unique(unlist(paths, use.names = FALSE))
    if (length(unique_paths) > 1) {
        stop("Assay data spread across multiple HDF5 files.")
    }
    dirs <- dirname(unlist(paths, use.names = FALSE))
    unique(dirs)
}

# Should return TRUE for BSseq object created with read.bismark() or saved with
# HDF5Array::saveHDF5SummarizedExperiment().
# TODO: Check dirname(paths[[1L]]) also contains 'se.rds'? It looks like dir
#       can contain other files besides these; check.
.isHDF5BackedBSseqUpdatable <- function(x) {
    stopifnot(is(x, "BSseq"))
    assay_class <- vapply(assays(x, withDimnames = FALSE), class, character(1L))
    if (!all(assay_class == "HDF5Matrix")) {
        return(FALSE)
    }
    paths <- vapply(assays(x, withDimnames = FALSE), path, character(1L))
    if (!all(paths == paths[[1L]]) || !all(basename(paths) == "assays.h5")) {
        return(FALSE)
    }
    TRUE
}
