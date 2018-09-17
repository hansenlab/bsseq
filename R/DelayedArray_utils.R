# Functions/methods that would be good to have in DelayedArray

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

.DelayedMatrix <- function(x) {
    x_name <- deparse(substitute(x))
    X <- try(DelayedArray(x), silent = TRUE)
    if (is(X, "try-error")) {
        stop("Could not construct DelayedMatrix from '", x_name, "'",
             call. = FALSE)
    }
    if (!is(X, "DelayedMatrix")) {
        stop("'", x_name, "' must be matrix-like", call. = FALSE)
    }
    X
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
    x_grid <- DelayedArray:::.normarg_grid(x_grid, x)
    sink_grid <- DelayedArray:::.normarg_grid(sink_grid, sink)
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
            write_block(x = sink, viewport = sink_viewport, block = block_ans)
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
    if (!identical(.getBSseqBackends(x), "HDF5Array")) {
        return(FALSE)
    }
    paths <- vapply(assays(x, withDimnames = FALSE), path, character(1L))
    if (all(paths == paths[[1L]]) && all(basename(paths) == "assays.h5")) {
        return(TRUE)
    }
    FALSE
}
