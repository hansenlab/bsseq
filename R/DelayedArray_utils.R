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

# Helper functions for setting up ArrayGrid instances --------------------------

# NOTE: Copy of minfi:::colGrid()
# TODO: Perhaps move this to DelayedMatrixStats?
colGrid <- function(x) {
    block_maxlen <- max(nrow(x),
                        DelayedArray:::get_default_block_maxlength(type(x)))
    spacings <- DelayedArray:::get_spacings_for_linear_capped_length_blocks(
        refdim = dim(x),
        block_maxlen = block_maxlen)
    RegularArrayGrid(dim(x), spacings)
}

# NOTE: Copy of minfi:::rowGrid()
# TODO: Perhaps move this to DelayedMatrixStats?
rowGrid <- function(x) {
    block_maxlen <- max(ncol(x),
                        DelayedArray:::get_default_block_maxlength(type(x)))
    spacings <- DelayedArray:::get_spacings_for_hypercube_capped_length_blocks(
        refdim = dim(x),
        block_maxlen = block_maxlen)
    RegularArrayGrid(dim(x), spacings)
}

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
        block <- DelayedArray:::extract_block(x, x_viewport)
        if (!is.array(block)) {
            block <- DelayedArray:::.as_array_or_matrix(block)
        }
        attr(block, "from_grid") <- x_grid
        attr(block, "block_id") <- b
        block_ans <- FUN(block, ...)
        # NOTE: This is the only part different from DelayedArray::blockApply()
        if (!is.null(sink)) {
            write_block_to_sink(block_ans, sink, sink_viewport)
            block_ans <- NULL
        }
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}

