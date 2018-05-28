# Internal generics ------------------------------------------------------------

# 'idx': A list of columns (MARGIN = 1) or rows (MARGIN = 2) to collapse over.
# E.g., `idx = list(1:3, 4:5)` with `MARGIN = 1` means collapse columns 1-3
# into a single columna dn columns 4-5 into a single column.
# MARGIN = 1: collapse along columns using rowSums2
# MARGIN = 2: collapse along rows using colSums2
# `...` are additional arguments passed to methods.
setGeneric(
    ".collapseMatrixLike",
    function(x, idx, MARGIN, ...) standardGeneric,
    signature = "x")

# Internal methods -------------------------------------------------------------

setMethod(".collapseMatrixLike", "matrix", function(x, idx, MARGIN) {
    if (MARGIN == 1L) {
        # TODO: Need colnames?
        return(do.call(cbind, lapply(idx, function(j) rowSums2(x, cols = j))))
    } else if (MARGIN == 2L) {
        # TODO: Need rownames?
        do.call(rbind, lapply(idx, function(i) colSums2(x, rows = i)))
    } else {
        stop("'MARGIN' must be 1 or 2")
    }
})

setMethod(
    ".collapseMatrixLike",
    "DelayedMatrix",
    function(x, idx, MARGIN, BPREDO = list(), BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink objects of appropriate
        # dimensions and type
        # NOTE: `type = "double"` because .collapseMatrixLike,matrix-method
        #       uses colSums2()/rowSums2(), which returns a numeric vector.
        # NOTE: This is ultimately coerced to the output DelayedMatrix
        #       object
        # Set up ArrayGrid instances over `x` as well as "parallel"
        # ArrayGrid instances over `sink`.
        if (MARGIN == 1L) {
            sink <- DelayedArray:::RealizationSink(
                dim = c(nrow(x), length(idx)),
                dimnames = list(rownames(x), names(idx)),
                type = "double")
            x_grid <- minfi:::colGrid(x)
            sink_grid <- RegularArrayGrid(
                refdim = dim(sink),
                spacings = c(nrow(sink), ncol(sink) / length(x_grid)))
        } else if (MARGIN == 2L) {
            # TODO: Check sink has correct dim and dimnames
            sink <- DelayedArray:::RealizationSink(
                dim = c(length(idx), ncol(x)),
                dimnames = list(names(idx), colnames(x)),
                type = "double")
            on.exit(close(sink))
            x_grid <- minfi:::rowGrid(x)
            sink_grid <- RegularArrayGrid(
                refdim = dim(sink),
                spacings = c(nrow(sink) / length(x_grid), ncol(sink)))
        } else {
            stop("'MARGIN' must be 1 or 2")
        }

        # Loop over blocks of 'x' and write to 'sink'.
        minfi:::blockApplyWithRealization(
            x = x,
            FUN = .collapseMatrixLike,
            idx = idx,
            MARGIN = MARGIN,
            sink = sink,
            x_grid = x_grid,
            sink_grid = sink_grid,
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix object
        as(sink, "DelayedArray")
    }
)

# Exported functions -----------------------------------------------------------

collapseBSseq <- function(BSseq, columns) {
    if (hasBeenSmoothed(BSseq)) {
        warning("Collapsing a smoothed BSseq object. You will need to ",
                "re-smooth using 'BSmooth()' on the returned object.")
    }

    # Construct index between current samples and collapsed samples
    stopifnot(is.character(columns))
    if (is.null(names(columns)) && length(columns) != ncol(BSseq)) {
        stop("if `columns' does not have names, it needs to be of the same ",
             "length as `BSseq` has columns (samples)")
    }
    if (!is.null(names(columns)) &&
        !all(names(columns) %in% sampleNames(BSseq))) {
        stop("if `columns` has names, they need to be sampleNames(BSseq)")
    }
    if (is.null(names(columns))) {
        columns.idx <- seq_len(ncol(BSseq))
    } else {
        columns.idx <- match(names(columns), sampleNames(BSseq))
    }
    idx <- split(columns.idx, columns)

    # Collapse 'M' and 'Cov' matrices
    M <- .collapseMatrixLike(
        x = assay(BSseq, "M", withDimnames = FALSE),
        idx = idx,
        MARGIN = 1L)
    Cov <- .collapseMatrixLike(
        x = assay(BSseq, "Cov", withDimnames = FALSE),
        idx = idx,
        MARGIN = 1L)

    # Construct BSseq object
    # TODO: Check sampleNames are preserved (could extract from names(idx) and
    #       pass down to constructors).
    se <- SummarizedExperiment(
        assays = SimpleList(M = M, Cov = Cov),
        rowRanges = rowRanges(BSseq))
    .BSseq(se)
}
