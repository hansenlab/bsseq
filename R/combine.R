# Internal functions -----------------------------------------------------------

# TODO: Option to specify that objects should be combined into DelayedArray
#       with a given backend (e.g., even if all inputs are matrix)?
.combineList <- function(x, x_rowRanges, ans_rowRanges) {
    x_is_matrix <- vapply(x, is.matrix, logical(1L))
    if (isTRUE(all(x_is_matrix))) {
        return(.combineList_matrix(x, x_rowRanges, ans_rowRanges))
    }
    .combineList_matrix_like(x, x_rowRanges, ans_rowRanges)
}

.combineList_matrix <- function(x, x_rowRanges, ans_rowRanges) {
    stopifnot(identical(length(x), length(x_rowRanges)))
    # Set up output matrix with appropriate dimension and type
    # TODO: Modify minfi:::.highestType() to accept a list. Can then do
    #       minfi:::.highestType(x)
    x_types <- vapply(x, DelayedArray::type, character(1L))
    ans_type <- typeof(do.call("c", lapply(x_types, vector)))
    ans <- matrix(
        data = .zero_type(ans_type),
        nrow = length(ans_rowRanges),
        ncol = sum(vapply(x, ncol, integer(1L))))

    # Fill output matrix
    col_offset <- 0L
    for (k in seq_along(x_rowRanges)) {
        ol <- findOverlaps(x_rowRanges[[k]], ans_rowRanges, type = "equal")
        # row_idx <- findOverlaps(
        #     query = x_rowRanges[[k]],
        #     subject = ans_rowRanges,
        #     type = "equal",
        #     select = "first")
        col_idx <- col_offset + seq_len(ncol(x[[k]]))
        ans[subjectHits(ol), col_idx] <- x[[k]][queryHits(ol), , drop = FALSE]
        col_offset <- col_offset + ncol(x[[k]])
    }
    # Return output matrix
    ans
}

.combineList_matrix_like <- function(x, x_rowRanges, ans_rowRanges) {
    # Argument checks
    stopifnot(identical(length(x), length(x_rowRanges)))

    # Construct data frame mapping each sample to an element of x and its
    # corresponding column number.
    x_ncol <- vapply(x, ncol, integer(1L))
    x_samples_df <- data.frame(
        sample = seq(1, sum(x_ncol)),
        x_idx = rep(seq_along(x), times = x_ncol),
        column_idx = unlist(lapply(x_ncol, function(end) seq(1, end))))

    # Set up intermediate RealizationSink object of appropriate dimensions
    # and type
    # TODO: Modify minfi:::.highestType() to accept a list. Can then do
    #       minfi:::.highestType(x)
    x_types <- vapply(x, DelayedArray::type, character(1L))
    ans_type <- typeof(do.call("c", lapply(x_types, vector)))
    sink <- DelayedArray:::RealizationSink(
        dim = c(length(ans_rowRanges), sum(vapply(x, ncol, integer(1L)))),
        type = ans_type)
    on.exit(close(sink))

    # Set up ArrayGrid instance over columns of `sink`.
    sink_grid <- minfi:::colGrid(sink)

    # Loop over column grid of 'sink', identify samples required for that
    # block, bring that data into memory, pass down to .combineList_matrix(),
    # and write result to sink.
    for (k in seq_along(sink_grid)) {
        sink_viewport <- sink_grid[[k]]
        block_samples <- seq(start(sink_viewport)[2L], end(sink_viewport)[2L])
        block_samples_df <- x_samples_df[
            x_samples_df[["sample"]] %in% block_samples, ]
        x_to_load <- unique(block_samples_df[["x_idx"]])
        block_x <- lapply(x_to_load, function(idx) {
            columns_to_load <- block_samples_df[
                block_samples_df[["x_idx"]] == idx, "column_idx"]
            as.matrix(x[[idx]][, columns_to_load, drop = FALSE])
        })
        block_rowRanges <- x_rowRanges[x_to_load]
        block_ans <- .combineList_matrix(
            x = block_x,
            x_rowRanges = block_rowRanges,
            ans_rowRanges = ans_rowRanges)
        write_block_to_sink(block_ans, sink, sink_viewport)
    }

    as(sink, "DelayedArray")
}

# Exported methods -------------------------------------------------------------

setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    combineList(list(x, y))
})

# Exported functions -----------------------------------------------------------

# TODO: Check sampleNames are disjoint
# TODO: Should this have a BACKEND argument?
combineList <- function(x, ..., BACKEND = NULL) {
    # Argument checks ----------------------------------------------------------

    # Check inputs are BSseq objects
    if (is(x, "BSseq")) x <- list(x, ...)
    x_is_BSseq <- vapply(x, is, logical(1L), "BSseq")
    stopifnot(isTRUE(all(x_is_BSseq)))
    # Check inputs are combinable
    x_trans <- lapply(x, getBSseq, "trans")
    x_has_same_trans <- vapply(
        x_trans,
        function(t) isTRUE(all.equal(t, x_trans[[1L]])),
        logical(1L))
    stopifnot(isTRUE(all(x_has_same_trans)))
    x_parameters <- lapply(x, getBSseq, "parameters")
    x_has_same_parameters <- vapply(
        x_parameters,
        function(p) isTRUE(all.equal(p, x_parameters[[1L]])),
        logical(1L))
    stopifnot(isTRUE(all(x_has_same_parameters)))
    # Check if all inputs have the same set of loci
    x_rowRanges <- lapply(x, rowRanges)
    ans_rowRanges <- Reduce(union, x_rowRanges)
    intersect_rowRanges <- Reduce(intersect, x_rowRanges)
    all_x_have_same_loci <- identical(ans_rowRanges, intersect_rowRanges)
    # Check if safe to combine smoothed representations (coef and se.coef).
    # It's only safe if all objects contain the same set of loci.
    x_has_been_smoothed <- vapply(x, hasBeenSmoothed, logical(1L))
    if (all_x_have_same_loci && all(x_has_been_smoothed)) {
        combine_smooths <- TRUE
        ans_trans <- x_trans[[1L]]
        ans_parameters <- x_parameters[[1L]]
    } else {
        if (any(x_has_been_smoothed)) {
            warning("Combining smoothed and unsmoothed BSseq objects. You ",
                    "will need to re-smooth using 'BSmooth()' on the ",
                    "returned object.")
        }
        if (!all_x_have_same_loci && any(x_has_been_smoothed)) {
            warning("Combining BSseq objects with different loci. You ",
                    "will need to re-smooth using 'BSmooth()' on the ",
                    "returned object.")
        }
        combine_smooths <- FALSE
        ans_coef <- NULL
        ans_se.coef <- NULL
        ans_trans <- function(x) NULL
        ans_parameters <- list()
    }

    # Combine elements of BSseq objects ----------------------------------------

    # Combine colData
    ans_colData <- as(
        Reduce(combine, lapply(x, function(xx) as.data.frame(colData(xx)))),
        "DataFrame")
    # Extract assays to be combined
    # NOTE: Use of `assay(..., withDimnames = FALSE)` is deliberate
    x_M <- lapply(x, assay, "M", withDimnames = FALSE)
    x_Cov <- lapply(x, assay, "Cov", withDimnames = FALSE)
    if (combine_smooths) {
        x_coef <- lapply(x, getBSseq, "coef")
        x_se.coef <- lapply(x, getBSseq, "se.coef")
    }
    # Combine assay data
    if (all_x_have_same_loci) {
        x_has_same_rowRanges_as_ans_rowRanges <- vapply(
            x_rowRanges,
            function(r) isTRUE(all.equal(r, ans_rowRanges)),
            logical(1L))
        if (all(x_has_same_rowRanges_as_ans_rowRanges)) {
            # If all rowRanges are identical then can just cbind assays.
            ans_M <- do.call(cbind, x_M)
            ans_Cov <- do.call(cbind, x_Cov)
            if (combine_smooths) {
                ans_coef <- do.call(cbind, x_coef)
                ans_se.coef <- do.call(cbind, x_se.coef)
            }
        } else {
            # If all rowRanges contain same loci but aren't identical then they
            # must differ by how loci are ordered. So order all inputs using a
            # common order and then cbind assays.
            x_order <- lapply(x_rowRanges, function(rr) {
                seqlevels(rr) <- seqlevels(ans_rowRanges)
                order(rr)
            })
            ans_M <- do.call(cbind, Map(extractROWS, x_M, x_order))
            ans_Cov <- do.call(cbind, Map(extractROWS, x_Cov, x_order))
            if (combine_smooths) {
                ans_M <- do.call(cbind, Map(extractROWS, x_M, x_order))
                ans_Cov <- do.call(cbind, Map(extractROWS, x_Cov, x_order))
            }
        }
    } else {
        ans_M <- .combineList(x_M, x_rowRanges, ans_rowRanges)
        ans_Cov <- .combineList(x_Cov, x_rowRanges, ans_rowRanges)
    }

    # Construct BSseq object ---------------------------------------------------

    ans_assays <- SimpleList(M = ans_M, Cov = ans_Cov)
    if (combine_smooths) {
        ans_assays <- c(ans_assays, SimpleList(coef = ans_coef))
        if (!is.null(ans_se.coef)) {
            ans_assays <- c(ans_assays, SimpleList(se.coef = ans_se.coef))
        }
    }
    se <- SummarizedExperiment(
        assays = ans_assays,
        rowRanges = ans_rowRanges,
        colData = ans_colData)
    .BSseq(se, parameters = ans_parameters, trans = ans_trans)
}
