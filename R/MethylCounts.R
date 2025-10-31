# Exported classes -------------------------------------------------------------

# NOTE: This create a 'classGeneratorFunction' (internal constructor), .BSseq()
.MethylCounts <- setClass(
    "MethylCounts",
    slots = representation(
        parameters = "list"),
    contains = "RangedSummarizedExperiment"
)

# Validity methods -------------------------------------------------------------

# TODO: Benchmark validity method
setValidity2("MethylCounts", function(object) {
    msg <- NULL

    if (identical(object, .BSseq())) {
        # No validity checks for object returned by internal constructor
        return(msg)
    }
    msg <- validMsg(msg, .checkAssayNames(object, c("M", "U")))

    ## TODO: Add a check like .checkMandCov
    ##       It should check that M, U, D (optionally H) are integers >= 0
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# Exported functions -----------------------------------------------------------

# TODO: BSseq() is arguably a bad constructor. It doesn't return a valid BSseq
#       object when called without any arguments. It also does some pretty
#       complicated parsing of the inputs. But we're stuck with it because it's
#       been around for a long time.

MethylCounts <- function(M = NULL, U = NULL, D = NULL, H = NULL,
                     parameters = NULL, colData = NULL, gr = NULL,
                     pos = NULL, chr = NULL, sampleNames = NULL,
                     rmZeroCov = FALSE) {

    # Argument checks ----------------------------------------------------------

    # Process assays.
    # NOTE: Nothing to do for 'coef', and 'se.coef'.
    if (is.null(M) || is.null(U)) {
        stop("Need 'M' and 'U'.")
    }
    # Process 'parameters'.
    if (is.null(parameters)) {
        parameters <- list()
    }
    # Process 'sampleNames' and 'colData'.
    if (is.null(sampleNames)) {
        if (is.null(colData)) {
            ## BSseq object will have no colnames.
            colData <- make_zero_col_DFrame(ncol(M))
        } else {
            # BSseq object will have 'sampleNames' as colnames.
            colData <- DataFrame(row.names = sampleNames)
        }
    } else {
        if (is.null(colData)) {
            # BSseq object will have 'sampleNames' as colnames.
            colData <- DataFrame(row.names = sampleNames)
        } else {
            if (is.null(rownames(colData))) {
                rownames(colData) <- sampleNames
            } else {
                stopifnot(identical(rownames(colData), sampleNames))
            }
        }
    }
    # Process 'gr', 'pos', and 'chr'.
    if (is.null(gr)) {
        if (is.null(pos) || is.null(chr)) {
            stop("Need 'pos' and 'chr' if 'gr' not supplied.")
        }
        gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1L))
    }
    if (!is(gr, "GRanges")) {
        stop("'gr' needs to be a GRanges.")
    }
    # Process 'rmZeroCov'.
    stopifnot(isTRUEorFALSE(rmZeroCov))

    # Collapse duplicate loci --------------------------------------------------

    is_duplicated <- duplicated(gr)
    if (any(is_duplicated)) {
        warning("Detected duplicate loci. Collapsing counts in 'M' and 'U'",
                "(optionally 'H' and 'D') at these positions.")
        loci <- gr[!is_duplicated]
        ol <- findOverlaps(gr, loci, type = "equal")
        M <- rowsum(x = M, group = subjectHits(ol), reorder = FALSE)
        rownames(M) <- NULL
        U <- rowsum(x = U, group = subjectHits(ol), reorder = FALSE)
        rownames(Cov) <- NULL
        if(!is.null(H)) {
            H <- rowsum(x = H, group = subjectHits(ol), reorder = FALSE)
            rownames(H) <- NULL
        }
        if(!is.null(D)) {
            D <- rowsum(x = D, group = subjectHits(ol), reorder = FALSE)
            rownames(D) <- NULL
        }
    } else {
        loci <- gr
    }

    # Optionally, remove positions with zero coverage --------------------------

    if (rmZeroCov) {
        Cov <- M + U
        if(!is.null(H)) {
            Cov <- Cov + H
        }
        loci_with_zero_cov <- rowAlls(Cov, value = 0)
        if (any(loci_with_zero_cov)) {
            loci_with_nonzero_cov <- !loci_with_zero_cov
            gr <- gr[loci_with_nonzero_cov]
            M <- M[loci_with_nonzero_cov, , drop = FALSE]
            U <- U[loci_with_nonzero_cov, , drop = FALSE]
            if(!is.null(H)) {
                H <- H[loci_with_nonzero_cov, , drop = FALSE]
            }
            if(!is.null(D)) {
                D <- D[loci_with_nonzero_cov, , drop = FALSE]
            }
        }
    }

    # Construct BSseq object ---------------------------------------------------

    se <- SummarizedExperiment(
        assays = SimpleListExcludeNULL(M = M, U = U, H = H, D = D),
        rowRanges = loci,
        colData = colData)
    .MethylCounts(se, parameters = parameters)
}

## Move to BSseq-utils?
getMethylCounts <- function(MethylCounts,
                        type = c("M", "U", "H", "D", "gr", "parameters"),
                        withDimnames = TRUE) {
    type <- match.arg(type)
    if (type %in% c("M", "U", "H", "D")) {
        return(assay(MethylCounts, type, withDimnames = withDimnames))
    }
    if (type == "parameters") {
        return(MethylCounts@parameters)
    }
    if (type == "gr") {
        return(MethylCounts@rowRanges)
    }
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "MethylCounts"), function(object) {
    cat("An object of type 'MethylCounts' with\n")
    cat(" ", nrow(object), "methylation loci\n")
    cat(" ", ncol(object), "samples\n")
    if (.isHDF5ArrayBacked(object)) {
        cat("Some assays are HDF5Array-backed\n")
    } else {
        cat("All assays are in-memory\n")
    }
})

