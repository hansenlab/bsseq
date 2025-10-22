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
                     parameters = NULL, pData = NULL, gr = NULL,
                     pos = NULL, chr = NULL, sampleNames = NULL,
                     rmZeroCov = FALSE) {

    # Argument checks ----------------------------------------------------------

    # Process assays.
    # NOTE: Nothing to do for 'coef', and 'se.coef'.
    if (is.null(M) || is.null(Cov)) {
        stop("Need 'M' and 'Cov'.")
    }
    # Process 'trans' and 'parameters'.
    if (is.null(trans)) {
        trans <- function() NULL
        environment(trans) <- emptyenv()
    }
    if (is.null(parameters)) {
        parameters <- list()
    }
    # Process 'sampleNames' and 'pData'.
    if (is.null(sampleNames)) {
        if (is.null(pData)) {
            # BSseq object will have no colnames.
            pData <- make_zero_col_DFrame(ncol(M))
        } else {
            # BSseq object will have 'sampleNames' as colnames.
            pData <- DataFrame(row.names = sampleNames)
        }
    } else {
        if (is.null(pData)) {
            # BSseq object will have 'sampleNames' as colnames.
            pData <- DataFrame(row.names = sampleNames)
        } else {
            if (is.null(rownames(pData))) {
                rownames(pData) <- sampleNames
            } else {
                stopifnot(identical(rownames(pData), sampleNames))
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
        warning("Detected duplicate loci. Collapsing counts in 'M' and 'Cov' ",
                "at these positions.")
        if (!is.null(coef) || !is.null(se.coef)) {
            stop("Cannot collapse when 'coef' or 'se.coef' are non-NULL.")
        }
        loci <- gr[!is_duplicated]
        ol <- findOverlaps(gr, loci, type = "equal")
        M <- rowsum(x = M, group = subjectHits(ol), reorder = FALSE)
        rownames(M) <- NULL
        Cov <- rowsum(x = Cov, group = subjectHits(ol), reorder = FALSE)
        rownames(Cov) <- NULL
        if (!is.null(Filtered)){
          Filtered <- rowsum(x = Filtered, group = subjectHits(ol), reorder = FALSE)
          rownames(Filtered) <- NULL}
    } else {
        loci <- gr
    }

    # Optionally, remove positions with zero coverage --------------------------

    if (rmZeroCov) {
        loci_with_zero_cov <- rowAlls(Cov, value = 0)
        if (any(loci_with_zero_cov)) {
            loci_with_nonzero_cov <- !loci_with_zero_cov
            gr <- gr[loci_with_nonzero_cov]
            M <- M[loci_with_nonzero_cov, , drop = FALSE]
            Cov <- Cov[loci_with_nonzero_cov, , drop = FALSE]
            coef <- coef[loci_with_nonzero_cov, , drop = FALSE]
            se.coef <- se.coef[loci_with_nonzero_cov, , drop = FALSE]
            Filtered <- Filtered[loci_with_nonzero_cov, , drop = FALSE]
        }
    }

    # Construct BSseq object ---------------------------------------------------

    assays <- SimpleList(M = M, Cov = Cov, Filtered = Filtered,
                         coef = coef, se.coef = se.coef)
    assays <- assays[!S4Vectors:::sapply_isNULL(assays)]
    se <- SummarizedExperiment(
        assays = assays,
        rowRanges = loci,
        colData = pData)
    .BSseq(se, trans = trans, parameters = parameters)
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

setMethod("pData", "MethylCounts", function(object) {
    object@colData
})

setReplaceMethod(
    "pData",
    signature = signature(object = "MethylCounts", value = "data.frame"),
    function(object, value) {
        colData(object) <- as(value, "DataFrame")
        object
    }
)

setReplaceMethod(
    "pData",
    signature = signature(object = "MethylCounts", value = "DataFrame"),
    function(object, value) {
        colData(object) <- value
        object
    }
)

setMethod("sampleNames", "MethylCounts", function(object) colnames(object))

setReplaceMethod(
    "sampleNames",
    signature = signature(object = "MethylCounts", value = "ANY"),
    function(object, value) {
        colnames(object) <- value
        object
    }
)
