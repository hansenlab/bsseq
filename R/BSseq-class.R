# Exported classes -------------------------------------------------------------

# NOTE: This create a 'classGeneratorFunction' (internal constructor), .BSseq()
.BSseq <- setClass(
    "BSseq",
    slots = representation(
        trans = "function",
        parameters = "list"),
    contains = "RangedSummarizedExperiment"
)

# Validity methods -------------------------------------------------------------

.checkAssayNames <- function(object, names) {
    if (all(names %in% assayNames(object))) {
        return(NULL)
    } else {
        paste0(
            "object of class ",
            class(object),
            " needs to have assay slots with names ",
            paste0(names, collapse = ", "))

    }
}

.checkMandCov <- function(M, Cov) {
    msg <- NULL
    validMsg(msg, .Call(cxx_check_M_and_Cov, M, Cov))
}

setValidity2("BSseq", function(object) {
    msg <- NULL

    if (identical(object, .BSseq())) {
        # No validity checks for object returned by internal constructor
        return(msg)
    }
    msg <- validMsg(msg, .checkAssayNames(object, c("Cov", "M")))

    # TODO: Are colnames strictly necessary?
    if (is.null(colnames(object))) {
        msg <- validMsg(msg, "colnames (aka sampleNames) need to be set")
    }

    assay_rownames <- lapply(assays(object), rownames)
    if (!all(S4Vectors:::sapply_isNULL(assay_rownames))) {
        msg <- validMsg(msg, "unnecessary rownames on one-or-more assays")
    }

    msg <- validMsg(msg,
                    .checkMandCov(assay(object, "M", withDimnames = FALSE),
                                  assay(object, "Cov", withDimnames = FALSE)))

    if (is.null(msg)) TRUE else msg
})

# Internal functions -----------------------------------------------------------

.oldTrans <- function(x) {
    y <- x
    ix <- which(x < 0)
    ix2 <- which(x > 0)
    y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
    y[ix2] <- 1/(1 + exp(-x[ix2]))
    y
}

# Exported functions -----------------------------------------------------------

# TODO: BSseq() is arguably a bad constructor. It doesn't return a valid BSseq
#       object when called without any arguments. It also does some pretty
#       complicated parsing of the inputs. Still, I think we're stuck with it
#       because it's been around for a long time.
BSseq <- function(M = NULL, Cov = NULL, coef = NULL, se.coef = NULL,
                  trans = NULL, parameters = NULL, pData = NULL,
                  gr = NULL, pos = NULL, chr = NULL, sampleNames = NULL,
                  rmZeroCov = FALSE) {

    # Process 'gr', or 'pos' and 'chr'
    if (is.null(gr)) {
        if (is.null(pos) || is.null(chr)) {
            stop("Need 'pos' and 'chr' if 'gr' not supplied.")
        }
        gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1L))
    }
    if (!is(gr, "GRanges")) stop("'gr' needs to be a GRanges.")
    if (any(width(gr) != 1L)) stop("'gr' needs to have widths of 1.")

    # Process 'M' and 'Cov'
    if (is.null(M) || is.null(Cov)) stop("Need 'M' and 'Cov'.")
    if (length(gr) != nrow(M) ||
        length(gr) != nrow(Cov) ||
        ncol(Cov) != ncol(M)) {
        stop("'gr', 'M' and 'Cov' need to have compatible dimensions.")
    }
    if (!is.null(rownames(M))) rownames(M) <- NULL
    if (!is.null(rownames(Cov))) rownames(Cov) <- NULL
    if (!is.null(names(gr))) names(gr) <- NULL

    # Process 'sampleNames' and 'pData'
    if (is.null(pData)) {
        if (is.null(sampleNames)) {
            if (!is.null(colnames(M))) {
                sampleNames <- colnames(M)
            } else if (!is.null(colnames(Cov))) {
                sampleNames <- colnames(Cov)
            } else {
                sampleNames <- paste0("V", seq_len(ncol(M)))
            }
        }
        pData <- DataFrame(row.names = sampleNames)
    } else {
        pData <- as(pData, "DataFrame")
    }
    if (is.null(sampleNames) && !is.null(pData) && !is.null(rownames(pData))) {
        sampleNames <- rownames(pData)
    }
    if (length(unique(sampleNames)) != ncol(M)) {
        stop("sampleNames need to be unique and of the right length.")
    }

    # TODO: Is this really necessary? It complicates things because HDF5Matrix
    #       will be degraded to a DelayedMatrix
    # Add column names to assays
    if (is.null(colnames(M)) || any(sampleNames != colnames(M))) {
        colnames(M) <- sampleNames
    }
    if (is.null(colnames(Cov)) || any(sampleNames != colnames(Cov))) {
        colnames(Cov) <- sampleNames
    }
    if (!is.null(coef)) {
        if (nrow(coef) != nrow(M) || ncol(coef) != ncol(M)) {
            stop("'coef' does not have the right dimensions")
        }
        if (is.null(colnames(coef)) || any(sampleNames != colnames(coef))) {
            colnames(coef) <- sampleNames
        }
        if (!is.null(rownames(coef))) {
            rownames(coef) <- NULL
        }
    }
    if (!is.null(se.coef)) {
        if (nrow(se.coef) != nrow(M) || ncol(se.coef) != ncol(M)) {
            stop("'se.coef' does not have the right dimensions")
        }
        if (is.null(colnames(se.coef)) ||
            any(sampleNames != colnames(se.coef))) {
            colnames(se.coef) <- sampleNames
        }
        if (!is.null(rownames(se.coef))) {
            rownames(se.coef) <- NULL
        }
    }

    # Optionally, remove positions with Cov == 0
    if (rmZeroCov) {
        loci_with_zero_cov <- rowAlls(Cov, value = 0)
        if (any(loci_with_zero_cov)) {
            gr <- gr[!loci_with_zero_cov]
            M <- M[!loci_with_zero_cov, , drop = FALSE]
            Cov <- Cov[!loci_with_zero_cov, , drop = FALSE]
        }
    }

    # Collapse duplicate loci
    if (any(duplicated(gr))) {
        warning("Detected duplicate loci. Collapsing counts in 'M' and 'Cov' ",
                "at these positions.")
        if (!is.null(coef) || !is.null(se.coef)) {
            stop("Cannot collapse when 'coef' or 'se.coef' are present.")
        }
        # NOTE: reduce() sorts the output
        unique_gr <- unique(gr)
        ol <- findOverlaps(unique_gr, gr)
        gr <- unique(gr)
        idx <- as.list(ol)
        M <- .collapseMatrixLike(M, idx, MARGIN = 2L)
        Cov <- .collapseMatrixLike(Cov, idx, MARGIN = 2L)
    }

    # Sort loci
    if (is.unsorted(gr)) {
        o <- order(gr)
        gr <- gr[o]
        M <- M[o, , drop = FALSE]
        Cov <- Cov[o, , drop = FALSE]
        coef <- coef[o, , drop = FALSE]
        se.coef <- se.coef[o, , drop = FALSE]
    }

    # Construct BSseq object
    assays <- SimpleList(M = M, Cov = Cov, coef = coef, se.coef = se.coef)
    assays <- assays[!S4Vectors:::sapply_isNULL(assays)]
    se <- SummarizedExperiment(assays = assays, rowRanges = gr, colData = pData)
    if (is.null(parameters)) {
        parameters <- list()
    }
    if (is.null(trans)) {
        trans <- function(x) NULL
    }
    .BSseq(se, trans = trans, parameters = parameters)
}

hasBeenSmoothed <- function(BSseq) "coef" %in% assayNames(BSseq)

# TODO: Document withDimnames
getBSseq <- function(BSseq,
                     type = c("Cov", "M", "gr", "coef", "se.coef", "trans",
                              "parameters"),
                     withDimnames = TRUE) {
    type <- match.arg(type)
    if (type %in% c("M", "Cov")) {
        return(assay(BSseq, type, withDimnames = withDimnames))
    }
    if (type %in% c("coef", "se.coef") && type %in% assayNames(BSseq)) {
        return(assay(BSseq, type, withDimnames = withDimnames))
    }
    if (type %in% c("coef", "se.coef")) return(NULL)
    if (type == "trans") return(BSseq@trans)
    if (type == "parameters") return(BSseq@parameters)
    if (type == "gr") return(BSseq@rowRanges)

}

strandCollapse <- function(BSseq, shift = TRUE) {
    if (all(runValue(strand(BSseq)) == "*")) {
        warning("All loci are unstranded; nothing to collapse")
        return(BSseq)
    }
    if (!(all(runValue(strand(BSseq)) %in% c("+", "-")))) {
        stop("'BSseq' object has a mix of stranded and unstranded loci.")
    }
    BS.forward <- BSseq[strand(BSseq) == "+"]
    strand(BS.forward) <- "*"
    BS.reverse <- BSseq[strand(BSseq) == "-"]
    strand(BS.reverse) <- "*"
    if (shift) rowRanges(BS.reverse) <- shift(rowRanges(BS.reverse), -1L)
    sampleNames(BS.reverse) <- paste0(sampleNames(BS.reverse), "_REVERSE")
    BS.comb <- combine(BS.forward, BS.reverse)
    newBSseq <- collapseBSseq(BS.comb, columns = rep(sampleNames(BSseq), 2))
    newBSseq
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "BSseq"), function(object) {
    cat("An object of type 'BSseq' with\n")
    cat(" ", nrow(object), "methylation loci\n")
    cat(" ", ncol(object), "samples\n")
    if (hasBeenSmoothed(object)) {
        cat("has been smoothed with\n")
        cat(" ", object@parameters$smoothText, "\n")
    } else {
        cat("has not been smoothed\n")
    }
    if (.isHDF5ArrayBacked(object)) {
        cat("Some assays are HDF5Array-backed\n")
    } else {
        cat("All assays are in-memory\n")
    }
})

setMethod("pData", "BSseq", function(object) {
    object@colData
})

setReplaceMethod(
    "pData",
    signature = signature(object = "BSseq", value = "data.frame"),
    function(object, value) {
        colData(object) <- as(value, "DataFrame")
        object
    }
)

setReplaceMethod(
    "pData",
    signature = signature(object = "BSseq", value = "DataFrame"),
    function(object, value) {
        colData(object) <- value
        object
    }
)

setMethod("sampleNames", "BSseq", function(object) colnames(object))

setReplaceMethod(
    "sampleNames",
    signature = signature(object = "BSseq", value = "ANY"),
    function(object, value) {
        colnames(object) <- value
        object
    }
)

setMethod("updateObject", "BSseq",
          function(object, ...) {
              # NOTE: identical() is too strong
              if (isTRUE(all.equal(getBSseq(object, "trans"), .oldTrans))) {
                  object@trans <- plogis
              }
              if (is(object, "SummarizedExperiment")) {
                  # NOTE: Call method for SummarizedExperiment objects
                  object <- callNextMethod()
                  return(object)
              } else {
                  BSseq(gr = object@gr,
                        M = object@M,
                        Cov = object@Cov,
                        coef = object@coef,
                        se.coef = object@se.coef,
                        trans = object@trans,
                        parameters = object@parameters,
                        pData = object@phenoData@data)
              }
          }
)
