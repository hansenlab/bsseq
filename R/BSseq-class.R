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
    validMsg(msg, .Call(cxx_check_M_and_Cov, initializeCpp(M), initializeCpp(Cov), 1))
}

# TODO: Benchmark validity method
setValidity2("BSseq", function(object) {
    msg <- NULL

    if (identical(object, .BSseq())) {
        # No validity checks for object returned by internal constructor
        return(msg)
    }
    msg <- validMsg(msg, .checkAssayNames(object, c("Cov", "M")))

    msg <- validMsg(
        msg = msg,
        result = .checkMandCov(
            M = assay(object, "M", withDimnames = FALSE),
            Cov = assay(object, "Cov", withDimnames = FALSE)))

    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
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
#       complicated parsing of the inputs. But we're stuck with it because it's
#       been around for a long time.

BSseq <- function(mc = NULL, M = NULL, Cov = NULL, coef = NULL, se.coef = NULL,
                  trans = NULL, parameters = NULL, pData = NULL, gr = NULL,
                  pos = NULL, chr = NULL, sampleNames = NULL,
                  rmZeroCov = FALSE, mods = c("5mC+5hmC", "5mC", "5hmC")) {

    # Add MethylCounts as input ------------------------------------------------
    if (inherits(mc, "MethylCounts")) {
        # Match the conversionType
        mods <- match.arg(mods)
        gr <- getMethylCounts(mc, type = "gr")
        pData <- colData(mc)
        sampleNames <- colnames(mc)
        Cov <- getMethylCounts(mc, type = "Cov")
        M <- getMethylCounts(mc, type = "M")
        H <- getMethylCounts(mc, type = "H")
        # Compute M and Cov based on conversionType
        if (mods == "5mC+5hmC") {
            M <- M + H
        } else if (mods == "5mC") {
            M <- M
        } else if (mods == "5hmC") {
            M <- H
        }
    pData$Mod <- rep(mods,ncol(mc))
    }

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
        }
    }

    # Construct BSseq object ---------------------------------------------------

    se <- SummarizedExperiment(
        assays = SimpleListExcludeNULL(M = M, Cov = Cov,
                         coef = coef, se.coef = se.coef),
        rowRanges = loci,
        colData = pData)
    .BSseq(se, trans = trans, parameters = parameters)
}

# Move to BSseq-utils?
hasBeenSmoothed <- function(BSseq) {
    "coef" %in% assayNames(BSseq)
}

# Move to BSseq-utils?
getBSseq <- function(BSseq,
                     type = c("Cov", "M", "gr","coef", "se.coef","trans",
                              "parameters"),
                     withDimnames = TRUE) {
    type <- match.arg(type)
    if (type %in% c("M", "Cov")) {
        return(assay(BSseq, type, withDimnames = withDimnames))
    }
    if (type %in% c("coef", "se.coef")) {
        if (type %in% assayNames(BSseq)) {
            return(assay(BSseq, type, withDimnames = withDimnames))
        } else {
            return(NULL)
        }
    }
    if (type == "trans") {
        return(BSseq@trans)
    }
    if (type == "parameters") {
        return(BSseq@parameters)
    }
    if (type == "gr") {
        return(BSseq@rowRanges)
    }
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "BSseq"), function(object) {
    cat("An object of type 'BSseq' with\n")
    cat(" ", nrow(object), "loci\n")
    cat(" ", ncol(object), "samples\n")

    # Check if Mod column in colData contains valid entries
    if (!is.null(colData(object)$Mod)) {
        cat(" ", unique(colData(object)$Mod), "values are stored in 'M'\n")
    }

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
              if (hasBeenSmoothed(object) &&
                  isTRUE(all.equal(getBSseq(object, "trans"), .oldTrans))) {
                  object@trans <- plogis
              }
              if (is(object, "SummarizedExperiment")) {
                  # NOTE: Call method for SummarizedExperiment objects
                  object <- callNextMethod()
                  return(object)
              } else {
                  BSseq(
                      gr = object@gr,
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
