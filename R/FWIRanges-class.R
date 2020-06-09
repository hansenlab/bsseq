# Internal classes -------------------------------------------------------------

# The FWIRanges class is a simple container for storing a vector of fixed-width
# integer ranges.
# NOTE: The intention is to make this a fully-fledged IntegerRanges subclass
#       that is part of the IRanges package. For now, this class is internal to
#       bsseq and only a subset of methods are properly implemented.
.FWIRanges <- setClass(
    "FWIRanges",
    contains = "IPosRanges",
    representation(
        start = "integer",
        width = "integer",
        NAMES = "character_OR_NULL"  # R doesn't like @names !!
    )
)

# # Validity methods -------------------------------------------------------------

# TODO: Some of these may be redundant
.valid.FWIRanges <- function(x) {
    msg <- NULL
    if (!is.integer(x@width)) {
        msg <- validMsg(msg, "'width' must be an integer.")
    }
    if (length(x@width) != 1L ) {
        msg <- validMsg(msg, "'width' must have length 1.")
    }
    if (x@width < 0) {
        msg <- validMsg(msg, "'width' must be non-negative.")
    }
}

setValidity2("FWIRanges", .valid.FWIRanges)

# Internal methods -------------------------------------------------------------

setMethod("parallel_slot_names", "FWIRanges",
          function(x) c("start", "NAMES", callNextMethod())
)

setMethod("start", "FWIRanges", function(x, ...) x@start)

setMethod("width", "FWIRanges", function(x) rep.int(x@width, length(x)))

setMethod("end", "FWIRanges", function(x) {
    if (x@width == 1L) {
        x@start
    } else {
        width(x) - 1L + start(x)
    }
})

setMethod("names", "FWIRanges", function(x) x@NAMES)

# TODO: Room for optmisation (e.g., write in C to reduce memory allocations).
.set_FWIRanges_start <- function(x, value, check = TRUE) {
    if (!isTRUEorFALSE(check)) stop("'check' must be TRUE or FALSE")
    old_start <- start(x)
    new_start <- S4Vectors:::numeric2integer(value)
    new_width <- x@width - new_start + old_start
    if (any(new_width != new_width[1L])) {
        x <- as(x, "IRanges")
        start(x) <- value
        return(x)
    }
    x@width <- new_width[1L]
    if (check) validObject(x)
    x
}

setReplaceMethod(
    "start",
    "FWIRanges",
    function(x, ..., value) .set_FWIRanges_start(x, value)
)

.set_FWIRanges_end <- function(x, value, check = TRUE) {
    if (!isTRUEorFALSE(check)) stop("'check' must be TRUE or FALSE")
    new_width <- x@width + S4Vectors:::numeric2integer(value) - end(x)
    if (any(new_width != new_width[1L])) {
        x <- as(x, "IRanges")
        end(x) <- value
        return(x)
    }
    x@width <- new_width[1L]
    if (check) validObject(x)
    x
}

setReplaceMethod(
    "end",
    "FWIRanges",
    function(x, ..., value) .set_FWIRanges_end(x, value)
)

.set_FWIRanges_width <- function(x, value, check = TRUE) {
    if (!isTRUEorFALSE(check)) stop("'check' must be TRUE or FALSE")
    new_width <- S4Vectors:::numeric2integer(value)
    if (any(new_width != new_width[1L])) {
        x <- as(x, "IRanges")
        width(x) <- value
        return(x)
    }
    x@width <- new_width[1L]
    if (check) validObject(x)
    x
}

setReplaceMethod(
    "width",
    "FWIRanges",
    function(x, ..., value) .set_FWIRanges_width(x, value)
)

set_FWIRanges_names <- function(x, value) {
    x@NAMES <- S4Vectors:::normarg_names(value, class(x), length(x))
    # NOTE: No need to validate an FWIRanges object after setting its names so
    #       this should be safe.
    x
}

setReplaceMethod("names", "FWIRanges", set_FWIRanges_names)

setMethod("replaceROWS", "FWIRanges", function(x, i, value) {
    x_width <- x@width
    value_width <- value@width
    if (!identical(x_width, value_width)) {
        x <- as(x, "IRanges")
        value <- as(value, "IRanges")
        return(replaceROWS(x, i, value))
    }
    i <- normalizeSingleBracketSubscript(i, x, as.NSBS = TRUE)
    ans_start <- replaceROWS(start(x), i, start(value))
    ans_width <- value@width
    ans_mcols <- replaceROWS(mcols(x), i, mcols(value))
    BiocGenerics:::replaceSlots(
        x,
        start = ans_start,
        width = ans_width,
        mcols = ans_mcols,
        check = FALSE)
})

# TODO: Follow shift,IRanges-method
setMethod("shift", "FWIRanges", function(x, shift = 0L, use.names = TRUE) {
    stopifnot(use.names)
    shift <- recycleIntegerArg(shift, "shift", length(x))
    new_start <- start(x) + shift
    # TODO: Use BiocGenerics:::replaceSlots()
    x@start <- new_start
    validObject(x)
    x
})
