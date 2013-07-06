setClass("BSseqTstat", contains = "hasGRanges", 
         representation(stats = "matrix",
                        parameters = "list")
         )
setValidity("BSseqTstat", function(object) {
    msg <- NULL
    if(length(object@gr) != nrow(object@stats))
        msg <- c(msg, "length of 'gr' is different from the number of rows of 'stats'")
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqTstat"),
          function(object) {
              cat("An object of type 'BSseqTstat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
              cat("with parameters\n")
              cat(" ", object@parameters$tstatText, "\n")
          })

setMethod("[", "BSseqTstat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x@stats <- x@stats[i,, drop = FALSE]
    x
})

BSseqTstat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqTstat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}

