## Add new union members as needed (e.g. "SparseMatrix").
setClassUnion("matrix_OR_DelayedMatrix", c("matrix", "DelayedMatrix"))

setClass("BSseqTstat", contains = "GRanges",
         representation(stats = "matrix_OR_DelayedMatrix",
                        parameters = "list")
         )
setValidity("BSseqTstat", function(object) {
    msg <- NULL
    if(length(object) != nrow(object@stats))
        msg <- c(msg, "'stats' must have one row per methylation locus")
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
    i <- normalizeSingleBracketSubscript(i, x)
    x@stats <- x@stats[i,, drop = FALSE]
    callNextMethod()
})

BSseqTstat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    new("BSseqTstat", gr, stats = stats, parameters = parameters)
}

summary.BSseqTstat <- function(object, ...) {
    quant <- .quantile(getStats(object)[, "tstat.corrected"],
                       prob = c(0.0001, 0.001, 0.01, 0.5, 0.99, 0.999, 0.9999))
    quant <- t(t(quant))
    colnames(quant) <- "quantiles"
    out <- list(quantiles = quant)
    class(out) <- "summary.BSseqTstat"
    out
}

print.summary.BSseqTstat <- function(x, ...) {
    print(as.matrix(x$quantiles))
}

plot.BSseqTstat <- function(x, y, ...) {
    tstat <- getStats(x)[, "tstat"]
    tstat <- as.array(tstat)
    plot(density(tstat), xlim = c(-10,10), col = "blue", main = "")
    if("tstat.corrected" %in% colnames(getStats(x))) {
        tstat.cor <- getStats(x)[, "tstat.corrected"]
        tstat.cor <- as.array(tstat.cor)
        lines(density(tstat.cor), col = "black")
        legend("topleft", legend = c("uncorrected", "corrected"), lty = c(1,1),
               col = c("blue", "black"))
    } else {
        legend("topleft", legend = c("uncorrected"), lty = 1,
               col = c("blue"))
    }
}

setMethod("updateObject", "BSseqTstat",
    function(object, ..., verbose = FALSE) {
        if (!.hasSlot(object, "gr")) {
            return(callNextMethod())
        }
        if (verbose)
            message("[updateObject] ", class(object), " object ",
                    "uses old internal representation from\n",
                    "[updateObject] bsseq <= 1.49.2. ",
                    "Updating it ... ", appendLF=FALSE)
        gr <- updateObject(object@gr)
        ans <- BSseqTstat(gr = gr, stats = object@stats,
                                   parameters = object@parameters)
        if (verbose)
            message("OK")
        ans
    }
)
