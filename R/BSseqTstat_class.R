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
          function(object, ...) {
              stats <- object@stats
              stats <- stats
              object@stats <- stats
              object
          }
)
