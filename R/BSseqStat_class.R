setClass("BSseqStat", contains = "hasGRanges",
         representation(stats = "list",
                        parameters = "list")
         )

setValidity("BSseqStat", function(object) {
    msg <- NULL
    if(is.null(names(object@stats)) || any(names(object@stats) == "") ||
       anyDuplicated(names(object@stats)))
        msg <- validMsg(msg, "the 'stats' list needs to be named with unique names.")
    for(name in c("rawSds", "smoothsSds", "stat", "rawTstats")) {
        if(name %in% names(object@stats) && isTRUE(nrow(object@stats[[name]]) != length(object@gr)))
            msg <- validMsg(msg, sprintf("component '%s' of slot 'stats' has to have the same number of rows as slot 'gr' is long", name))
    }
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqStat"),
          function(object) {
              cat("An object of type 'BSseqStat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
          })

setMethod("[", "BSseqStat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    statnames <- names(x@stats)
    names(statnames) <- statnames
    x@stats <- lapply(statnames, function(nam) {
        if(nam %in% c("rawTstats", "modelCoefficients", "rawSds", "smoothSds",
                      "stat")) {
            return(x@stats[[nam]][i,,drop=FALSE])
        }
        x@stats[[nam]]
    })
    x
})

BSseqStat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqStat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}

# TODO: updateObject() to use ordinary matrix instead of DelayedMatrix with
#       in-memory seed.
# setMethod("updateObject", "BSseqStat",
#           function(object, ...) {
#               not_delayed_matrices <- c("cor.coefficients", "stat.type")
#               delayed_matrices <- setdiff(names(stats), not_delayed_matrices)
#               stats <- object@stats
#               stats[delayed_matrices] <- endoapply(stats[delayed_matrices],
#                                                    .DelayedMatrix)
#               object@stats <- stats
#               object
#           }
# )

## summary.BSseqStat <- function(object, ...) {
##     quant <- quantile(getStats(object)[, "tstat.corrected"],
##                       prob = c(0.0001, 0.001, 0.01, 0.5, 0.99, 0.999, 0.9999))
##     quant <- t(t(quant))
##     colnames(quant) <- "quantiles"
##     out <- list(quantiles = quant)
##     class(out) <- "summary.BSseqStat"
##     out
## }

## print.summary.BSseqStat <- function(x, ...) {
##     print(as.matrix(x$quantiles))
## }

## plot.BSseqStat <- function(x, y, ...) {
##     tstat <- getStats(x)[, "tstat"]
##     plot(density(tstat), xlim = c(-10,10), col = "blue", main = "")
##     if("tstat.corrected" %in% colnames(getStats(x))) {
##         tstat.cor <- getStats(x)[, "tstat.corrected"]
##         lines(density(tstat.cor), col = "black")
##         legend("topleft", legend = c("uncorrected", "corrected"), lty = c(1,1),
##                col = c("blue", "black"))
##     } else {
##         legend("topleft", legend = c("uncorrected"), lty = 1,
##                col = c("blue"))
##     }
## }

