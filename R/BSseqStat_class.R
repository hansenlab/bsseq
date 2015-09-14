setClass("BSseqStat", contains = "hasGRanges", 
         representation(stats = "list",
                        parameters = "list")
         )
setValidity("BSseqStat", function(object) {
    msg <- NULL
    statsClasses <- sapply(object@stats, function(xx) is.vector(xx) || is.matrix(xx))
    if(!all(statsClasses))
        msg <- validMsg(msg, sprintf("All classes in the `stats` list has to be either vectors or matrices"))
    lengths <- sapply(object@stats, function(xx) {
        if(is.matrix(xx))
            return(nrow(xx))
        if(is.vector(xx))
            return(length(xx))
    })
    if(!all(length(object@gr) == statsLengths))
        msg <- validMsg(msg, sprintf("All components in the `stats` list needs to have the same length or number of rows as the `gr` slot"))
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqStat"),
          function(object) {
              cat("An object of type 'BSseqStat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
              cat("with parameters\n")
              cat(" ", object@parameters$tstatText, "\n")
          })

## setMethod("[", "BSseqStat", function(x, i, ...) {
##     if(missing(i))
##         stop("need [i] for subsetting")
##     if(missing(i))
##         return(x)
##     x@gr <- x@gr[i]
##     x@stats <- x@stats[i,, drop = FALSE]
##     x
## })

BSseqStat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqStat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}

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

## getStats <- function(BSseqStat, regions = NULL, stat = "tstat.corrected") {
##     stopifnot(is(BSseqStat, "BSseqStat"))
##     if(is.null(regions))
##         return(BSseqStat@stats)
##     if(class(regions) == "data.frame")
##         regions <- data.frame2GRanges(regions)
##     stopifnot(stat %in% colnames(BSseqStat@stats))
##     stopifnot(length(stat) == 1)
##     stopifnot(is(regions, "GenomicRanges"))
##     ov <- findOverlaps(BSseqStat, regions)
##     ov.sp <- split(queryHits(ov), subjectHits(ov))
##     getRegionStats <- function(idx) {
##         mat <- BSseqStat@stats[idx,, drop=FALSE]
##         areaStat <- sum(mat[, stat])
##         maxStat <- max(mat[, stat])
##         c(areaStat, maxStat)
##     }
##     stats <- matrix(NA, ncol = 2, nrow = length(regions))
##     colnames(stats) <- c("areaStat", "maxStat")
##     tmp <- lapply(ov.sp, getRegionStats)
##     stats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
##     out <- as.data.frame(stats)
##     if(! stat %in% c("tstat.corrected", "tstat"))
##         return(out)
##     getRegionStats_ttest <- function(idx) {
##         mat <- BSseqStat@stats[idx,, drop=FALSE]
##         group1.mean <- mean(mat[, "group1.means"])
##         group2.mean <- mean(mat[, "group2.means"])
##         meanDiff <- mean(mat[, "group1.means"] - mat[, "group2.means"])
##         tstat.sd <- mean(mat[, "tstat.sd"])
##         c(meanDiff, group1.mean, group2.mean, tstat.sd)
##     }
##     stats_ttest<- matrix(NA, ncol = 4, nrow = length(regions))
##     colnames(stats_ttest) <- c("meanDiff", "group1.mean", "group2.mean", "tstat.sd")
##     tmp <- lapply(ov.sp, getRegionStats_ttest)
##     stats_ttest[as.integer(names(tmp)),] <- do.call(rbind, tmp)
##     out <- cbind(out, as.data.frame(stats_ttest))
##     out
## }
