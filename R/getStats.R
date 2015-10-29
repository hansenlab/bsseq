getStats <- function(bstat, regions = NULL, ...) {
    stopifnot(is(bstat, "BSseqTstat") || is(bstat, "BSseqStat"))
    if(is(bstat, "BSseqTstat"))
        return(getStats_BSseqTstat(bstat, regions = regions, ...))
    getStats_BSseqStat(bstat, regions = regions, ...)
}

getStats_BSseqStat <- function(BSseqStat, regions = NULL, what = NULL) {
    stopifnot(is(BSseqStat, "BSseqStat"))
    if(!is.null(what)) {
        stopifnot(what %in% names(BSseqStat@stats))
        return(BSseqStat@stats[[what]])
    }
    if(is.null(regions))
        return(BSseqStat@stats)
    ## Now we have regions and no what
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    ov <- findOverlaps(BSseqStat, regions)
    ov.sp <- split(queryHits(ov), subjectHits(ov))
    ## We need areaStat and maxStat
    ## Could need averageMeth in each group?
    ## Need to have a specific design for that
    ## Could at least get contrast coefficient
    getRegionStats <- function(idx) {
        areaStat <- sum(BSseqStat@stats$stat[idx], na.rm = TRUE)
        maxStat <- max(BSseqStat@stats$stat[idx], na.rm = TRUE)
        c(areaStat, maxStat)
    }
    regionStats <- matrix(NA, ncol = 2, nrow = length(regions))
    colnames(regionStats) <- c("areaStat", "maxStat")
    tmp <- lapply(ov.sp, getRegionStats)
    regionStats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    regionStats
}

getStats_BSseqTstat <- function(BSseqTstat, regions = NULL, stat = "tstat.corrected") {
    stopifnot(is(BSseqTstat, "BSseqTstat"))
    if(is.null(regions))
        return(BSseqTstat@stats)
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(stat %in% colnames(BSseqTstat@stats))
    stopifnot(length(stat) == 1)
    stopifnot(is(regions, "GenomicRanges"))
    ov <- findOverlaps(BSseqTstat, regions)
    ov.sp <- split(queryHits(ov), subjectHits(ov))
    getRegionStats <- function(idx) {
        mat <- BSseqTstat@stats[idx,, drop=FALSE]
        areaStat <- sum(mat[, stat])
        maxStat <- max(mat[, stat])
        c(areaStat, maxStat)
    }
    stats <- matrix(NA, ncol = 2, nrow = length(regions))
    colnames(stats) <- c("areaStat", "maxStat")
    tmp <- lapply(ov.sp, getRegionStats)
    stats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    out <- as.data.frame(stats)
    if(! stat %in% c("tstat.corrected", "tstat"))
        return(out)
    getRegionStats_ttest <- function(idx) {
        mat <- BSseqTstat@stats[idx,, drop=FALSE]
        group1.mean <- mean(mat[, "group1.means"])
        group2.mean <- mean(mat[, "group2.means"])
        meanDiff <- mean(mat[, "group1.means"] - mat[, "group2.means"])
        tstat.sd <- mean(mat[, "tstat.sd"])
        c(meanDiff, group1.mean, group2.mean, tstat.sd)
    }
    stats_ttest<- matrix(NA, ncol = 4, nrow = length(regions))
    colnames(stats_ttest) <- c("meanDiff", "group1.mean", "group2.mean", "tstat.sd")
    tmp <- lapply(ov.sp, getRegionStats_ttest)
    stats_ttest[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    out <- cbind(out, as.data.frame(stats_ttest))
    out
}
