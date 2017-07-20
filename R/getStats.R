getStats <- function(bstat, regions = NULL, ...) {
    stopifnot(is(bstat, "BSseqTstat") || is(bstat, "BSseqStat"))
    if(is(bstat, "BSseqTstat"))
        return(getStats_BSseqTstat(bstat, regions = regions, ...))
    getStats_BSseqStat(bstat, regions = regions, ...)
}

# NOTE: Realises in memory a matrix with nrow = length(hits) and ncol = 1
.getRegionStats <- function(stat, hits, na.rm = FALSE) {
    stat_by_region <- split(as.array(stat[queryHits(hits), ]),
                            subjectHits(hits))
    data.frame(areaStat = vapply(stat_by_region, sum, numeric(1),
                                 na.rm = na.rm, USE.NAMES = FALSE),
               maxStat = vapply(stat_by_region, max, numeric(1),
                                na.rm = na.rm, USE.NAMES = FALSE))
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
    hits <- findOverlaps(BSseqStat, regions)
    regionStats <- .getRegionStats(stat = BSseqStat@stats$stat,
                                   hits = hits)
    regionStats
}

# NOTE: Realises in memory a matrix with nrow = length(hits) and ncol = 3
.getRegionStats_ttest <- function(stats, hits, na.rm = FALSE) {
    stat_names <- c("group1.means", "group2.means", "tstat.sd")
    stats <- as.array(stats[queryHits(hits), stat_names])
    stats_by_region <- lapply(seq_along(stat_names), function(j) {
        split(stats[, j], subjectHits(hits))
    })
    names(stats_by_region) <- stat_names
    meanDiff <- mapply(function(g1, g2, na.rm) {
        mean(g1 - g2, na.rm = na.rm)
    }, g1 = stats_by_region[["group1.means"]],
    g2 = stats_by_region[["group2.means"]],
    MoreArgs = list(na.rm = na.rm), SIMPLIFY = TRUE, USE.NAMES = FALSE)
    group1.mean <- vapply(stats_by_region[["group1.means"]], mean, numeric(1),
                          na.rm = na.rm, USE.NAMES = FALSE)
    group2.mean <- vapply(stats_by_region[["group2.means"]], mean, numeric(1),
                          na.rm = na.rm, USE.NAMES = FALSE)
    tstat.sd <- vapply(stats_by_region[["tstat.sd"]], mean, numeric(1),
                       na.rm = na.rm, USE.NAMES = FALSE)
    data.frame(meanDiff = meanDiff,
               group1.mean = group1.mean,
               group2.mean = group2.mean,
               tstat.sd = tstat.sd)
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
    hits <- findOverlaps(BSseqTstat, regions)
    out <- .getRegionStats(stat = BSseqTstat@stats[, stat, drop = FALSE],
                           hits = hits)
    if(! stat %in% c("tstat.corrected", "tstat"))
        return(out)
    stats_ttest <- .getRegionStats_ttest(stats = BSseqTstat@stats,
                                         hits = hits)
    out <- cbind(out, stats_ttest)
    out
}
