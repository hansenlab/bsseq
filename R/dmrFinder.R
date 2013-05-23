dmrFinder <- function(BSseqTstat, cutoff = NULL, qcutoff = c(0.025, 0.975),
                             maxGap = 300, column = c("tstat.corrected", "tstat"), verbose = TRUE) {
    column <- match.arg(column)
    if(! column %in% colnames(BSseqTstat@stats))
        stop("'column' needs to be a column of 'BSseqTstat@stats'")
    if(is.null(cutoff))
        cutoff <- quantile(BSseqTstat@stats[, column], qcutoff)
    direction <- as.integer(BSseqTstat@stats[, column] >= cutoff[2])
    direction[BSseqTstat@stats[, column] <= cutoff[1]] <- -1L
    direction[is.na(direction)] <- 0L
    chrs <- as.character(seqnames(BSseqTstat@gr))
    positions <- start(BSseqTstat)
    regions <- bsseq:::regionFinder3(direction, chr = chrs, pos = positions, maxGap = maxGap, verbose = verbose)
    if(verbose) cat("creating dmr data.frame\n")
    regions <- do.call(rbind, regions)
    rownames(regions) <- NULL
    stats <- getStats(BSseqTstat, regions, column = column)
    regions <- cbind(regions, stats)
    regions$direction <- ifelse(regions$meanDiff > 0, "hyper", "hypo")
    regions$width <- regions$end - regions$start + 1
    regions$invdensity <- regions$width / regions$n
    regions <- regions[order(abs(regions$areaStat), decreasing = TRUE),]
    regions
}

clusterMaker <- function(chr, pos, order.it=TRUE, maxGap=300){
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for(i in seq(along = Indexes)){
        Index <- Indexes[[i]]
        x <- pos[Index]
        if(order.it){
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}


regionFinder3 <- function(x, chr, positions, keep, maxGap = 300, verbose = TRUE) {
    getSegments3 <- function(x, cid, verbose = TRUE){
        if(verbose) cat("segmenting\n")
        segments <- cumsum(c(1, diff(x) != 0)) +
            cumsum(c(1, diff(cid) != 0))
        names(segments) <- NULL
        if(verbose) cat("splitting\n")
        out <- list(up = split(which(x > 0), segments[x > 0]),
                    down = split(which(x < 0), segments[x < 0]),
                    nochange = split(which(x == 0), segments[x == 0]))
        names(out[[1]]) <- NULL
        names(out[[2]]) <- NULL
        names(out[[3]]) <- NULL
        out
    }
    cid0 <- clusterMaker(chr, positions, maxGap = maxGap)
    segments <- getSegments3(x = x, cid = cid0, verbose = verbose)
    out <- vector("list", 2)
    names(out) <- c("up", "down")
    for(ii in 1:2) {
        idx <- segments[[ii]]
        if(length(idx) == 0) {
            out[[ii]] <- NULL
        } else {
            out[[ii]] <- data.frame(chr = sapply(idx, function(jj) chr[jj[1]]),
                                    start = sapply(idx, function(jj) min(positions[jj])),
                                    end = sapply(idx, function(jj) max(positions[jj])),
                                    idxStart = sapply(idx, function(jj) min(jj)),
                                    idxEnd = sapply(idx, function(jj) max(jj)),
                                    cluster = sapply(idx, function(jj) cid0[jj[1]]),
                                    n = sapply(idx, length))
        }
    }
    out
}

