BSmooth.tstat2 <- function(BSseq, group1, group2, estimate.var = c("same", "paired", "group2"),
                           maxGap = NULL, qSd = 0.75, k = 101, mc.cores = 1, verbose = TRUE){
    smoothSd <- function(Sds, k) {
        k0 <- floor(k/2)
        if(all(is.na(Sds))) return(Sds)
        thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
        addSD <- rep(median(Sds, na.rm = TRUE), k0)
        sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
        sSds
    }
    estimate.var <- match.arg(estimate.var)
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))
    if(is.character(group1)) {
        stopifnot(all(group1 %in% sampleNames(BSseq)))
        group1 <- match(group1, sampleNames(BSseq))
    }
    if(is.numeric(group1)) {
        stopifnot(min(group1) >= 1 & max(group1) <= ncol(BSseq))
    } else stop("problems with argument 'group1'")
    if(is.character(group2)) {
        stopifnot(all(group2 %in% sampleNames(BSseq)))
        group2 <- match(group2, sampleNames(BSseq))
    }    
    if(is.numeric(group2)) {
        stopifnot(min(group2) >= 1 & max(group2) <= ncol(BSseq))
    } else stop("problems with argument 'group2'")
    stopifnot(length(intersect(group1, group2)) == 0)
    stopifnot(length(group1) > 0)
    stopifnot(length(group2) > 0)
    stopifnot(length(group1) + length(group2) >= 3)
    if(estimate.var == "paired")
        stopifnot(length(group1) == length(group2))
    
    if(any(rowSums(getCoverage(BSseq)[, c(group1, group2)]) == 0))
        warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
    
    if(is.null(maxGap))
        maxGap <- BSseq@parameters$maxGap
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    
    if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
        
    if(verbose) cat("[BSmooth.tstat] computing stats within groups ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    group1.means <- rowMeans(allPs[, group1, drop = FALSE], na.rm = TRUE)
    group2.means <- rowMeans(allPs[, group2, drop = FALSE], na.rm = TRUE)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(verbose) cat("[BSmooth.tstat] computing stats across groups ... ")
    ptime1 <- proc.time()
    switch(estimate.var,
           "group2" = {
               rawSds <- rowSds(allPs[, group2, drop = FALSE], na.rm = TRUE)
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1) + 1/length(group2))
               tstat.sd <- smoothSds * scale
           },
           "same" = {
               rawSds <- sqrt( ((length(group1) - 1) * rowVars(allPs[, group1, drop = FALSE]) +
                                (length(group2) - 1) * rowVars(allPs[, group2, drop = FALSE])) /
                              (length(group1) + length(group2) - 2))
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1) + 1/length(group2))
               tstat.sd <- smoothSds * scale
           },
           "paired" = {
               rawSds <- rowSds(allPs[, group1, drop = FALSE] - allPs[, group2, drop = FALSE])
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1))
               tstat.sd <- smoothSds * scale
           })
    tstat <- (group1.means - group2.means) / tstat.sd
    is.na(tstat)[tstat.sd == 0] <- TRUE
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    stats <- cbind(rawSds, tstat.sd, group2.means, group1.means, tstat)
    colnames(stats) <- c("rawSds", "tstat.sd", "group2.means", "group1.means", "tstat")
    parameters <- c(BSseq@parameters,
                    list(tstatText = sprintf("BSmooth.tstat2 (local.correct = %s, maxGap = %d)",
                                             FALSE, maxGap),
                         group1 = group1, group2 = group2, k = k, qSd = qSd,
                         local.correct = FALSE, maxGap = maxGap))
    out <- granges(BSseq)
    out <- BSseqTstat(gr = granges(BSseq), stats = stats, parameters = parameters)
    out
}


tstat.localCorrect <- function(BSseqTstat, threshold = c(-15,15), mc.cores = 1, verbose = TRUE) {
    compute.correction <- function(idx) {
        xx <- start(BSseqTstat)[idx]
        yy <- tstat[idx]
        if(!is.null(threshold)) {
            stopifnot(is.numeric(threshold) && length(threshold) == 2)
            stopifnot(threshold[1] < 0 && threshold[2] > 0)
            yy[yy < threshold[1]] <- threshold[1]
            yy[yy > threshold[2]] <- threshold[2]
        }
        suppressWarnings({
            drange <- diff(range(xx, na.rm = TRUE))
        })
        if(drange <= 25000)
            return(yy)
        tstat.function <- approxfun(xx, yy)
        xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
        yy.reg <- tstat.function(xx.reg)
        fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0),
                      family = "huber", maxk = 50000) 
        correction <- predict(fit, newdata = data.frame(xx.reg = xx))
        yy - correction 
    }
    maxGap <- BSseqTstat@parameters$maxGap
    if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- bsseq:::makeClusters(BSseqTstat, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    tstat <- BSseqTstat@stats[, "tstat"]
    tstat.corrected <- do.call(c, mclapply(clusterIdx, compute.correction, mc.cores = mc.cores))
    BSseqTstat@stats <- cbind(BSseqTstat@stats, "tstat.corrected" = tstat.corrected)
    BSseqTstat@parameters$local.local <- TRUE
    BSseqTstat
}


