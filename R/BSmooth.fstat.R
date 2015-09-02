BSmooth.fstat <- function(BSseq, groups, designMat, contrasts,
                          weights = NULL, maxGap = NULL, qSd = 0.75,
                          k = 101, mc.cores = 1, verbose = TRUE){
    smoothSd <- function(Sds, k) {
        k0 <- floor(k/2)
        if(all(is.na(Sds))) return(Sds)
        thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
        addSD <- rep(median(Sds, na.rm = TRUE), k0)
        sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
        sSds
    }
    
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))
    
    ## for(i in 1:length(groups)){
    ##     if(is.character(groups[[i]])){
    ##         stopifnot(all(groups[[i]] %in% sampleNames(BSseq)))
    ##         groups[[i]] <- match(groups[[i]], sampleNames(BSseq))
    ##     }
    ##     if(is.numeric(groups[[i]])){
    ##         group <- groups[[i]]
    ##         stopifnot(min(group) >=1 & max(group) <= ncol(BSseq))
    ##     }
    ##     else{
    ##         stop(paste("Problems with argument group", i))
    ##     }
        
    ## }
    
    ## for(i in 1: (length(groups) - 1)){
    ##     for(j in (i + 1):length(groups)){
    ##         stopifnot(length(intersect(groups[[i]], groups[[j]])) == 0)
    ##     }
    ## }
    
    ## for(i in 1:length(groups)){
    ##     stopifnot(length(groups[[i]]) > 0)
    ## }
    
    ## n_i <- unlist((lapply(groups, function(x) length(x))))
    ## N <- sum(n_i)
    ## stopifnot(N >= length(groups) + 1)
    
    
    ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
    ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
    
    if(is.null(maxGap))
        maxGap <- BSseq@parameters$maxGap
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    
    if(verbose) cat("[BSmooth.fstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(verbose) cat("[BSmooth.fstat] computing variation between groups ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    
    contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
    fit <- lmFit(allPs, design, weights = weights)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    rawSds <- fit2$sigma
    smoothSds <- do.call("c", mclapply(clusterIdx, function(idx) {
                                  smoothSd(rawSds[idx], k = k)
                              }, mc.cores = mc.cores))
    MS_w <- smoothSds^2
    ## Following is tstat with no smoothing or local correction
    ## rawTstat <- fit2$coefficients / fit2$stdev.unscaled / fit2$sigma
    ## limma:classifyTestsF has code for representing an F stat as a transform of t-stats
    ## This is really an fstat
    fstat <- fit2$coefficients / fit2$stdev.unscaled / smoothSds
    
    
    ## fit2b <- eBayes(fit2, trend = T)
    ## MS_b <- fit2b$F * fit2b$s2.post
    ## fstat <- MS_b/MS_w
    is.na(fstat)[MS_w == 0] <- TRUE

    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    stats <- cbind(rawSds, MS_w, fstat)
    colnames(stats) <- c("rawSds", "smoothSdsSquared", "fstat")
    parameters <- c(BSseq@parameters,
                    list(fstatText = sprintf("BSmooth.fstat (maxGap = %d)", maxGap),
                         k = k, qSd = qSd, maxGap = maxGap))
    out <- BSseqFstat(gr = granges(BSseq), stats = stats, parameters = parameters)
    out
}




setClass("BSseqFstat", contains = "hasGRanges", 
         representation(stats = "matrix",
                        parameters = "list")
         )
setValidity("BSseqFstat", function(object) {
    msg <- NULL
    if(length(object@gr) != nrow(object@stats))
        msg <- c(msg, "length of 'gr' is different from the number of rows of 'stats'")
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqFstat"),
          function(object) {
    cat("An object of type 'BSseqFstat' with\n")
    cat(" ", length(object), "methylation loci\n")
    cat("based on smoothed data:\n")
    cat(" ", object@parameters$smoothText, "\n")
    cat("with parameters\n")
    cat(" ", object@parameters$FstatText, "\n")
})

setMethod("[", "BSseqFstat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x@stats <- x@stats[i,, drop = FALSE]
    x
})

BSseqFstat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqFstat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}

summary.BSseqFstat <- function(object, ...) {
    quant <- quantile(geFstats(object)[, "Fstat.corrected"],
                      prob = c(0.0001, 0.001, 0.01, 0.5, 0.99, 0.999, 0.9999))
    quant <- t(t(quant))
    colnames(quant) <- "quantiles"
    out <- list(quantiles = quant)
    class(out) <- "summary.BSseqFstat"
    out
}

print.summary.BSseqFstat <- function(x, ...) {
    print(as.matrix(x$quantiles))
}

plot.BSseqFstat <- function(x, y, ...) {
    Fstat <- geFstats(x)[, "Fstat"]
    plot(density(Fstat), xlim = c(-10,10), col = "blue", main = "")
    if("Fstat.corrected" %in% colnames(geFstats(x))) {
        Fstat.cor <- geFstats(x)[, "Fstat.corrected"]
        lines(density(Fstat.cor), col = "black")
        legend("topleft", legend = c("uncorrected", "corrected"), lty = c(1,1),
               col = c("blue", "black"))
    } else {
        legend("topleft", legend = c("uncorrected"), lty = 1,
               col = c("blue"))
    }
}

geFstats <- function(BSseqFstat, regions = NULL, stat = "Fstat.corrected") {
    stopifnot(is(BSseqFstat, "BSseqFstat"))
    if(is.null(regions))
        return(BSseqFstat@stats)
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(stat %in% colnames(BSseqFstat@stats))
    stopifnot(length(stat) == 1)
    stopifnot(is(regions, "GenomicRanges"))
    ov <- findOverlaps(BSseqFstat, regions)
    ov.sp <- split(queryHits(ov), subjectHits(ov))
    getRegionStats <- function(idx) {
        mat <- BSseqFstat@stats[idx,, drop=FALSE]
        areaStat <- sum(mat[, stat])
        maxStat <- max(mat[, stat])
        c(areaStat, maxStat)
    }
    stats <- matrix(NA, ncol = 2, nrow = length(regions))
    colnames(stats) <- c("areaStat", "maxStat")
    tmp <- lapply(ov.sp, getRegionStats)
    stats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    out <- as.data.frame(stats)
    if(! stat %in% c("Fstat.corrected", "Fstat"))
        return(out)
    getRegionStats_ttest <- function(idx) {
        mat <- BSseqFstat@stats[idx,, drop=FALSE]
        group1.mean <- mean(mat[, "group1.means"])
        group2.mean <- mean(mat[, "group2.means"])
        meanDiff <- mean(mat[, "group1.means"] - mat[, "group2.means"])
        Fstat.sd <- mean(mat[, "Fstat.sd"])
        c(meanDiff, group1.mean, group2.mean, Fstat.sd)
    }
    stats_ttest<- matrix(NA, ncol = 4, nrow = length(regions))
    colnames(stats_ttest) <- c("meanDiff", "group1.mean", "group2.mean", "Fstat.sd")
    tmp <- lapply(ov.sp, getRegionStats_ttest)
    stats_ttest[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    out <- cbind(out, as.data.frame(stats_ttest))
    out
}
