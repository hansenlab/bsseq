BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE){
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))
        
    ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
    ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
    
    if(verbose) cat("[BSmooth.fstat] fitting linear models ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    fit <- lmFit(allPs, design)
    fitC <- contrasts.fit(fit, contrasts)
    ## Need
    ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
    ## actuall just need
    ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    ##   rawSds <- fitC$sigma
    ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
    ## might also need means to get at mean vs variation relationship
    rawSds <- fitC$sigma
    cor.coefficients <- cov2cor(fitC$cov.coefficients)
    rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    names(dimnames(rawTstats)) <- NULL
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    parameters <- c(BSseq@parameters,
                    list(design = design, contrasts = contrasts))    
    out <- list(gr = granges(BSseq),
                rawSds = rawSds,
                cor.coefficients = cor.coefficients,
                rawTstats = rawTstats,
                parameters = parameters)
    out
}

extractStat <- function(fstat.list, coef = NULL) {
    if(is.null(coef)) {
        coef <- 1:ncol(fstat.list$rawTstats)
    }
    tstats <- fstat.list$rawTstats[, coef, drop = FALSE]
    tstats <- tstats * fstat.list$rawSds / fstat.list$smoothSds
    if(length(coef) > 1) {
        cor.coefficients <- fstat.list$cor.coefficients[coef,coef]
        stat <- classifyTestsF(tstats, cor.coefficients, fstat.only = TRUE)
    }
    else
        stat <- tstats
    fstat.list$stat <- stat
    fstat.list
}

localCorrectStat <- function(fstat.list, threshold = c(-15,15), mc.cores = 1, verbose = TRUE) {
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
    maxGap <- fstat.list$parameters$maxGap
    if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- bsseq:::makeClusters(fstat.list$gr, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    stat <- fstat.list$stat
    stat.corrected <- do.call(c, mclapply(clusterIdx, compute.correction, mc.cores = mc.cores))
    BSseqTstat@stats <- cbind(BSseqTstat@stats, "tstat.corrected" = tstat.corrected)
    BSseqTstat@parameters$local.local <- TRUE
    BSseqTstat
}




smoothSds <- function(fstat.list, k = 101, qSd = 0.75, mc.cores = 1,
                      maxGap = 10^8, verbose = TRUE) {
    smoothSd <- function(Sds, k, qSd) {
        k0 <- floor(k/2)
        if(all(is.na(Sds))) return(Sds)
        thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
        addSD <- rep(median(Sds, na.rm = TRUE), k0)
        sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
        sSds
    }
    if(is.null(maxGap))
        maxGap <- fstat.list$parameters$maxGap
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    if(verbose) cat("[smoothSds] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- bsseq:::makeClusters(fstat.list$gr, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    smoothSds <- do.call("c", mclapply(clusterIdx, function(idx) {
                                  smoothSd(fstat.list$rawSds[idx], k = k, qSd = qSd)
                              }, mc.cores = mc.cores))
    fstat.list$smoothSds <- smoothSds
    fstat.list
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
