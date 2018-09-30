BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE){
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))

    ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
    ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")

    if(verbose) cat("[BSmooth.fstat] fitting linear models ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    allPs <- as.array(allPs)
    fit <- lmFit(allPs, design)
    fitC <- contrasts.fit(fit, contrasts)
    ## Need
    ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
    ## actuall just need
    ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    ##   rawSds <- fitC$sigma
    ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
    rawSds <- as.matrix(fitC$sigma)
    cor.coefficients <- cov2cor(fitC$cov.coefficients)
    rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    names(dimnames(rawTstats)) <- NULL
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    parameters <- c(BSseq@parameters,
                    list(design = design, contrasts = contrasts))
    stats <- list(rawSds = rawSds,
                  cor.coefficients = cor.coefficients,
                  rawTstats = rawTstats)
    out <- BSseqStat(gr = granges(BSseq),
                     stats = stats,
                     parameters = parameters)
    out
}

smoothSds <- function(BSseqStat, k = 101, qSd = 0.75, mc.cores = 1,
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
        maxGap <- BSseqStat@parameters[["maxGap"]]
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    if(verbose) cat("[smoothSds] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(granges(BSseqStat), maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    smoothSds <- do.call("c",
                         mclapply(clusterIdx, function(idx) {
                             rawSds <- getStats(BSseqStat,
                                                what = "rawSds")[idx, ]
                             rawSds <- as.array(rawSds)
                             smoothSd(rawSds, k = k, qSd = qSd)
                         }, mc.cores = mc.cores))
    smoothSds <- as.matrix(smoothSds)
    if("smoothSds" %in% names(getStats(BSseqStat)))
        BSseqStat@stats[["smoothSds"]] <- smoothSds
    else
        BSseqStat@stats <- c(getStats(BSseqStat), list(smoothSds = smoothSds))
    BSseqStat
}

# Quieten R CMD check
globalVariables("tstat")

# NOTE: Realises in memory a matrix with nrow = length(BSseqStat) and
#       ncol = length(coef)
computeStat <- function(BSseqStat, coef = NULL) {
    stopifnot(is(BSseqStat, "BSseqStat"))
    if(is.null(coef)) {
        coef <- 1:ncol(getStats(BSseqStat, what = "rawTstats"))
    }
    raw_tstats <- getStats(BSseqStat, what = "rawTstats")[, coef, drop = FALSE]
    scaled_sds <- getStats(BSseqStat, what = "rawSds") /
        getStats(BSseqStat, what = "smoothSds")
    scaled_sds_matrix <- matrix(
        rep(scaled_sds, ncol(raw_tstats)),
        ncol = ncol(raw_tstats))
    tstats <- raw_tstats * scaled_sds_matrix
    if(length(coef) > 1) {
        cor.coefficients <- getStats(BSseqStat,
                                     what = "cor.coefficients")[coef,coef]
        # NOTE: classifyTestsF() calls as.matrix(tstats) and so realises this
        #       array
        stat <- as.matrix(
            classifyTestsF(tstats, cor.coefficients, fstat.only = TRUE))
        stat.type <- "fstat"
    } else {
        stat <- tstats
        stat.type <- "tstat"
    }
    if("stat" %in% names(getStats(BSseqStat))) {
        BSseqStat@stats[["stat"]] <- stat
        BSseqStat@stats[["stat.type"]] <- stat.type
    } else {
        BSseqStat@stats <- c(getStats(BSseqStat),
                             list(stat = stat, stat.type = stat.type))
    }
    BSseqStat
}

localCorrectStat <- function(BSseqStat, threshold = c(-15,15), mc.cores = 1, verbose = TRUE) {
    compute.correction <- function(idx) {
        xx <- start(BSseqTstat)[idx]
        yy <- tstat[idx] ## FIXME
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
    maxGap <- BSseqStat$parameters$maxGap
    if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseqStat$gr, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    stat <- BSseqStat$stat
    stat.corrected <- do.call(c, mclapply(clusterIdx, compute.correction,
                                          mc.cores = mc.cores))
    BSseqTstat@stats <- cbind(getStats(BSseqTstat),
                              "tstat.corrected" = stat.corrected)
    BSseqTstat@parameters$local.local <- TRUE
    BSseqTstat
}

fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                           type = "dmrs", mc.cores = 1) {
    type <- match.arg(type, c("dmrs", "blocks"))
    bstat <- BSmooth.fstat(BSseq = BSseq, design = design,
                           contrasts = contrasts)
    bstat <- smoothSds(bstat)
    bstat <- computeStat(bstat, coef = coef)
    dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
    if (is.null(dmrs)) {
        stop("No DMRs identified. Consider reducing the 'cutoff' from (",
             paste0(cutoff, collapse = ", "), ")")
    }
    idxMatrix <- permuteAll(nperm, design)
    nullDist <- getNullDistribution_BSmooth.fstat(BSseq = BSseq,
                                                  idxMatrix = idxMatrix,
                                                  design = design,
                                                  contrasts = contrasts,
                                                  coef = coef,
                                                  cutoff = cutoff,
                                                  maxGap.sd = maxGap.sd,
                                                  maxGap.dmr = maxGap.dmr,
                                                  mc.cores = mc.cores)
    fwer <- getFWER.fstat(null = c(list(dmrs), nullDist), type = type)
    dmrs$fwer <- fwer
    meth <- getMeth(BSseq, dmrs, what = "perRegion")
    meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
    dmrs <- cbind(dmrs, meth)
    dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
    list(bstat = bstat, dmrs = dmrs, idxMatrix = idxMatrix, nullDist = nullDist)
}

fstat.comparisons.pipeline <- function(BSseq, design, contrasts, cutoff, fac,
                                       maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                                       verbose = TRUE) {
    bstat <- BSmooth.fstat(BSseq = BSseq,
                           design = design,
                           contrasts = contrasts,
                           verbose = verbose)
    bstat <- smoothSds(bstat, maxGap = maxGap.sd, verbose = verbose)
    # NOTE: Want to keep the bstat object corresponding to the original fstat
    #       and not that from the last t-tsts in the following lapply()
    bstat_f <- bstat
    if (verbose) {
        cat(paste0("[fstat.comparisons.pipeline] making ", ncol(contrasts),
                   " comparisons ... \n"))
    }
    l <- lapply(seq_len(ncol(contrasts)), function(coef) {
        if (verbose) {
            cat(paste0("[fstat.comparisons.pipeline] ",
                       colnames(contrasts)[coef], "\n"))
        }
        bstat <- computeStat(bstat, coef = coef)
        dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr,
                          verbose = verbose)
        if (!is.null(dmrs)) {
            meth <- getMeth(BSseq, dmrs, what = "perRegion")
            meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
            dmrs <- cbind(dmrs, meth)
            dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
        }
        dmrs
    })
    names(l) <- colnames(contrasts)
    list(bstat = bstat, dmrs = l)
}

