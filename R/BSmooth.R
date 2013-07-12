makeClusters <- function(hasGRanges, maxGap = 10^8, mc.cores = 1) {
    chrOrder <- as.character(runValue(seqnames(hasGRanges)))
    if(anyDuplicated(chrOrder))
        stop("argument 'hasGRanges' is not properly order")
    grBase <- granges(hasGRanges)
    clusters <- reduce(resize(grBase, width = 2*maxGap + 1, fix = "center"))
    start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
    clusters.sp <- split(clusters, seqnames(clusters))
    stopifnot(all(sapply(clusters.sp, function(cluster.gr) {
        if(length(cluster.gr) <= 1) return(TRUE)
        all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
    }))) # are the clusters ordered within the chromosome? This is probably guranteed
    clusters <- Reduce(c, clusters.sp[chrOrder])
    stopifnot(all(chrOrder == runValue(seqnames(clusters))))
    ov <- bsseq:::findOverlaps_mclapply(grBase, clusters, mc.cores = mc.cores)
    clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
    names(clusterIdx) <- NULL
    clusterIdx
}
    

BSmooth <- function(BSseq, ns = 70, h = 1000, maxGap = 10^8, parallelBy = c("sample", "chromosome"),
                    mc.preschedule = FALSE, mc.cores = 1, keep.se = FALSE, verbose = TRUE) {
    smooth <- function(idxes, sampleIdx) {
        ## Assuming that idxes is a set of indexes into the BSseq object
        ## sampleIdx is a single character
        this_sample_chr <- c(sampleNames(BSseq)[sampleIdx],
                             as.character(seqnames(BSseq)[idxes[1]]))
        if(verbose >= 2)
            cat(sprintf("[BSmooth]   smoothing start: sample:%s, chr:%s, nLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes)))
        Cov <- getCoverage(BSseq, type = "Cov")[idxes, sampleIdx]
        M <- getCoverage(BSseq, type = "M")[idxes, sampleIdx]
        pos <- start(BSseq)[idxes]
        stopifnot(all(diff(pos) > 0))
        wh <- which(Cov != 0)
        nn <- ns / length(wh)
        if(length(wh) <= ns) {
            if(keep.se)
                se.coef <- rep(NA_real_, length(Cov))
            else
                se.coef <- NULL
            return(list(coef = rep(NA_real_, length(Cov)),
                        se.coef = se.coef,
                        trans = NULL, h = h, nn = nn))
        }
        sdata <- data.frame(pos = pos[wh],
                            M = pmin(pmax(M[wh], 0.01), Cov[wh] - 0.01),
                            Cov = Cov[wh])
        fit <- locfit(M ~ lp(pos, nn = nn, h = h), data = sdata,
                      weights = Cov, family = "binomial", maxk = 10000)
        pp <- preplot(fit, where = "data", band = "local",
                      newdata = data.frame(pos = pos))
        if(keep.se) {
            se.coef <- pp$se.fit
        } else {
            se.coef <- NULL
        }
        if(verbose >= 2)
            cat(sprintf("[BSmooth]   smoothing end: sample:%s, chr:%s, nLoci:%s, nCoveredLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes), nrow(sdata)))
        return(list(coef = pp$fit, se.coef = se.coef,
                    trans = pp$trans, h = h, nn = nn))
    }
    stopifnot(class(BSseq) == "BSseq")
    parallelBy <- match.arg(parallelBy)
    if(verbose) cat("[BSmooth] preprocessing ... ")
    stime <- system.time({
        clusterIdx <- makeClusters(BSseq, maxGap = maxGap, mc.cores = mc.cores)
    })[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames

    stimeAll <- system.time({
        switch(parallelBy, "sample" = {
            if(verbose) cat(sprintf("[BSmooth] smoothing by 'sample' (mc.cores = %d, mc.preschedule = %s)\n",
                                    mc.cores, mc.preschedule))
            out <- mclapply(seq(along = sampleNames), function(sIdx) {
                stime <- system.time({
                    tmp <- lapply(clusterIdx, function(jj) {
                        try(smooth(idxes = jj, sampleIdx = sIdx))
                    })
                    coef <- do.call(c, lapply(tmp, function(xx) xx$coef))
                    se.coef <- do.call(c, lapply(tmp, function(xx) xx$se.coef))
                })[3]
                if(verbose) {
                    cat(sprintf("[BSmooth] sample %s (out of %d), done in %.1f sec\n",
                                sampleNames[sIdx], length(sampleNames), stime))
                }
                return(list(coef = coef, se.coef = se.coef))
            }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
            if(any(sapply(out, is, class2 = "try-error")))
                stop("BSmooth encountered smoothing errors")
            coef <- do.call(cbind, lapply(out, function(xx) xx$coef))
            se.coef <- do.call(cbind, lapply(out, function(xx) xx$se.coef))
        }, "chromosome" = {
            if(verbose) cat(sprintf("[BSmooth] smoothing by 'chromosome' (mc.cores = %d, mc.preschedule = %s)\n",
                                    mc.cores, mc.preschedule))
            out <- mclapply(1:length(clusterIdx), function(ii) {
                stime <- system.time({
                    tmp <- lapply(seq(along = sampleNames), function(sIdx) {
                        smooth(idxes = clusterIdx[[ii]], sampleIdx = sIdx)
                    })
                    coef <- do.call(cbind, lapply(tmp, function(xx) xx$coef))
                    se.coef <- do.call(cbind, lapply(tmp, function(xx) xx$se.coef))
                })[3]
                if(verbose)
                    cat(sprintf("[BSmooth] chr idx %d (out of %d), done in %.1f sec\n",
                                ii, length(clusterIdx), stime))
                return(list(coef = coef, se.coef = se.coef))
            }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
            if(any(sapply(out, is, class2 = "try-error")))
                stop("BSmooth encountered smoothing errors")
            coef <- do.call(rbind, lapply(out, function(xx) xx$coef))
            se.coef <- do.call(rbind, lapply(out, function(xx) xx$se.coef))
        })
    })[3]
    if(verbose)
        cat(sprintf("[BSmooth] smoothing done in %.1f sec\n", stimeAll))

    rownames(coef) <- NULL
    colnames(coef) <- sampleNames(BSseq)
    if(!is.null(se.coef)) {
        rownames(se.coef) <- NULL
        colnames(se.coef) <- sampleNames(BSseq)
    }

    if(!is.null(coef))
        assay(BSseq, "coef") <- coef
    if(!is.null(se.coef))
        assay(BSseq, "se.coef") <- se.coef
    mytrans <- function(x) {
        y <- x
        ix <- which(x < 0)
        ix2 <- which(x > 0)
        y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
        y[ix2] <- 1/(1 + exp(-x[ix2]))
        y
    }
    environment(mytrans) <- baseenv()
    BSseq@trans <- mytrans
    parameters <- list(smoothText = sprintf("BSmooth (ns = %d, h = %d, maxGap = %d)", ns, h, maxGap),
                       ns = ns, h = h, maxGap = maxGap)
    BSseq@parameters <- parameters
    BSseq
}


