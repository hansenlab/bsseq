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
    ov <- findOverlaps_mclapply(grBase, clusters, mc.cores = mc.cores)
    clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
    names(clusterIdx) <- NULL
    clusterIdx
}
    

BSmooth <- function(BSseq, ns = 70, h = 1000, maxGap = 10^8, parallelBy = c("sample", "chromosome"),
                    mc.preschedule = FALSE, mc.cores = 1, keep.se = FALSE, verbose = TRUE) {
    smooth <- function(idxes, sname) {
        ## Assuming that idxes is a set of indexes into the BSseq object
        ## sname is a single character
        if(verbose >= 3)
            cat(sprintf("    beginning: sample:%s, chr:%s, nLoci:%s\n",
                        sname, as.character(seqnames(BSseq)[idxes[1]]),
                        length(idxes)))
        Cov <- getCoverage(BSseq, type = "Cov")[idxes, sname]
        M <- getCoverage(BSseq, type = "M")[idxes, sname]
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
            cat(sprintf("    sample:%s, chr:%s, nLoci:%s, nCoveredLoci:%s, done\n",
                        sname, as.character(seqnames(BSseq)[idxes[1]]),
                        length(idxes), nrow(sdata)))
        return(list(coef = pp$fit, se.coef = se.coef,
                    trans = pp$trans, h = h, nn = nn))
    }
    stopifnot(class(BSseq) == "BSseq")
    parallelBy <- match.arg(parallelBy)
    if(verbose) cat("preprocessing ... ")
    stime <- system.time({
        clusterIdx <- makeClusters(BSseq, maxGap = maxGap, mc.cores = mc.cores)
    })[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames

    stimeAll <- system.time({
        switch(parallelBy, "sample" = {
            if(verbose) cat(sprintf("smoothing by 'sample' (mc.cores = %d, mc.preschedule = %s)\n",
                                    mc.cores, mc.preschedule))
            out <- mclapply(sampleNames, function(nam) {
                stime <- system.time({
                    tmp <- lapply(clusterIdx, function(jj) {
                        smooth(idxes = jj, sname = nam)
                    })
                    coef <- do.call(c, lapply(tmp, function(xx) xx$coef))
                    se.coef <- do.call(c, lapply(tmp, function(xx) xx$se.coef))
                })[3]
                if(verbose) {
                    cat(sprintf("  sample %s (out of %d), done in %.1f sec\n",
                                nam, length(sampleNames), stime))
                }
                return(list(coef = coef, se.coef = se.coef))
            }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
            coef <- do.call(cbind, lapply(out, function(xx) xx$coef))
            se.coef <- do.call(cbind, lapply(out, function(xx) xx$se.coef))
        }, "chromosome" = {
            if(verbose) cat(sprintf("smoothing by 'chromosome' (mc.cores = %d, mc.preschedule = %s)\n",
                                    mc.cores, mc.preschedule))
            out <- mclapply(1:length(clusterIdx), function(ii) {
                stime <- system.time({
                    tmp <- lapply(sampleNames, function(nam) {
                        smooth(idxes = clusterIdx[[ii]], sname = nam)
                    })
                    coef <- do.call(cbind, lapply(tmp, function(xx) xx$coef))
                    se.coef <- do.call(cbind, lapply(tmp, function(xx) xx$se.coef))
                })[3]
                if(verbose)
                    cat(sprintf("  chr idx %d (out of %d), done in %.1f sec\n",
                                ii, length(clusterIdx), stime))
                return(list(coef = coef, se.coef = se.coef))
            }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
            coef <- do.call(rbind, lapply(out, function(xx) xx$coef))
            se.coef <- do.call(rbind, lapply(out, function(xx) xx$se.coef))
        })
    })[3]
    if(verbose)
        cat(sprintf("smoothing done in %.1f sec\n", stimeAll))

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


BSmooth.tstat <- function(BSseq, group1, group2, estimate.var = c("same", "paired", "group2"),
                          local.correct = TRUE, maxGap = NULL, qSd = 0.75, k = 101, mc.cores = 1, verbose = TRUE){
    smoothSd <- function(Sds, k) {
        k0 <- floor(k/2)
        if(all(is.na(Sds))) return(Sds)
        thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
        addSD <- rep(median(Sds, na.rm = TRUE), k0)
        sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
        sSds
    }
    compute.correction <- function(idx, qSd = 0.75) {
        xx <- start(BSseq)[idx]
        yy <- tstat[idx]
        suppressWarnings({
            drange <- diff(range(xx, na.rm = TRUE))
        })
        if(drange <= 25000)
            return(yy)
        tstat.function <- approxfun(xx, yy)
        xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
        yy.reg <- tstat.function(xx.reg)
        fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0),
                      family = "huber", maxk = 10000) 
        correction <- predict(fit, newdata = data.frame(xx.reg = xx))
        yy - correction 
    }

    estimate.var <- match.arg(estimate.var)
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))
    if(is.numeric(group1)) {
        stopifnot(min(group1) >= 1 & max(group1) <= ncol(BSseq))
        group1 <- sampleNames(BSseq)[group1]
    }
    if(is.numeric(group2)) {
        stopifnot(min(group2) >= 1 & max(group2) <= ncol(BSseq))
        group2 <- sampleNames(BSseq)[group2]
    }
    stopifnot(length(intersect(group1, group2)) == 0)
    if(estimate.var == "paired")
        stopifnot(length(group1) == length(group2))
    stopifnot(all(c(group1, group2) %in% sampleNames(BSseq)))

    if(any(rowSums(getCoverage(BSseq)[, c(group1, group2)]) == 0))
        warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
    
    if(is.null(maxGap))
        maxGap <- BSseq@parameters$maxGap
    
    if(verbose) cat("preprocessing ... ")
    stime <- system.time({
        clusterIdx <- makeClusters(BSseq, maxGap = maxGap, mc.cores = mc.cores)
    })[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
        
    if(verbose) cat("computing stats within groups ... ")
    stime <- system.time({
        allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                         confint = FALSE)[, c(group1, group2)]
        group1.means <- rowMeans(allPs[, group1, drop = FALSE], na.rm = TRUE)
        group2.means <- rowMeans(allPs[, group2, drop = FALSE], na.rm = TRUE)
        })[3]                               
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(verbose) cat("computing stats across groups ... ")
    stime <- system.time({
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
        if(local.correct) {
            tstat.corrected <- do.call(c, mclapply(clusterIdx,
                                                   compute.correction, qSd = qSd,
                                                   mc.cores = mc.cores))
        }
    })[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(local.correct) {
        stats <- cbind(rawSds, tstat.sd, group2.means, group1.means,
                       tstat, tstat.corrected)
        colnames(stats) <- c("rawSds", "tstat.sd",
                             "group2.means", "group1.means", "tstat",
                             "tstat.corrected")
 
    } else {
        stats <- cbind(rawSds, tstat.sd, group2.means, group1.means,
                       tstat)
        colnames(stats) <- c("rawSds", "tstat.sd",
                             "group2.means", "group1.means", "tstat")
    }
    
    parameters <- c(BSseq@parameters,
                    list(tstatText = sprintf("BSmooth.tstat (local.correct = %s, maxGap = %d)",
                         local.correct, maxGap),
                         group1 = group1, group2 = group2, k = k, qSd = qSd,
                         local.correct = local.correct, maxGap = maxGap))
    out <- BSseqTstat(gr = granges(BSseq), stats = stats, parameters = parameters)
    out
}
