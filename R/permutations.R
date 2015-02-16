getNullPermutation <- function(BSseq, idxMatrix1, idxMatrix2,
                               estimate.var, local.correct,
                               cutoff, stat, maxGap, mc.cores = 1) {
    stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
    nullDist <- mclapply(1:nrow(idxMatrix1), function(ii) {
        BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var,
                                  group1 = idxMatrix1[ii,],
                                  group2 = idxMatrix2[ii,],
                                  local.correct = local.correct, maxGap = 10^8)
        dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = maxGap)
        dmrs0
    }, mc.cores = min(nrow(idxMatrix1), mc.cores), mc.preschedule = FALSE)
    nullDist
}

getNullBlocks <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                          mc.cores = 1) {
    getNullPermutation(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = FALSE,
                       estimate.var = estimate.var,
                       cutoff = c(-2,2), stat = "tstat", maxGap = 10000,
                       mc.cores = mc.cores)
}

getNullDmrs <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                        mc.cores = 1) {
    getNullPermutation(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = TRUE,
                       estimate.var = estimate.var,
                       cutoff = c(-4.6,4.6), stat = "tstat.corrected", maxGap = 300,
                       mc.cores = mc.cores)
}

subsetDmrs <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    subset(xx, n >= 3 & abs(meanDiff) > 0.1 & invdensity <= 300)
}

subsetBlocks <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    subset(xx, width >= 10000)
}

getFWER <- function(null, type = "blocks") {
    reference <- null[[1]]
    null <- null[-1]
    better <- sapply(1:nrow(reference), function(ii) {
        meanDiff <- reference$meanDiff[ii]
        width <- reference$width[ii]
        n <- reference$n[ii]
        if(meanDiff < 0 && type == "blocks") {
            out <- sapply(null, function(nulldist) {
                any(nulldist$meanDiff <= meanDiff &
                    nulldist$width >= width)
            })
        }
        if(meanDiff > 0 && type == "blocks") {
            out <- sapply(null, function(nulldist) {
                any(nulldist$meanDiff <= meanDiff &
                    nulldist$width >= width)
            })
        }
        if(meanDiff < 0 && type == "dmrs") {
            out <- sapply(null, function(nulldist) {
                any(nulldist$meanDiff <= meanDiff &
                    nulldist$n >= n)
            })
        }
        if(meanDiff > 0 && type == "dmrs") {
            out <- sapply(null, function(nulldist) {
                any(nulldist$meanDiff <= meanDiff &
                    nulldist$n >= n)
            })
        }
        sum(out)
    })
    better
}
       
makeIdxMatrix <- function(group1, group2) {
    stopifnot(length(group2) >= length(group1))
    stopifnot(length(group1) > 0 && length(group1) <= 3)
    stopifnot(length(group2) > 0 && length(group2) <= 3)
    groupBoth <- c(group1, group2)
    if(length(group1) == 1) {
        idxMatrix1 <- as.matrix(c(group1, group2))
    }
    if(length(group1) == 2) {
        idxMatrix1 <- rbind(group1,
                           cbind(group1[1], group2),
                           cbind(group1[2], group2))
        if(length(group2) == 3)
            idxMatrix1 <- rbind(idxMatrix1,
                                group2[c(1,2)],
                                group2[c(1,3)],
                                group2[c(2,3)])
    }
    if(length(group1) == 3) {
        idxMatrix1 <- rbind(group1,
                            c(group1[1:2], group2[1]),
                            c(group1[1:2], group2[2]),
                            c(group1[1:2], group2[3]),
                            c(group1[c(1,3)], group2[1]),
                            c(group1[c(1,3)], group2[2]),
                            c(group1[c(1,3)], group2[3]),
                            c(group1[c(2,3)], group2[1]),
                            c(group1[c(2,3)], group2[2]),
                            c(group1[c(2,3)], group2[3]))
    }
    idxMatrix2 <- do.call(rbind, lapply(1:nrow(idxMatrix1), function(ii) {
        setdiff(groupBoth, idxMatrix1[ii,])
    }))
    return(list(idxMatrix1 = idxMatrix1, idxMatrix2 = idxMatrix2))
}
