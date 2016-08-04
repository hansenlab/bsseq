getNullDistribution_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2,
                                              estimate.var, local.correct,
                                              cutoff, stat, maxGap, mc.cores = 1) {
    stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
    message(sprintf("[getNullDistribution_BSmooth.tstat] performing %d permutations\n", nrow(idxMatrix1)))
    nullDist <- mclapply(1:nrow(idxMatrix1), function(ii) {
        ptime1 <- proc.time()
        BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var,
                                  group1 = idxMatrix1[ii,],
                                  group2 = idxMatrix2[ii,],
                                  local.correct = local.correct, maxGap = 10^8,
                                  verbose = FALSE)
        dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = maxGap)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        message(sprintf("[getNullDistribution_BSmooth.tstat] completing permutation %d in %.1f sec\n", ii, stime))
        dmrs0
    }, mc.cores = min(nrow(idxMatrix1), mc.cores), mc.preschedule = FALSE)
    nullDist
}

getNullDistribution_BSmooth.fstat <- function(BSseq,
                                              idxMatrix,
                                              design, contrasts,
                                              coef = NULL, cutoff,
                                              maxGap.sd, maxGap.dmr,
                                              mc.cores = 1) {

    message(sprintf("[getNullDistribution_BSmooth.fstat] performing %d permutations\n",
                    nrow(idxMatrix)))
    # NOTE: Using mc.preschedule = TRUE
    # TODO: Need some protection for when a core(s) error, otherwise nullDist
    #       contains a mixture of data frame, NULL, and try-error objects
    #       (which will subesequently break when passed to getFWER.fstat())
    nullDist <- mclapply(seq_len(nrow(idxMatrix)), function(ii) {
        ptime1 <- proc.time()
        # NOTE: More efficient to permute design matrix using idxMatrix[ii, ] than
        #       permute the raw data with BSseq[, idxMatrix[ii, ]]
        bstat <- BSmooth.fstat(BSseq, design = design[idxMatrix[ii, ], ],
                               contrasts = contrasts)
        bstat <- smoothSds(bstat, maxGap = maxGap.sd)
        bstat <- computeStat(bstat, coef = coef)
        dmrs0 <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        message(sprintf("[getNullDistribution_BSmooth.fstat] completing permutation %d in %.1f sec\n", ii, stime))
        dmrs0
    }, mc.cores = min(nrow(idxMatrix), mc.cores), mc.preschedule = TRUE)
    nullDist
}

permuteAll <- function(nperm, design) {
    message(sprintf("[permuteAll] performing %d unrestricted permutations of the design matrix\n", nperm))

    CTRL <- how(nperm = nperm)
    # NOTE: shuffleSet() returns a nperm * nrow(design) matrix of permutations
    idxMatrix <- shuffleSet(n = nrow(design), control = CTRL)
}

getNullBlocks_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                          mc.cores = 1) {
    getNullDistribution_BSmooth.tstat(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = FALSE,
                       estimate.var = estimate.var,
                       cutoff = c(-2,2), stat = "tstat", maxGap = 10000,
                       mc.cores = mc.cores)
}

getNullDmrs_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                        mc.cores = 1) {
    getNullDistribution_BSmooth.tstat(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = TRUE,
                       estimate.var = estimate.var,
                       cutoff = c(-4.6,4.6), stat = "tstat.corrected", maxGap = 300,
                       mc.cores = mc.cores)
}




subsetDmrs <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    out <- xx[ xx[,"n"] >= 3 & abs(xx[, "meanDiff"]) > 0.1 &
                  xx[, "invdensity"] <= 300, ]
    if(nrow(out) == 0)
        return(NULL)
    out
}

subsetBlocks <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    out <- subset(xx, width >= 10000)
    if(nrow(out) == 0)
        return(NULL)
    out
}

getFWER <- function(null, type = "blocks") {
    reference <- null[[1]]
    null <- null[-1]
    null <- null[!sapply(null, is.null)]
    better <- sapply(1:nrow(reference), function(ii) {
        # meanDiff <- abs(reference$meanDiff[ii])
        areaStat <- abs(reference$areaStat[ii])
        width <- reference$width[ii]
        n <- reference$n[ii]
        if (type == "blocks") {
            out <- sapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                        # nulldist$width >= width)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$width >= width)
            })
        }
        if (type == "dmrs") {
            out <- sapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                #     nulldist$n >= n)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$n >= n)
            })
        }
        sum(out)
    })
    better
}

# NOTE: Identical to getFWER() except uses areaStat rather than meanDiff
#       to compare regions.
getFWER.fstat <- function(null, type = "blocks") {
    reference <- null[[1]]
    null <- null[-1]
    null <- null[!sapply(null, is.null)]
    # TODO: Will break if null == list(), which can occur in practice (although
    #       rarely).
    better <- sapply(seq_len(nrow(reference)), function(ii) {
        # meanDiff <- abs(reference$meanDiff[ii])
        areaStat <- abs(reference$areaStat[ii])
        width <- reference$width[ii]
        n <- reference$n[ii]
        if (type == "blocks") {
            out <- vapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                # nulldist$width >= width)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$width >= width)
            }, logical(1L))
        }
        if (type == "dmrs") {
            out <- vapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                #     nulldist$n >= n)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$n >= n)
            }, logical(1L))
        }
        sum(out)
    })
    better
}

# TODO: Simplify makeIdxMatrix() by using permute package
makeIdxMatrix <- function(group1, group2, testIsSymmetric = TRUE, includeUnbalanced = TRUE) {
    groupBoth <- c(group1, group2)
    idxMatrix1 <- NULL
    subsetByMatrix <- function(vec, mat) {
        apply(mat, 2, function(xx) vec[xx])
    }
    combineMat <- function(mat1, mat2) {
        tmp <- lapply(1:nrow(mat1), function(ii) {
            t(apply(mat2, 1, function(xx) { c(mat1[ii,], xx) }))
        })
        do.call(rbind, tmp)
    }
    if(length(group1) == 1 && length(group1) == 1) {
        if(testIsSymmetric)
            idxMatrix1 <- as.matrix(group1)
        else
            idxMatrix1 <- as.matrix(c(group1, group2))
    }
    if(length(group1) == 2 && length(group2) == 2) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(matrix(group1[1], ncol = 1),
                                           matrix(group2, ncol = 1)))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(matrix(group1, ncol = 1), matrix(group2, ncol = 1)))
        }
    }
    if(length(group1) == 3 && length(group1) == 3) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           as.matrix(group2, ncol = 1)))

        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           as.matrix(group2, ncol = 1)),
                                combineMat(as.matrix(group1, ncol = 1),
                                           subsetByMatrix(group2, combinations(3,2))))
        }
    }
    if(length(group1) == 4 && length(group1) == 4) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           subsetByMatrix(group2, combinations(4,2))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(4,2)),
                                           subsetByMatrix(group2, combinations(4,2))))
        }
        if(includeUnbalanced) {
            newMatrix <- combineMat(subsetByMatrix(group1, combinations(4,3)),
                                    as.matrix(group2, ncol = 1))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix)
        }
        if(includeUnbalanced && !testIsSymmetric) {
            newMatrix <- combineMat(as.matrix(group1, ncol = 1),
                                    subsetByMatrix(group2, combinations(4,3)))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix)
        }
    }
    if(length(group1) == 5 && length(group1) == 5) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(5, 3)),
                                           subsetByMatrix(group2, combinations(5, 2))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(5, 3)),
                                           subsetByMatrix(group2, combinations(5, 2))),
                                combineMat(subsetByMatrix(group1, combinations(5, 2)),
                                           subsetByMatrix(group2, combinations(5, 3))))
        }
        if(includeUnbalanced) {
            idxMatrix1 <- rbind(idxMatrix1,
                                combineMat(subsetByMatrix(group1, combinations(5,4)),
                                           as.matrix(group2, ncol = 1)))
        }
        if(includeUnbalanced && !testIsSymmetric) {
            idxMatrix1 <- rbind(idxMatrix1,
                                combineMat(as.matrix(group1, ncol = 1),
                                           subsetByMatrix(group2, combinations(5,4))))
        }
    }
    if(length(group1) == 6 && length(group1) == 6) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(5,3)),
                                           subsetByMatrix(group2, combinations(6,3))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(6,3)),
                                           subsetByMatrix(group2, combinations(6,3))))
        }
        if(includeUnbalanced) {
            newMatrix1 <- combineMat(subsetByMatrix(group1, combinations(6,4)),
                                    subsetByMatrix(group2, combinations(6,2)))
            newMatrix2 <- combineMat(subsetByMatrix(group1, combinations(6,5)),
                                    as.matrix(group2, ncol = 1))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix1, newMatrix2)
        }
        if(includeUnbalanced && !testIsSymmetric) {
            newMatrix1 <- combineMat(subsetByMatrix(group1, combinations(6,2)),
                                     subsetByMatrix(group2, combinations(6,4)))
            newMatrix2 <- combineMat(as.matrix(group1, ncol = 1),
                                     subsetByMatrix(group2, combinations(6,5)))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix1, newMatrix2)
        }
    }
    if(is.null(idxMatrix1))
        stop("unable to handle this combination of 'group1', 'group2' and 'testIsSymmetric'")
    rownames(idxMatrix1) <- NULL
    idxMatrix2 <- do.call(rbind, lapply(1:nrow(idxMatrix1), function(ii) {
        setdiff(groupBoth, idxMatrix1[ii,])
    }))
    return(list(idxMatrix1 = idxMatrix1, idxMatrix2 = idxMatrix2))
}
