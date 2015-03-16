test_makeIdxMatrix <- function() {
    grp1 <- paste0("A", 1:5)
    grp2 <- paste0("B", 1:5)

    res <- bsseq:::makeIdxMatrix(grp1[1], grp2[1], testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0L)
    checkEquals(nrow(res), choose(2,1)) # Fail

    res <- bsseq:::makeIdxMatrix(grp1[1], grp2[1], testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0L)
    checkEquals(nrow(res), choose(2,1) / 2)

    res <- bsseq:::makeIdxMatrix(grp1[1:2], grp2[1:2], testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(4,2))

    res <- bsseq:::makeIdxMatrix(grp1[1:2], grp2[1:2], testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(4,2) / 2)
    
    res <- bsseq:::makeIdxMatrix(grp1[1:3], grp2[1:3], testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(6,3))

    res <- bsseq:::makeIdxMatrix(grp1[1:3], grp2[1:3], testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(6,3) / 2)

    res <- bsseq:::makeIdxMatrix(grp1[1:4], grp2[1:4],
                                 testIsSymmetric = FALSE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(8,4))

    res <- bsseq:::makeIdxMatrix(grp1[1:4], grp2[1:4],
                                 testIsSymmetric = TRUE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(8,4) / 2) # Fail

    res <- bsseq:::makeIdxMatrix(grp1[1:5], grp2[1:5],
                                 testIsSymmetric = FALSE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(10,5))

    res <- bsseq:::makeIdxMatrix(grp1[1:5], grp2[1:5],
                                 testIsSymmetric = TRUE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(nrow(res), choose(10,5) / 2)
}
