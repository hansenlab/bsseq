test_makeIdxMatrix_perm1 <- function() {
    grp1 <- paste0("A", 1)
    grp2 <- paste0("B", 1)

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0L)
    checkEquals(res[1,], grp1)
    checkEquals(res[2,], grp2)
    checkEquals(nrow(res), choose(2,1)) # Fail

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0L)
    checkEquals(res[1,], grp1)
    checkEquals(nrow(res), choose(2,1) / 2)
}

test_makeIdxMatrix_perm2 <- function() {
    grp1 <- paste0("A", 1:2)
    grp2 <- paste0("B", 1:2)
    
    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(res[2,], grp2)
    checkEquals(nrow(res), choose(4,2))
    
    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(nrow(res), choose(4,2) / 2)
}

test_makeIdxMatrix_perm3 <- function() {
    grp1 <- paste0("A", 1:3)
    grp2 <- paste0("B", 1:3)

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(res[2,], grp2)
    checkEquals(nrow(res), choose(6,3))

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(nrow(res), choose(6,3) / 2)
}
    
test_makeIdxMatrix_perm4 <- function() {
    grp1 <- paste0("A", 1:4)
    grp2 <- paste0("B", 1:4)

    res <- bsseq:::makeIdxMatrix(grp1, grp2,
                                 testIsSymmetric = FALSE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(res[2,], grp2)
    checkEquals(nrow(res), choose(8,4))

    res <- bsseq:::makeIdxMatrix(grp1, grp2,
                                 testIsSymmetric = TRUE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(nrow(res), choose(8,4) / 2) # Fail
}
    
test_makeIdxMatrix_perm5 <- function() {
    grp1 <- paste0("A", 1:5)
    grp2 <- paste0("B", 1:5)

    res <- bsseq:::makeIdxMatrix(grp1, grp2,
                                 testIsSymmetric = FALSE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(res[2,], grp2)
    checkEquals(nrow(res), choose(10,5))

    res <- bsseq:::makeIdxMatrix(grp1, grp2,
                                 testIsSymmetric = TRUE, includeUnbalanced = TRUE)[[1]]
    checkEquals(anyDuplicated(res), 0)
    checkEquals(res[1,], grp1)
    checkEquals(nrow(res), choose(10,5) / 2)
}

