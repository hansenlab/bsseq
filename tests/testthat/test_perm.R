context("Permutations")

test_that("makeIdxMatrix() n = 1", {
    grp1 <- paste0("A", 1)
    grp2 <- paste0("B", 1)

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    expect_equal(anyDuplicated(res), 0L)
    expect_equal(res[1, ], grp1)
    expect_equal(res[2, ], grp2)
    expect_equal(nrow(res), choose(2, 1)) # Fail

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0L)
    expect_equal(res[1, ], grp1)
    expect_equal(nrow(res), choose(2, 1) / 2)
})

test_that("makeIdxMatrix() n = 2", {
    grp1 <- paste0("A", 1:2)
    grp2 <- paste0("B", 1:2)

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(res[2,], grp2)
    expect_equal(nrow(res), choose(4, 2))

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(nrow(res), choose(4, 2) / 2)
})

test_that("makeIdxMatrix() n = 3", {
    grp1 <- paste0("A", 1:3)
    grp2 <- paste0("B", 1:3)

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = FALSE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(res[2,], grp2)
    expect_equal(nrow(res), choose(6, 3))

    res <- bsseq:::makeIdxMatrix(grp1, grp2, testIsSymmetric = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(nrow(res), choose(6, 3) / 2)
})

test_that("makeIdxMatrix() n = 4", {
    grp1 <- paste0("A", 1:4)
    grp2 <- paste0("B", 1:4)

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = FALSE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(res[2,], grp2)
    expect_equal(nrow(res), choose(8,4))

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = TRUE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(nrow(res), choose(8,4) / 2) # Fail
})

test_that("makeIdxMatrix() n = 5", {

    grp1 <- paste0("A", 1:5)
    grp2 <- paste0("B", 1:5)

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = FALSE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(res[2,], grp2)
    expect_equal(nrow(res), choose(10,5))

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = TRUE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(nrow(res), choose(10,5) / 2)
})

test_that("makeIdxMatrix() n = 6", {
    grp1 <- paste0("A", 1:6)
    grp2 <- paste0("B", 1:6)

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = FALSE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(res[2,], grp2)
    expect_equal(nrow(res), choose(12,6))

    res <- bsseq:::makeIdxMatrix(
        grp1,
        grp2,
        testIsSymmetric = TRUE,
        includeUnbalanced = TRUE)[[1]]
    expect_equal(anyDuplicated(res), 0)
    expect_equal(res[1,], grp1)
    expect_equal(nrow(res), choose(12,6) / 2)
})
