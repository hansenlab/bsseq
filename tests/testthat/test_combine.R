context("combine()")

# TODO: .subassignRowsDelayedMatrix() is no longer needed in bsseq; test can be
#       removed.
# test_that(".subassignRowsDelayedMatrix()", {
#     nrow <- 1000L
#     ncol <- 10L
#     x <- realize(matrix(seq_len(nrow * ncol), ncol = ncol), "HDF5Array")
#     x_i <- seq(1L, 2L * nrow, 2L)
#     y <- realize(matrix(seq(-1L, -nrow * ncol, -1L), ncol = 10), "HDF5Array")
#     y_i <- seq(2L, nrow, 2L)
#
#     z1 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
#                                               i = x_i,
#                                               nrow = 2L * nrow,
#                                               fill = NA_integer_,
#                                               BACKEND = NULL)
#     z2 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
#                                               i = x_i,
#                                               nrow = 2L * nrow,
#                                               fill = NA_integer_,
#                                               BACKEND = "HDF5Array",
#                                               by_row = FALSE)
#     z3 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
#                                               i = x_i,
#                                               nrow = 2L * nrow,
#                                               fill = NA_integer_,
#                                               BACKEND = "HDF5Array",
#                                               by_row = TRUE)
#     expect_identical(as.array(z1), as.array(z2))
#     expect_identical(as.array(z1), as.array(z3))
# })

# TODO: .combineListOfDelayedMatrixObjects() is no longer needed in bsseq; test
#        can be removed.
# test_that(".combineListOfDelayedMatrixObjects()", {
#     nrow <- 10
#     ncol <- 4
#     x <- matrix(seq_len(nrow),
#                 ncol = ncol / 2,
#                 dimnames = list(NULL, letters[1:2]))
#     y <- matrix(100L + seq_len(nrow),
#                 ncol = ncol / 2,
#                 dimnames = list(NULL, letters[3:4]))
#     x_i <- seq(1, nrow, ncol / 2)
#     y_i <- seq(2, nrow, ncol / 2)
#     fill <- NA_integer_
#
#     # The expected output
#     z <- matrix(fill,
#                 nrow = nrow,
#                 ncol = ncol,
#                 dimnames = list(NULL, letters[seq_len(ncol)]))
#     # NOTE: as.array(x) is a no-op if x is a matrix and realises a
#     #       DelayedMtrix in memory
#     z[x_i, seq(1, ncol(x))] <- x
#     z[y_i, seq(ncol(x) + 1, ncol(x) + ncol(y))] <- y
#
#     # # Test with in-memory DelayedMatrix objects
#     X <- bsseq:::.DelayedMatrix(x)
#     Y <- bsseq:::.DelayedMatrix(y)
#
#     Z <- bsseq:::.combineListOfDelayedMatrixObjects(
#         X = list(X, Y),
#         I = list(x_i, y_i),
#         nrow = nrow,
#         ncol = ncol,
#         dimnames = list(NULL, c(colnames(X), colnames(Y))),
#         fill = fill,
#         BACKEND = NULL)
#     expect_identical(z, as.array(Z))
#     expect_true(!bsseq:::.isHDF5ArrayBacked(Z))
#
#
#     # Test with HDF5Array-backed DelayedMatrix objects
#     hdf5_X <- realize(X, BACKEND = "HDF5Array")
#     hdf5_Y <- realize(Y, BACKEND = "HDF5Array")
#
#     hdf5_Z <- bsseq:::.combineListOfDelayedMatrixObjects(
#         X = list(hdf5_X, hdf5_Y),
#         I = list(x_i, y_i),
#         nrow = nrow,
#         ncol = ncol,
#         dimnames = list(NULL, c(colnames(hdf5_X), colnames(hdf5_Y))),
#         fill = fill,
#         BACKEND = "HDF5Array")
#     expect_identical(z, as.array(hdf5_Z))
#     expect_true(bsseq:::.isHDF5ArrayBacked(hdf5_Z))
# })

checkBSseqAssaysIdentical <- function(x, y) {
    stopifnot(is(x, "BSseq") && is(y, "BSseq"))
    assay_names <- c("M", "Cov", "coef", "se.coef")
    check_identical <- vapply(assay_names, function(an) {
        if (!is.null(getBSseq(x, an))) {
            identical(as.array(getBSseq(x, an)), as.array(getBSseq(y, an)))
        } else {
            identical(getBSseq(x, an), getBSseq(y, an))
        }
    }, logical(1L))
    checkTrue(all(check_identical))
}

checkBSseqIdentical <- function(x, y) {
    checkTrue(identical(rowRanges(x), rowRanges(y)) &&
                  identical(getBSseq(x, "trans"), getBSseq(y, "trans")) &&
                  identical(getBSseq(x, "parameters"),
                            getBSseq(y, "parameters")) &&
                  checkBSseqAssaysIdentical(x, y))
}

test_that("combine()", {
    bsseq_fit <- BSmooth(bsseq_test)
    BSSEQ_TEST <- realize(bsseq_test, "HDF5Array")
    BSSEQ_FIT <- realize(bsseq_fit, "HDF5Array")

    ai <- 1:100
    bi <- 51:150
    a <- bsseq_test[ai, 1]
    b <- bsseq_test[bi, 2]
    ab <- combine(a, b)
    A <- BSSEQ_TEST[ai, 1]
    B <- BSSEQ_TEST[bi, 2]
    AB <- combine(A, B)
    expect_equivalent_SE(ab, AB)
    aB <- combine(a, B)
    expect_equivalent_SE(ab, aB)
    Ab <- combine(A, b)
    expect_equivalent_SE(ab, Ab)

    z <- combine(bsseq_fit[, 1], bsseq_fit[, 2])
    expect_equivalent_SE(z, bsseq_fit)
    Z <- combine(BSSEQ_FIT[, 1], BSSEQ_FIT[, 2])
    expect_equivalent_SE(Z, BSSEQ_FIT)
})

test_that("combineList()", {
    bsseq_fit <- BSmooth(bsseq_test)
    BSSEQ_TEST <- realize(x = bsseq_test, "HDF5Array")
    BSSEQ_FIT <- realize(x = bsseq_fit, "HDF5Array")

    ai <- 1:100
    bi <- 51:150
    ci <- 201:300
    a <- bsseq_test[ai, 1]
    b <- bsseq_test[bi, 2]
    c <- bsseq_test[ci, 1]
    colnames(c) <- "r3"
    abc <- combineList(list(a, b, c))
    A <- BSSEQ_TEST[ai, 1]
    B <- BSSEQ_TEST[bi, 2]
    C <- BSSEQ_TEST[ci, 1]
    colnames(C) <- "r3"
    ABC <- combineList(list(A, B, C))
    expect_equivalent_SE(abc, ABC)
    expect_equivalent_SE(combineList(list(a, b, c)),
                         combineList(a, b, c))
    expect_equivalent_SE(combineList(list(a, B, c)),
                         combineList(A, b, C))

    expect_equivalent_SE(combine(a, b),
                         combineList(a, b))
    expect_equivalent_SE(combine(A, B),
                         combineList(A, B))

    z <- combineList(bsseq_fit[, 1], bsseq_fit[, 2])
    expect_equivalent_SE(z, bsseq_fit)
    Z <- combineList(BSSEQ_FIT[, 1], BSSEQ_FIT[, 2])
    expect_equivalent_SE(Z, BSSEQ_FIT)
})
test_that(
    "Test bug fix reported in https://github.com/hansenlab/bsseq/pull/54/", {
        M <- matrix(0:8, 3, 3)
        Cov <- matrix(1:9, 3, 3)

        BS1 <- BSseq(chr = c("chr1", "chr2", "chr1"),
                     pos = c(1, 2, 3),
                     M = M,
                     Cov = Cov,
                     sampleNames = c("A","B", "C"))
        BS2 <- BSseq(chr = c("chr2", "chr1", "chr2"),
                     pos = c(1, 2, 3),
                     M = M,
                     Cov = Cov, sampleNames = c("D", "E", "F"))
        # NOTE: Using combineList()
        combined <- combineList(list(BS1, BS2))
        over <- findOverlaps(BS2, combined, type = "equal")
        methAlone <- getMeth( BS2[queryHits(over), c("D", "E", "F")], type = "raw")
        methCombined <- getMeth(combined[subjectHits(over), c("D", "E", "F")],
                                type = "raw")
        expect_equal(as.matrix(methAlone), as.matrix(methCombined))
        # NOTE: Using combine()
        combined2 <- combine(BS1, BS2)
        over2 <- findOverlaps(BS2, combined2, type = "equal")
        methAlone2 <- getMeth( BS2[queryHits(over2), c("D", "E", "F")],
                               type = "raw")
        methCombined2 <- getMeth(combined2[subjectHits(over2), c("D", "E", "F")],
                                 type = "raw")
        expect_equal(as.matrix(methAlone2), as.matrix(methCombined2))

        # Extended version with 3 BSseq objects
        BS3 <- BSseq(chr = c("chr2", "chr1", "chr2"),
                     pos = c(3, 2, 1),
                     M = M,
                     Cov = Cov,
                     sampleNames = c("G", "H", "I")
        )
        # NOTE: Using combineList()
        combined <- combineList(BS2, BS1, BS3)
        over <- findOverlaps(BS3, combined, type = "equal")
        methAlone <- getMeth(BS3[queryHits(over), c("G", "H", "I")], type = "raw")
        methCombined <- getMeth(combined[subjectHits(over), c("G", "H", "I")],
                                type = "raw")
        expect_equal(as.matrix(methAlone), as.matrix(methCombined))
        # NOTE: Using combine()
        combined2 <- combine(BS2, BS1, BS3)
        over2 <- findOverlaps(BS3, combined2, type = "equal")
        methAlone2 <- getMeth(BS3[queryHits(over2), c("G", "H", "I")], type = "raw")
        methCombined2 <- getMeth(combined2[subjectHits(over2), c("G", "H", "I")],
                                 type = "raw")
        expect_equal(as.matrix(methAlone), as.matrix(methCombined))
    }
)
