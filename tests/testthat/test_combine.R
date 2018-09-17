context("combine()")

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
