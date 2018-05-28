context("Test fixes")

test_that("test_fixes", {
    gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 1))
    gr2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 2:3, width = 1))
    ov <- findOverlaps(gr1, gr2)
    expect_equal(queryHits(ov), 2)
    expect_equal(subjectHits(ov), 1)
    expect_equal(subsetByOverlaps(gr1, gr2), gr1[2])
})

test_that("BSseq", {
    M <- matrix(1:9, 3, 3)
    colnames(M) <- c("A1", "A2", "A3")
    BStest <- BSseq(pos = 1:3, chr = rep("chr1", 3), M = M, Cov = M + 2)
    expect_equal(seqnames(BStest), Rle(factor(rep("chr1", 3))))
    expect_equal(start(BStest), c(1, 2, 3))
    expect_equal(end(BStest), c(1, 2, 3))
    expect_equal(width(BStest), rep(1, 3))
    expect_equal(strand(BStest),
                 Rle(factor(rep("*", 3), levels = c("+", "-", "*"))))
    expect_equal(dim(BStest), c(3, 3))
    expect_equal(nrow(BStest), 3)
    expect_equal(ncol(BStest), 3)
    expect_equal(unname(as.array(getCoverage(BStest, type = "M"))), unname(M))
    expect_equal(unname(as.array(getCoverage(BStest, type = "Cov"))),
                 unname(M + 2))
    expect_equal(sampleNames(BStest), colnames(M))

    BStest2 <- BSseq(pos = 3:1, chr = rep("chr1", 3), M = M[3:1, ],
                     Cov = M[3:1, ] + 2)
    expect_equal(BStest, BStest2)
    M2 <- rbind(M[3:1, ], c(1, 1, 1))
    Cov2 <- rbind(M[3:1, ] + 2, c(1, 1, 1))
    M2[3, ] <- M2[3, ] - 1
    Cov2[3, ] <- Cov2[3, ] - 1
    BStest3 <- suppressWarnings(
        BSseq(pos = c(3:1,1), chr = rep("chr1", 4), M = M2, Cov = Cov2))
    expect_equal(BStest, BStest3)
})

test_that("BSseq with HDF5Backend", {
    M <- matrix(1:9, 3, 3)
    colnames(M) <- c("A1", "A2", "A3")
    Cov <- M + 2L
    BStest <- BSseq(pos = 1:3, chr = rep("chr1", 3), M = M, Cov = Cov)
    hdf5_M <- realize(M, BACKEND = "HDF5Array")
    hdf5_Cov <- realize(Cov, BACKEND = "HDF5Array")
    hdf5_BStest <- BSseq(pos = 1:3, chr = rep("chr1", 3), M = hdf5_M,
                         Cov = hdf5_Cov)
    expect_identical(as.array(getCoverage(BStest)),
                     as.array(getCoverage(hdf5_BStest)))
    expect_identical(as.array(getCoverage(BStest, type = "M")),
                     as.array(getCoverage(hdf5_BStest, type = "M")))
})

test_overlaps <- function() {
    M <- matrix(1:9, 3, 3)
    colnames(M) <- c("A1", "A2", "A3")
    Cov <- M + 2
    gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:3, width = 1))
    gr2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 3:5, width = 1))
    BS1 <- BSseq(gr = gr1, M = M, Cov = Cov)
    BS2 <- BSseq(gr = gr2, M = M, Cov = Cov)
    BSres <- BSseq(gr = gr1[3],
                   M = M[3, ,drop = FALSE],
                   Cov = Cov[3, , drop = FALSE])
    expect_equal(length(BS1), 3L)
    expect_equal(seqlevels(BS1), seqlevels(gr1))
    expect_equal(seqlengths(BS1), seqlengths(gr1))
    expect_equal(seqnames(BS1), seqnames(gr1))
    expect_equal(granges(BS1), gr1)
    expect_equal(start(BS1), start(gr1))
    expect_equal(end(BS1), end(gr1))
    expect_equal(width(BS1), width(gr1))
    expect_equal(strand(BS1), strand(gr1))
    expect_equal(findOverlaps(BS1, BS2), findOverlaps(gr1, gr2))
    expect_equal(findOverlaps(gr1, BS2), findOverlaps(gr1, gr2))
    expect_equal(findOverlaps(BS1, gr2), findOverlaps(gr1, gr2))
    expect_equal(subsetByOverlaps(BS1, gr2), BSres)
    expect_equal(subsetByOverlaps(BS1, BS2), BSres)
    expect_equal(subsetByOverlaps(gr1, BS2), gr1[3])
}
