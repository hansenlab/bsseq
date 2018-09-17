context("BSseq-class")

test_that("BSseq() allows loci of width != 1", {
    bsseq <- BSseq(
        M = matrix(1:10),
        Cov = matrix(1:10),
        gr = GRanges("chr1", IRanges(1:10, width = 10), strand = "*"))
    expect_true(validObject(bsseq))
})

test_that("strandCollapse() does nothing for unstranded data", {
    bsseq <- BSseq(
        M = matrix(1:10),
        Cov = matrix(1:10),
        gr = GRanges("chr1", IRanges(1:10, width = 1), strand = "*"))
    expect_warning(strandCollapse(bsseq))
})

test_that("strandCollapse() works on good input", {
    nrow <- 20
    ncol <- 5
    M <- matrix(sample(0:10, size = nrow * ncol, replace = TRUE), ncol = ncol)
    Cov <- M + sample(0:10, size = nrow * ncol, replace = TRUE)
    gr_pos <- GRanges(
        seqnames = c(rep(1, nrow / 2), rep(2, nrow / 2)),
        ranges = seq(1, 2 * nrow, by = 2),
        strand = "+")
    gr_neg <- invertStrand(shift(gr_pos, 1L))
    bsseq_pos <- BSseq(M = M, Cov = Cov, gr = gr_pos)
    bsseq_neg <- BSseq(M = M, Cov = Cov, gr = gr_neg)
    bsseq <- rbind(bsseq_pos, bsseq_neg)

    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_pos, type = "integer"),
        SE2 = unstrand(bsseq_pos))
    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_neg, shift = FALSE, type = "integer"),
        SE2 = unstrand(bsseq_neg))
    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_neg, type = "integer"),
        SE2 = unstrand(shift(bsseq_neg, -1L)))
    bsseq_strand_collapsed <- strandCollapse(bsseq, type = "integer")
    expect_identical(
        object = rowRanges(bsseq_strand_collapsed),
        expected = unstrand(rowRanges(bsseq_pos)))
    expect_identical(
        object = assay(bsseq_strand_collapsed, "M"),
        expected = assay(bsseq_pos, "M") + assay(bsseq_neg, "M"))
    expect_identical(
        object = assay(bsseq_strand_collapsed, "Cov"),
        expected = assay(bsseq_pos, "Cov") + assay(bsseq_neg, "Cov"))

    bsseq_strand_collapsed2 <- strandCollapse(bsseq[sample(nrow(bsseq))])
    expect_equivalent_SE(bsseq_strand_collapsed, bsseq_strand_collapsed2)

    bsseq_pos <- realize(bsseq_pos, "HDF5Array")
    bsseq_neg <- realize(bsseq_neg, "HDF5Array")
    bsseq <- realize(bsseq, "HDF5Array")
    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_pos, type = "integer"),
        SE2 = unstrand(bsseq_pos))
    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_neg, shift = FALSE),
        SE2 = unstrand(bsseq_neg))
    expect_equivalent_SE(
        SE1 = strandCollapse(bsseq_neg),
        SE2 = unstrand(shift(bsseq_neg, -1L)))
    bsseq_strand_collapsed <- strandCollapse(bsseq, type = "integer")
    expect_identical(
        object = rowRanges(bsseq_strand_collapsed),
        expected = unstrand(rowRanges(bsseq_pos)))
    expect_identical(
        object = as.matrix(assay(bsseq_strand_collapsed, "M")),
        expected = as.matrix(assay(bsseq_pos, "M") + assay(bsseq_neg, "M")))
    expect_identical(
        object = as.matrix(assay(bsseq_strand_collapsed, "Cov")),
        expected = as.matrix(assay(bsseq_pos, "Cov") + assay(bsseq_neg, "Cov")))

    bsseq_strand_collapsed2 <- strandCollapse(
        BSseq = bsseq[sample(nrow(bsseq))],
        type = "integer")
    expect_equivalent_SE(bsseq_strand_collapsed, bsseq_strand_collapsed2)
})

test_that("strandCollapse() errors on bad input", {
    bsseq <- BSseq(
        M = matrix(1:3),
        Cov = matrix(1:3),
        gr = GRanges(
            seqnames = "chr1",
            ranges = IRanges(1:3, width = 1),
            strand = c("+", "-", "*")))
    expect_error(strandCollapse(bsseq))
})

test_that("strandCollapse() will unstrand loci and may re-order them", {
    bsseq <- BSseq(
        M = matrix(1:10),
        Cov = matrix(1:10),
        gr = GRanges("chr1", IRanges(10:1, width = 1), strand = "+"))
    expect_true(all(strand(strandCollapse(bsseq)) == "*"))
    expect_false(expect_equivalent_SE(strandCollapse(bsseq), bsseq))
    expect_equivalent_SE(strandCollapse(bsseq), bsseq[10:1])
})
