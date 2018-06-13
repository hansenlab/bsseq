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
                  all.equal(getBSseq(x, "trans"), getBSseq(y, "trans")) &&
                  identical(getBSseq(x, "parameters"),
                            getBSseq(y, "parameters")) &&
                  checkBSseqAssaysIdentical(x, y))
}

test_combine <- function() {
    data(BS.chr22)
    bsseq <- BS.chr22[1:1000, ]
    bsseq_fit <- BSmooth(bsseq, verbose = FALSE)
    BSSEQ <- saveHDF5SummarizedExperiment(x = bsseq, dir = tempfile())
    BSSEQ_FIT <- saveHDF5SummarizedExperiment(x = bsseq_fit, dir = tempfile())

    ai <- 1:100
    bi <- 51:150
    a <- bsseq[ai, 1]
    b <- bsseq[bi, 2]
    ab <- combine(a, b)
    A <- BSSEQ[ai, 1]
    B <- BSSEQ[bi, 2]
    AB <- combine(A, B)
    checkBSseqIdentical(ab, AB)
    aB <- combine(a, B)
    checkBSseqIdentical(ab, aB)
    Ab <- combine(A, b)
    checkBSseqIdentical(ab, Ab)

    z <- combine(bsseq_fit[, 1], bsseq_fit[, 2])
    checkBSseqIdentical(z, bsseq_fit)
    Z <- combine(BSSEQ_FIT[, 1], BSSEQ_FIT[, 2])
    checkBSseqIdentical(Z, BSSEQ_FIT)
}

test_combineList <- function() {
    data(BS.chr22)
    bsseq <- BS.chr22[1:1000, ]
    bsseq_fit <- BSmooth(bsseq, verbose = FALSE)
    BSSEQ <- saveHDF5SummarizedExperiment(x = bsseq, dir = tempfile())
    BSSEQ_FIT <- saveHDF5SummarizedExperiment(x = bsseq_fit, dir = tempfile())

    ai <- 1:100
    bi <- 51:150
    ci <- 201:300
    a <- bsseq[ai, 1]
    b <- bsseq[bi, 2]
    c <- bsseq[ci, 1]
    colnames(c) <- "r3"
    abc <- combineList(list(a, b, c))
    A <- BSSEQ[ai, 1]
    B <- BSSEQ[bi, 2]
    C <- BSSEQ[ci, 1]
    colnames(C) <- "r3"
    ABC <- combineList(list(A, B, C))
    checkBSseqIdentical(abc, ABC)
    checkBSseqIdentical(combineList(list(a, b, c)),
                        combineList(a, b, c))
    checkBSseqIdentical(combineList(list(a, B, c)),
                        combineList(A, b, C))

    checkBSseqIdentical(combine(a, b),
                        combineList(a, b))
    checkBSseqIdentical(combine(A, B),
                        combineList(A, B))

    z <- combineList(bsseq_fit[, 1], bsseq_fit[, 2])
    checkBSseqIdentical(z, bsseq_fit)
    Z <- combineList(BSSEQ_FIT[, 1], BSSEQ_FIT[, 2])
    checkBSseqIdentical(Z, BSSEQ_FIT)
}

# Test bug fix reported in pull request 54
# (https://github.com/kasperdanielhansen/bsseq/pull/54/files)
test_PR54 <- function() {
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
    checkEquals(as.matrix(methAlone), as.matrix(methCombined))
    # NOTE: Using combine()
    combined2 <- combine(BS1, BS2)
    over2 <- findOverlaps(BS2, combined2, type = "equal")
    methAlone2 <- getMeth( BS2[queryHits(over2), c("D", "E", "F")],
                           type = "raw")
    methCombined2 <- getMeth(combined2[subjectHits(over2), c("D", "E", "F")],
                            type = "raw")
    checkEquals(as.matrix(methAlone2), as.matrix(methCombined2))

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
    checkEquals(as.matrix(methAlone), as.matrix(methCombined))
    # NOTE: Using combine()
    combined2 <- combine(combine(BS2, BS1), BS3)
    over2 <- findOverlaps(BS3, combined2, type = "equal")
    methAlone2 <- getMeth(BS3[queryHits(over2), c("G", "H", "I")], type = "raw")
    methCombined2 <- getMeth(combined2[subjectHits(over2), c("G", "H", "I")],
                            type = "raw")
    checkEquals(as.matrix(methAlone), as.matrix(methCombined))
}
