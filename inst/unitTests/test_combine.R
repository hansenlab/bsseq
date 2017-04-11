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

test_.subassignRowsDelayedMatrix <- function() {
    nrow <- 1000L
    ncol <- 10L
    x <- writeHDF5Array(matrix(seq_len(nrow * ncol), ncol = ncol))
    x_i <- seq(1L, 2L * nrow, 2L)
    y <- writeHDF5Array(matrix(seq(-1L, -nrow * ncol, -1L), ncol = 10))
    y_i <- seq(2L, nrow, 2L)

    z1 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
                                              i = x_i,
                                              nrow = 2L * nrow,
                                              fill = NA_integer_,
                                              BACKEND = NULL)
    z2 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
                                              i = x_i,
                                              nrow = 2L * nrow,
                                              fill = NA_integer_,
                                              BACKEND = "HDF5Array",
                                              by_row = FALSE)
    z3 <- bsseq:::.subassignRowsDelayedMatrix(x = x,
                                              i = x_i,
                                              nrow = 2L * nrow,
                                              fill = NA_integer_,
                                              BACKEND = "HDF5Array",
                                              by_row = TRUE)
    checkIdentical(as.array(z1), as.array(z2))
    checkIdentical(as.array(z1), as.array(z3))
}

test_.combineListOfDelayedMatrixObjects <- function() {
    nrow <- 10
    ncol <- 4
    x <- matrix(seq_len(nrow),
                ncol = ncol / 2,
                dimnames = list(NULL, letters[1:2]))
    y <- matrix(100L + seq_len(nrow),
                ncol = ncol / 2,
                dimnames = list(NULL, letters[3:4]))
    x_i <- seq(1, nrow, ncol / 2)
    y_i <- seq(2, nrow, ncol / 2)
    fill <- NA_integer_

    # The expected output
    z <- matrix(fill,
                nrow = nrow,
                ncol = ncol,
                dimnames = list(NULL, letters[seq_len(ncol)]))
    # NOTE: as.array(x) is a no-op if x is a matrix and realises a
    #       DelayedMtrix in memory
    z[x_i, seq(1, ncol(x))] <- x
    z[y_i, seq(ncol(x) + 1, ncol(x) + ncol(y))] <- y

    # # Test with in-memory DelayedMatrix objects
    X <- bsseq:::.DelayedMatrix(x)
    Y <- bsseq:::.DelayedMatrix(y)

    Z <- bsseq:::.combineListOfDelayedMatrixObjects(
        X = list(X, Y),
        I = list(x_i, y_i),
        nrow = nrow,
        ncol = ncol,
        dimnames = list(NULL, c(colnames(X), colnames(Y))),
        fill = fill,
        BACKEND = NULL)
    checkIdentical(z, as.array(Z))
    checkTrue(!bsseq:::.isHDF5ArrayBacked(Z))


    # Test with HDF5Array-backed DelayedMatrix objects
    hdf5_X <- realize(X, BACKEND = "HDF5Array")
    hdf5_Y <- realize(Y, BACKEND = "HDF5Array")

    hdf5_Z <- bsseq:::.combineListOfDelayedMatrixObjects(
        X = list(hdf5_X, hdf5_Y),
        I = list(x_i, y_i),
        nrow = nrow,
        ncol = ncol,
        dimnames = list(NULL, c(colnames(hdf5_X), colnames(hdf5_Y))),
        fill = fill,
        BACKEND = "HDF5Array")
    checkIdentical(z, as.array(hdf5_Z))
    checkTrue(bsseq:::.isHDF5ArrayBacked(hdf5_Z))
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
