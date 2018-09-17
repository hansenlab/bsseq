context(".checkMandCov")

test_that(".checkMandCov() works", {
    # Create data
    nrow <- 10L
    ncol <- 2L
    M <- matrix(seq.int(nrow * ncol), nrow = nrow, ncol = ncol)
    Cov <- M + 1L
    M_with_NA_at_start <- M
    M_with_NA_at_start[1L] <- NA_integer_
    Cov_with_NA_at_start <- M
    Cov_with_NA_at_start[1L] <- NA_integer_
    M_with_NA_at_end <- M
    M_with_NA_at_end[nrow * ncol] <- NA_integer_
    Cov_with_NA_at_end <- M
    Cov_with_NA_at_end[nrow * ncol] <- NA_integer_
    M_with_impossible_value_at_start <- M
    M_with_impossible_value_at_start[1L] <- Cov[1L] + 2L
    M_with_impossible_value_at_end <- M
    M_with_impossible_value_at_end[nrow * ncol] <- Cov[nrow * ncol] + 2L
    Cov_with_impossible_value_at_start <- Cov
    Cov_with_impossible_value_at_start[1L] <- M[1L] - 2L
    Cov_with_impossible_value_at_end <- Cov
    Cov_with_impossible_value_at_end[nrow * ncol] <- M[nrow * ncol] - 2L

    # Returns NULL on valid input
    expect_identical(.checkMandCov(M, Cov), NULL)
    expect_identical(.checkMandCov(DelayedArray(M), DelayedArray(Cov)), NULL)
    # NOTE: No check is made that input is strictly integer.
    expect_identical(.checkMandCov(M + 0.1, Cov), NULL)
    expect_identical(.checkMandCov(M, Cov + 0.1), NULL)
    # NOTE: Recall that Matrix::Matrix() stores all data using doubles
    expect_identical(.checkMandCov(Matrix::Matrix(M), Cov), NULL)
    expect_identical(.checkMandCov(M, Matrix::Matrix(Cov)), NULL)

    # Returns msg if anyNA(M) || anyNA(Cov)
    # NOTE: Returns a character vector rather than message() because
    #       .checkMandCov() is used within a validity method.
    expect_identical(.checkMandCov(M_with_NA_at_start, Cov),
                     "'M' must not contain NAs.")
    expect_identical(.checkMandCov(DelayedArray(M_with_NA_at_start), Cov),
                     "'M' must not contain NAs.")
    expect_identical(.checkMandCov(M_with_NA_at_end, Cov),
                     "'M' must not contain NAs.")
    expect_identical(.checkMandCov(DelayedArray(M_with_NA_at_end), Cov),
                     "'M' must not contain NAs.")
    expect_identical(.checkMandCov(M, Cov_with_NA_at_start),
                     "'Cov' must not contain NAs.")
    expect_identical(.checkMandCov(M, DelayedArray(Cov_with_NA_at_start)),
                     "'Cov' must not contain NAs.")
    expect_identical(.checkMandCov(M, Cov_with_NA_at_end),
                     "'Cov' must not contain NAs.")
    expect_identical(.checkMandCov(M, DelayedArray(Cov_with_NA_at_end)),
                     "'Cov' must not contain NAs.")

    # Returns msg if any(M > Cov)
    msg <- "All values of 'M' must be less than or equal to the corresponding value of 'Cov'."
    expect_identical(.checkMandCov(M_with_impossible_value_at_start, Cov), msg)
    expect_identical(
        .checkMandCov(DelayedArray(M_with_impossible_value_at_start), Cov),
        msg)
    expect_identical(.checkMandCov(M_with_impossible_value_at_end, Cov), msg)
    expect_identical(
        .checkMandCov(DelayedArray(M_with_impossible_value_at_end), Cov),
        msg)
    expect_identical(.checkMandCov(M, Cov_with_impossible_value_at_start), msg)
    expect_identical(
        .checkMandCov(M, DelayedArray(Cov_with_impossible_value_at_start)),
        msg)
    expect_identical(.checkMandCov(M, Cov_with_impossible_value_at_end), msg)
    expect_identical(
        .checkMandCov(M, DelayedArray(Cov_with_impossible_value_at_end)),
        msg)
})
