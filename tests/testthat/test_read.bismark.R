context("read.bismark")

# TODO: Re-factor read.bismark() and update tests accordingly
test_that("read.bismark() works for 'coverage' file", {
    infile <- system.file("extdata", "test_data.fastq_bismark.bismark.cov.gz",
                          package = "bsseq")
    bsseq <- read.bismark(files = infile,
                          colData = DataFrame(row.names = "test_data"),
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          verbose = FALSE)
    expect_is(bsseq, "BSseq")
})

test_that("read.bismark() works for 'genome wide cytosine report' file", {
    # Test that read.bismark() works on good input
    infile <- system.file("extdata", "test_data.cytosineReport.gz",
                          package = "bsseq")
    bsseq <- read.bismark(files = infile,
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          verbose = FALSE)
    expect_is(bsseq, "BSseq")
    # Check that strandCollapse = FALSE works
    expect_equal(nrow(bsseq), 100L)
    expect_true(all(strand(bsseq) %in% c("+", "-")))

    # Check that strandCollapse = TRUE works
    bsseq <- read.bismark(files = infile,
                          rmZeroCov = FALSE,
                          strandCollapse = TRUE,
                          verbose = FALSE)
    expect_equal(nrow(bsseq), 50L)
    expect_true(all(strand(bsseq) == "*"))
})
