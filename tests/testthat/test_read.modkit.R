context("read.modkit")

# TODO: Re-factor read.modkit() and update tests accordingly
test_that("read.modkit() works for BED files without 5hmc", {
    infile <- system.file("extdata", "modkit/chr21.chr22.HG002.top1000.bed.gz",
                          package = "bsseq")
    bsseq <- read.modkit(files = infile,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE)

    expect_is(bsseq, "BSseq")
})

test_that("read.modkit() works for BED files with 5hmc", {
    infile <- system.file("extdata", "modkit/Hypo1.first50Bed.txt",
                          package = "bsseq")
    bsseq <- read.modkit(files = infile,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE)

    expect_is(bsseq, "BSseq")
})
