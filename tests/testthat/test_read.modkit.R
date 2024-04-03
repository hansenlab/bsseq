context("read.modkit")

# TODO: Re-factor read.modkit() and update tests accordingly
test_that("read.modkit() works for BED files without 5hmc", {
    infile <- system.file("extdata", "modkit/chr21.chr22.HG002.top1000.other_mod.bed.gz",
                          package = "bsseq")
    bsseq <- read.modkit(files = infile,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE)

    lapply(bsseq, function(x) {expect_is(x, "BSseq")})
})

test_that("read.modkit() works for BED files with 5hmc", {
    infile <- system.file("extdata", "modkit/Hypo1.first40Bed.txt",
                          package = "bsseq")
    bsseq <- read.modkit(files = infile,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE)

    lapply(bsseq, function(x) {expect_is(x, "BSseq")})
})
