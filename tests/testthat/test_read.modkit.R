context("read.modkit")

# TODO: Re-factor read.modkit() and update tests accordingly
test_that("read.modkit() works for BED files", {
    infile <- system.file("extdata", "modkit/chr22.HG002.top1000.other_mod.bed.gz",
                          package = "bsseq")
    bsseq <- read.modkit(files = infile,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE)
    
    lapply(bsseq, function(x) {expect_is(x, "BSseq")})
})
