context("read.modbam2bed")

# TODO: Re-factor read.modbam2bed() and update tests accordingly
test_that("read.modbam2bed() works for BED files", {
  infile <- system.file("extdata", "modbam2bed/Ctr1.bed.gz",package = "bsseq")
  bsseq <- read.modbam2bed(files = infile,
                        colData = NULL,
                        rmZeroCov = FALSE,
                        strandCollapse = TRUE)
  expect_is(bsseq, "BSseq")
})
