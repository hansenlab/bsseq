context("read.bedMethyl")

test_that("read.bedMethyl() works for bedMethyl file", {
  infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
    package = "bsseq"),
  system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
    package = "bsseq"))
  bsseq <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
     rmZeroCov = TRUE,
     strandCollapse = TRUE,
     verbose = FALSE)
    expect_is(bsseq, "BSseq")
    expect_equal(nrow(bsseq), 1111L)
    expect_true(all(strand(bsseq) == "*"))
    #test of get CpGs
    bsseq_nano<-bsseq[,1]
    expect_equal(length(getCpGs(bsseq_nano, type="allCpG")), 828L)
    CpGMatrix<-getCpGMatrix(bsseq, allCpG = TRUE)
    MaxLikelihoodMatrix<-getMaxLikelihoodMatrix(bsseq, allCpG = TRUE)
    expect_equal(length(bsseq[which(rowAlls(CpGMatrix == 0) 
                                      & rowMins(MaxLikelihoodMatrix) > 0.99)]), 827L)
})
