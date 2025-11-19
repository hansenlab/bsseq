context("read.bedMethyl")

test_that("read.bedMethyl() works for bedMethyl files", {
#load test data
    infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
    package = "bsseq"),
               system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
    package = "bsseq"))

#read in files as bsseq object (Default: rmZeroCov=T, strandCollapse=T)
  bs <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
     rmZeroCov = TRUE,
     strandCollapse = TRUE,
     verbose = FALSE,
     output = "BSseq"
     )

#test class
  expect_is(bs, "BSseq")
#test dimensions
  expect_equal(dim(bs), c(1111L, 2L))
#test strand collapse
  expect_true(all(strand(bs) == "*"))
#test rm zeroCov
  cov <- getCoverage(bs, type="Cov")
  expect_false(any(rowAlls(cov == 0)))
#assay sums
  expect_equal(sum(getBSseq(bs, type="M")), 26632L)
  expect_equal(sum(getBSseq(bs, type="Cov")), 35708L)
getMeth(bs, type = "raw")

#read in files as bsseq object (rmZeroCov=F, strandCollapse=T)
bs <- read.bedMethyl(files = infiles,
                     colData = DataFrame(
                         row.names = c("test_nanopore","test_pacbio")),
                     rmZeroCov = FALSE,
                     strandCollapse = TRUE,
                     verbose = FALSE,
                     output = "BSseq"
)

#test class
expect_is(bs, "BSseq")
#test dimensions
expect_equal(dim(bs), c(1111L, 2L))
#test strand collapse
expect_true(all(strand(bs) == "*"))
#test rm zeroCov
cov <- getCoverage(bs, type="Cov")
expect_false(any(rowAlls(cov == 0)))
#assay sums
expect_equal(sum(getBSseq(bs, type="M")), 26632L)
expect_equal(sum(getBSseq(bs, type="Cov")), 35708L)

#read in files as bsseq object (rmZeroCov=T, strandCollapse=F)
bs <- read.bedMethyl(files = infiles,
                       colData = DataFrame(
                           row.names = c("test_nanopore","test_pacbio")),
                       rmZeroCov = TRUE,
                       strandCollapse = FALSE,
                       verbose = FALSE,
                       output = "BSseq"
  )

  #test class
  expect_is(bs, "BSseq")
  #test dimensions
  expect_equal(dim(bs), c(1939L, 2L))
  #test rm zeroCov
  cov <- getCoverage(bs, type="Cov")
  expect_false(any(rowAlls(cov == 0)))
  #assay sums
  expect_equal(sum(getBSseq(bs, type="M")), 26632L)
  expect_equal(sum(getBSseq(bs, type="Cov")), 35708L)

  #read in files as MethylCount object (Default: rmZeroCov=T, strandCollapse=T)
  mc <- read.bedMethyl(files = infiles,
                       colData = DataFrame(
                           row.names = c("test_nanopore","test_pacbio")),
                       rmZeroCov = TRUE,
                       strandCollapse = TRUE,
                       verbose = FALSE,
                       output = "MethylCounts"
  )

  #test class
  expect_is(mc, "MethylCounts")
  #test dimensions
  expect_equal(dim(mc), c(1111L, 2L))
  #test strand collapse
  expect_true(all(strand(mc) == "*"))
  #test rm zeroCov
  cov <- getMethylCounts(mc, type="Cov")
  expect_false(any(rowAlls(cov == 0)))
  #assay sums
  expect_equal(sum(getMethylCounts(mc, type="M")), 26632L)
  expect_equal(sum(getMethylCounts(mc, type="Cov")), 35708L)

  #read in files as MethylCount object (rmZeroCov=F, strandCollapse=T)
  mc <- read.bedMethyl(files = infiles,
                       colData = DataFrame(
                           row.names = c("test_nanopore","test_pacbio")),
                       rmZeroCov = FALSE,
                       strandCollapse = TRUE,
                       verbose = FALSE,
                       output = "MethylCounts"
  )

  #test class
  expect_is(mc, "MethylCounts")
  #test dimensions
  expect_equal(dim(mc), c(1111L, 2L))
  #test strand collapse
  expect_true(all(strand(mc) == "*"))
  #test rm zeroCov
  cov <- getMethylCounts(mc, type="Cov")
  expect_false(any(rowAlls(cov == 0)))
  #assay sums
  expect_equal(sum(getMethylCounts(mc, type="M")), 26632L)
  expect_equal(sum(getMethylCounts(mc, type="Cov")), 35708L)

  #read in files as MethylCount object (rmZeroCov=T, strandCollapse=F)
  mc <- read.bedMethyl(files = infiles,
                       colData = DataFrame(
                           row.names = c("test_nanopore","test_pacbio")),
                       rmZeroCov = TRUE,
                       strandCollapse = FALSE,
                       verbose = FALSE,
                       output = "MethylCounts"
  )

  #test class
  expect_is(mc, "MethylCounts")
  #test dimensions
  expect_equal(dim(mc), c(1939L, 2L))
  #test rm zeroCov
  cov <- getMethylCounts(mc, type="Cov")
  expect_false(any(rowAlls(cov == 0)))
  #assay sums
  expect_equal(sum(getMethylCounts(mc, type="M")), 26632L)
  expect_equal(sum(getMethylCounts(mc, type="Cov")), 35708L)
})
