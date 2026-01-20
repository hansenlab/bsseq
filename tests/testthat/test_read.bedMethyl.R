context("read.bedMethyl")

test_that("read.bedMethyl() works for bedMethyl files", {
    #load test data
    infiles <- c(system.file("extdata/HG002_nanopore_test.bedMethyl.gz",
    package = "bsseq"),
               system.file("extdata/HG002_pacbio_test.bedMethyl.gz",
    package = "bsseq"))

  #read in files as bsseq object (Default: strandCollapse=T)
  bs <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
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
  #assay sums
  expect_equal(sum(getBSseq(bs, type="M")), 26632L)
  expect_equal(sum(getBSseq(bs, type="Cov")), 35708L)

  #read in files as bsseq object (strandCollapse=F)
  bs <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
     strandCollapse = FALSE,
     verbose = FALSE,
     output = "BSseq"
     )

  #test class
  expect_is(bs, "BSseq")
  #test dimensions
  expect_equal(dim(bs), c(1939L, 2L))
  #assay sums
  expect_equal(sum(getBSseq(bs, type="M")), 26632L)
  expect_equal(sum(getBSseq(bs, type="Cov")), 35708L)

  #read in files as MethylCount object (Default strandCollapse=T)
  mc <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
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
  #assay sums
  expect_equal(sum(getMethylCounts(mc, type="M")), 26632L)
  expect_equal(sum(getMethylCounts(mc, type="Cov")), 35708L)

  #read in files as MethylCount object (strandCollapse=T)
  mc <- read.bedMethyl(files = infiles,
     colData = DataFrame(row.names = c("test_nanopore","test_pacbio")),
     strandCollapse = FALSE,
     verbose = FALSE,
     output = "MethylCounts"
     )

  #test class
  expect_is(mc, "MethylCounts")
  #test dimensions
  expect_equal(dim(mc), c(1939L, 2L))
  #assay sums
  expect_equal(sum(getMethylCounts(mc, type="M")), 26632L)
  expect_equal(sum(getMethylCounts(mc, type="Cov")), 35708L)
})
