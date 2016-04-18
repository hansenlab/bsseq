test_read.bismark_coverage <- function() {
    # Test that read.bismark() works on good input
    infile <- system.file("extdata", "test_data.fastq_bismark.bismark.cov.gz",
                          package = "bsseq")
    bsseq <- read.bismark(files = infile,
                          sampleNames = "test_data",
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          fileType = "cov",
                          verbose = FALSE)
    checkTrue(is(bsseq, "BSseq"))
    # Regression check
    # Should not fail if sampleNames have names()
    bsseq <- read.bismark(files = infile,
                          sampleNames = c(test = "test_data"),
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          fileType = "cov",
                          verbose = FALSE)
    checkTrue(is(bsseq, "BSseq"))

    # Should also work because the "cov" fileType is the same as the
    # "oldBedGraph" fileType
    bsseq <- read.bismark(files = infile,
                          sampleNames = "test_data",
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          fileType = "oldBedGraph",
                          verbose = FALSE)
    checkTrue(is(bsseq, "BSseq"))

    # Test that read.bismark() fails if given incorrect fileType.
    checkException(read.bismark(files = infile,
                                sampleNames = "test_data",
                                rmZeroCov = FALSE,
                                strandCollapse = FALSE,
                                fileType = "cytosineReport",
                                verbose = FALSE))
}

test_read.bismark_cytosineReport <- function() {
    # Test that read.bismark() works on good input
    infile <- system.file("extdata", "test_data.cytosineReport.gz",
                          package = "bsseq")
    bsseq <- read.bismark(files = infile,
                          sampleNames = "test_data",
                          rmZeroCov = FALSE,
                          strandCollapse = FALSE,
                          fileType = "cytosineReport",
                          verbose = FALSE)
    checkTrue(is(bsseq, "BSseq"))
    # Check that strandCollapse = FALSE works
    checkEquals(length(bsseq), 100L)
    checkTrue("+" %in% strand(bsseq))

    # Check that strandCollapse = TRUE works
    bsseq <- read.bismark(files = infile,
                          sampleNames = "test_data",
                          rmZeroCov = FALSE,
                          strandCollapse = TRUE,
                          fileType = "cytosineReport",
                          verbose = FALSE)
    checkEquals(length(bsseq), 50L)
    checkTrue(all(strand(bsseq) == "*"))
}

