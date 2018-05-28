# BSseq object used in tests ---------------------------------------------------

data(BS.chr22)
bsseq_test <- BS.chr22[1:1000, ]
seqlevels(bsseq_test) <- "chr21"
bsseq_test <- suppressWarnings(rbind(bsseq_test, BS.chr22[1:1000, ]))

# Helper function used in tests ------------------------------------------------

expect_equivalent_SE <- function(SE1, SE2) {
    stopifnot(require(SummarizedExperiment))
    stopifnot(is(SE1, "SummarizedExperiment"),
              is(SE2, "SummarizedExperiment"))
    assays(SE1) <- endoapply(assays(SE1), as.matrix)
    assays(SE2) <- endoapply(assays(SE2), as.matrix)
    if (isTRUE(all.equal(SE1, SE2))) {
        return(all.equal(assays(SE1), assays(SE2)))
    }
    FALSE
}
