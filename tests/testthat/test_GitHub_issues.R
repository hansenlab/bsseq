context("GitHub issues")

test_that("Issue 74 is fixed", {

    # Can combine two BSseq objects with same loci but different names.
    # (https://github.com/hansenlab/bsseq/issues/74#issue-341014110)
    nloci <- 68164
    nsamples <- 12
    M <- matrix(rpois(nloci * nsamples, 10), ncol = nsamples)
    Cov <- M + rpois(nloci * nsamples, 2)
    gr <- sort(GRanges(
        seqnames = sample(1:10, nloci, replace = TRUE),
        ranges = IRanges(sample(1:10^7, nloci, replace = FALSE), width = 1)))
    BS <- BSseq(M = M, Cov = Cov, gr = gr, sampleNames = paste0("V", 1:nsamples))
    BS.avg <- collapseBSseq(BS, group = paste0("W", c(rep(1, 5), 2:8)))
    BS.all <- combine(BS, BS.avg)
    colnames(BS.all)

    # As above, but when one of the BSseq objects contains additional colData.
    # https://github.com/hansenlab/bsseq/issues/74#issue-341014110
    BS2 <- BS
    BS2$rep <- paste0("A", c(rep(1, 5), 2:8))
    BS2.avg <- collapseBSseq(BS, group = paste0("W", c(rep(1, 5), 2:8)))
    BS2.all <- combine(BS2,BS2.avg)
    colnames(BS2.all)
    colData(BS2.all)
})
