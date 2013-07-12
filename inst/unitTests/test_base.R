test_fixes <- function() {
    gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 1))
    gr2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 2:3, width = 1))
    ov <- findOverlaps(gr1, gr2)
    checkEquals(queryHits(ov), 2)
    checkEquals(subjectHits(ov), 1)
    checkEquals(subsetByOverlaps(gr1, gr2), gr1[2])
}

test_BSseq <- function() {
    M <- matrix(1:9, 3,3)
    colnames(M) <- c("A1", "A2", "A3")
    BStest <- BSseq(pos = 1:3, chr = rep("chr1", 3), M = M, Cov = M + 2)
    checkEquals(seqnames(BStest), Rle(factor(rep("chr1", 3))))
    checkEquals(start(BStest), c(1,2,3))
    checkEquals(end(BStest), c(1,2,3))
    checkEquals(width(BStest), rep(1,3))
    checkEquals(width(BStest), rep(1,3))
    checkEquals(strand(BStest), Rle(factor(rep("*", 3), levels = c("+", "-", "*"))))
    checkEquals(dim(BStest), c(3,3))
    checkEquals(nrow(BStest), 3)
    checkEquals(ncol(BStest), 3)
    checkEquals(getCoverage(BStest, type = "M"), unname(M))
    checkEquals(getCoverage(BStest, type = "Cov"), unname(M + 2))
    checkEquals(sampleNames(BStest), colnames(M))

    BStest2 <- BSseq(pos = 3:1, chr = rep("chr1", 3), M = M[3:1,],
                     Cov = M[3:1,] + 2)
    checkEquals(BStest, BStest2)
    M2 <- rbind(M[3:1,], c(1,1,1))
    Cov2 <- rbind(M[3:1,] + 2, c(1,1,1))
    M2[3,] <- M2[3,] - 1
    Cov2[3,] <- Cov2[3,] - 1
    suppressWarnings(BStest3 <- BSseq(pos = c(3:1,1), chr = rep("chr1", 4), M = M2, Cov = Cov2))
    checkEquals(BStest, BStest3) 
}

test_overlaps <- function() {
    M <- matrix(1:9, 3,3)
    colnames(M) <- c("A1", "A2", "A3")
    Cov <- M + 2
    gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:3, width = 1))
    gr2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 3:5, width = 1))
    BS1 <- BSseq(gr = gr1, M = M, Cov = Cov)
    BS2 <- BSseq(gr = gr2, M = M, Cov = Cov)
    BSres <- BSseq(gr = gr1[3], M = M[3,,drop = FALSE], Cov = Cov[3,,drop = FALSE])
    checkEquals(length(BS1), 3L)
    checkEquals(seqlevels(BS1), seqlevels(gr1))
    checkEquals(seqlengths(BS1), seqlengths(gr1))
    checkEquals(seqnames(BS1), seqnames(gr1))
    checkEquals(granges(BS1), gr1)
    checkEquals(start(BS1), start(gr1))
    checkEquals(end(BS1), end(gr1))
    checkEquals(width(BS1), width(gr1))
    checkEquals(strand(BS1), strand(gr1))
    checkEquals(findOverlaps(BS1, BS2), findOverlaps(gr1, gr2))
    checkEquals(findOverlaps(gr1, BS2), findOverlaps(gr1, gr2))
    checkEquals(findOverlaps(BS1, gr2), findOverlaps(gr1, gr2))
    checkEquals(subsetByOverlaps(BS1, gr2), BSres)
    checkEquals(subsetByOverlaps(BS1, BS2), BSres)
    checkEquals(subsetByOverlaps(gr1, BS2), gr1[3])
}





## phenoData(BStest)
## pData(BStest)

## BStest[1:2,]
## BStest[,1:2]
## BStest[1,1]

## getMeth(BStest, type = "raw")
## getMeth(BStest, type = "raw", what = "perBase", confint = TRUE)
## getMeth(BStest, regions = gr1, type = "raw", what = "perRegion")
## getMeth(BStest, regions = gr1, type = "raw", what = "perBase")

## getCoverage(BStest, type = "Cov", what = "perBase")
## getCoverage(BStest, type = "M", what = "perBase")
## getCoverage(BStest, regions = gr1, type = "Cov", what = "perRegionAverage")
## getCoverage(BStest, regions = gr1, type = "M", what = "perRegionAverage")

## getCoverage(BStest, regions = gr1, type = "Cov", what = "perRegionTotal")
## getCoverage(BStest, regions = gr1, type = "M", what = "perRegionTotal")

## combine(BStest, BStest2)
