## findOverlaps_mclapply <- function (query, subject, maxgap = 0L, minoverlap = 1L, 
##                                    type = c("any", "start", "end", "within", "equal"),
##                                    select = c("all", "first"), ignore.strand = FALSE,
##                                    mc.cores = 1, mc.preschedule = TRUE, verbose = FALSE, ...) {
##     if(!is(query, "GenomicRanges") || !is(subject, "GenomicRanges"))
##         stop("findOverlaps_mclapply needs 'query' and 'subject' to be 'GenomicRanges'")
##     if (!isSingleNumber(maxgap) || maxgap < 0) 
##         stop("'maxgap' must be a non-negative integer")
##     type <- match.arg(type)
##     select <- match.arg(select)
##     seqinfo <- merge(seqinfo(query), seqinfo(subject))
##     DIM <- c(length(query), length(subject))
##     if (min(DIM) == 0L) {
##         matchMatrix <- matrix(integer(0), nrow = 0L, ncol = 2L, 
##                               dimnames = list(NULL, c("queryHits", "subjectHits")))
##     }
##     else {
##         querySeqnames <- seqnames(query)
##         querySplitRanges <- splitRanges(querySeqnames)
##         uniqueQuerySeqnames <- names(querySplitRanges)[sapply(querySplitRanges, length) > 0] # FIX: only keep seqnames with ranges
##         subjectSeqnames <- seqnames(subject)
##         subjectSplitRanges <- splitRanges(subjectSeqnames)
##         uniqueSubjectSeqnames <- names(subjectSplitRanges)[sapply(subjectSplitRanges, length) > 0] # FIX: only keep seqnames with ranges
##         commonSeqnames <- intersect(uniqueQuerySeqnames, 
##                                     uniqueSubjectSeqnames)
##         if (ignore.strand) {
##             queryStrand <- rep.int(1L, length(query))
##             subjectStrand <- rep.int(1L, length(subject))
##         }
##         else {
##             queryStrand <- strand(query)
##             levels(queryStrand) <- c("1", "-1", "0")
##             queryStrand@values <- as.integer(as.character(runValue(queryStrand)))
##             queryStrand <- as.vector(queryStrand)
##             subjectStrand <- strand(subject)
##             levels(subjectStrand) <- c("1", "-1", "0")
##             subjectStrand@values <- as.integer(as.character(runValue(subjectStrand)))
##             subjectStrand <- as.vector(subjectStrand)
##         }
##         queryRanges <- unname(ranges(query))
##         subjectRanges <- unname(ranges(subject))
##         matchMatrix <- do.call(rbind, mclapply(commonSeqnames, 
##                                                function(seqnm) {
##                                                    if(verbose) cat(seqnm, "\n") # FIX : added verbosity
##                                                    if (isCircular(seqinfo)[seqnm] %in% TRUE) 
##                                                        circle.length <- seqlengths(seqinfo)[seqnm]
##                                                    else circle.length <- NA
##                                                    qIdxs <- querySplitRanges[[seqnm]]
##                                                    sIdxs <- subjectSplitRanges[[seqnm]]
##                                                    ## FIX: added ::: tpo get .findOverlaps.circle
##                                                    overlaps <- GenomicRanges:::.findOverlaps.circle(circle.length, 
##                                                                                                     seqselect(queryRanges, qIdxs), seqselect(subjectRanges, 
##                                                                                                                                              sIdxs), maxgap, minoverlap, type)
##                                                    qHits <- queryHits(overlaps)
##                                                    sHits <- subjectHits(overlaps)
##                                                    matches <- cbind(queryHits = as.integer(qIdxs)[qHits], 
##                                                                     subjectHits = as.integer(sIdxs)[sHits])
##                                                    matches[which(seqselect(queryStrand, qIdxs)[qHits] * 
##                                                                  seqselect(subjectStrand, sIdxs)[sHits] != 
##                                                                  -1L), , drop = FALSE]
##                                                }, mc.cores = mc.cores, mc.preschedule = mc.preschedule))
##         if (is.null(matchMatrix)) {
##                 matchMatrix <- matrix(integer(0), nrow = 0L, 
##                                       ncol = 2L, dimnames = list(NULL, c("queryHits", 
##                                                  "subjectHits")))
##             }
##         matchMatrix <- matchMatrix[IRanges:::orderIntegerPairs(matchMatrix[, 
##                                                                            1L], matchMatrix[, 2L]), , drop = FALSE]
##     }
##     if (select == "all") {
##         new("Hits", queryHits = unname(matchMatrix[, 1L]), 
##             subjectHits = unname(matchMatrix[, 2L]), queryLength = DIM[1], 
##             subjectLength = DIM[2])
##     }
##     else {
##         IRanges:::.hitsMatrixToVector(matchMatrix, length(query))
##     }
## }
