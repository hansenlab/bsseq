# Internal functions -----------------------------------------------------------

.strandCollapseBSseqDT <- function(bsseq_dt, by_ref = FALSE) {
    if (bsseq_dt[, all(strand == "*")]) return(bsseq_dt)
    if (!by_ref) {
        bsseq_dt <- copy(bsseq_dt)
    }
    # Shift loci on negative strand by 1 to the left
    bsseq_dt[strand == "-", start := start - 1L]
    # Unstrand all loci
    bsseq_dt[, strand := strand("*")]
    # Aggregate counts
    bsseq_dt[, .(M = sum(M), U = sum(U)), by = c("seqnames", "strand", "start")]
}


.readBismarkAsBSseqDT <- function(file, rmZeroCov = FALSE,
                                  strandCollapse = FALSE, check = TRUE,
                                  verbose = FALSE) {
    dt <- .readBismarkAsDT(
        file = file,
        col_spec = "BSseq",
        check = check,
        verbose = verbose)
    # Data is unstranded if none is provided
    # TODO: What's the data.table way to check if column exists and if create
    #       it if it doesn't exist?
    if (is.null(dt[["strand"]])) {
        dt[, strand := strand("*")]
    }
    if (strandCollapse) {
        dt <- .strandCollapseBSseqDT(dt)
    }
    if (rmZeroCov) {
        return(dt[(M + U) > 0])
    }
    dt
}

# TODOs ------------------------------------------------------------------------

# TODO: Document internal functions for my own sanity. Also, some may be useful
#       to users of bsseq (although I still won't export these for the time
#       being).
# TODO: (long term) Formalise 'lociDT' and 'BSseqDT' concepts as a S4 classes
#       (e.g., lociDT could be a concrete subclass of a GenomicRanges). Also
#       seems re-visiting the FWGenomicRanges idea (fixed-width GenomicRanges
#       class), of which 'lociDT' would be a concrete subclass.
