data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                     "DataFrame")
        elementMetadata(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

.checkAssayNames <- function(object, names) {
    nms <- names(assays(object, withDimnames = FALSE))
    if(!all(names %in% nms))
        return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                       class(object), paste0(names, collapse = ", ")))
    else
        NULL
}
