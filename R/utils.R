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
        mcols(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

.checkAssayNames <- function(object, names) {
    nms <- assayNames(object)
    if(!all(names %in% nms))
        return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                       class(object), paste0(names, collapse = ", ")))
    else
        NULL
}

.checkAssayClasses <- function(object, names) {
    nms <- intersect(assayNames(object), names)
    is_DelayedMatrix <- vapply(assays(object, withDimnames = FALSE)[nms],
                               function(assay) {
                                   if (is.null(assay)) {
                                       return(TRUE)
                                   }
                                   is(assay, "DelayedMatrix")
                               }, logical(1L))
    if (!all(is_DelayedMatrix)) {
        return(paste0("assay slots '", paste0(nms, collapse = "', '"),
                      "' of object of class '", class(object),
                      "' need be DelayedMatrix objects"))
    } else {
        NULL
    }
}

.oldTrans <- function(x) {
    y <- x
    ix <- which(x < 0)
    ix2 <- which(x > 0)
    y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
    y[ix2] <- 1/(1 + exp(-x[ix2]))
    y
}
