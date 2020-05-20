setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("hasGRanges",
         representation(gr = "GRanges"))

setMethod("seqnames", signature(x = "hasGRanges"), function(x) {
    seqnames(x@gr)
})
setReplaceMethod("seqnames", "hasGRanges", function(x, value) {
    seqnames(x@gr) <- value
    x
})

setMethod("seqlevels", signature(x = "hasGRanges"), function(x) {
    seqlevels(x@gr)
})
setReplaceMethod("seqlevels", "hasGRanges", function(x, value) {
    seqlevels(x@gr) <- value
    x
})

setMethod("seqlengths", signature(x = "hasGRanges"), function(x) {
    seqlengths(x@gr)
})
setReplaceMethod("seqlengths", "hasGRanges", function(x, value) {
    seqlengths(x@gr) <- value
    x
})

setMethod("granges", signature(x = "hasGRanges"),
    function(x, use.names = TRUE, use.mcols = FALSE, ...)
        granges(x@gr, use.names = use.names, use.mcols = use.mcols, ...)
)

## FIXME: might want a granges replacement function

setMethod("start", "hasGRanges", function(x, ...) {
    start(x@gr, ...)
})
setReplaceMethod("start", "hasGRanges", function(x, check = TRUE, value) {
    start(x@gr, check = check) <- value
    x
})

setMethod("end", "hasGRanges", function(x, ...) {
    end(x@gr, ...)
})
setReplaceMethod("end", "hasGRanges", function(x, check = TRUE, value) {
    end(x@gr, check = check) <- value
    x
})

setMethod("width", "hasGRanges", function(x) {
    width(x@gr)
})
setReplaceMethod("width", "hasGRanges", function(x, check = TRUE, value) {
    width(x@gr, check = check) <- value
    x
})

setMethod("strand", "hasGRanges", function(x) {
    strand(x@gr)
})
setReplaceMethod("strand", "hasGRanges", function(x, value) {
    strand(x@gr) <- value
    x
})

setMethod("length", "hasGRanges", function(x) length(x@gr))

setMethod("[", "hasGRanges", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x
})

setMethod("findOverlaps",
          signature(query = "hasGRanges", subject = "GenomicRanges"),
          function (query, subject, maxgap = -1L, minoverlap = 0L,
                    type = c("any", "start", "end", "within", "equal"),
                    select = c("all", "first", "last", "arbitrary"),
                    ignore.strand = FALSE, ...) {
              findOverlaps(query = query@gr, subject = subject,
                           maxgap = maxgap, minoverlap = minoverlap,
                           type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
          })

setMethod("findOverlaps",
          signature(query = "GenomicRanges", subject = "hasGRanges"),
          function (query, subject, maxgap = -1L, minoverlap = 0L,
                    type = c("any", "start", "end", "within", "equal"),
                    select = c("all", "first", "last", "arbitrary"),
                    ignore.strand = FALSE, ...) {
              findOverlaps(query = query, subject = subject@gr,
                           maxgap = maxgap, minoverlap = minoverlap,
                           type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
          })

setMethod("findOverlaps",
          signature(query = "hasGRanges", subject = "hasGRanges"),
          function (query, subject, maxgap = -1L, minoverlap = 0L,
                    type = c("any", "start", "end", "within", "equal"),
                    select = c("all", "first", "last", "arbitrary"),
                    ignore.strand = FALSE, ...) {
              findOverlaps(query = query@gr, subject = subject@gr,
                           maxgap = maxgap, minoverlap = minoverlap,
                           type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
          })

setMethod("overlapsAny", c("hasGRanges", "GenomicRanges"),
    IRanges:::default_overlapsAny
)
setMethod("overlapsAny", c("GenomicRanges", "hasGRanges"),
    IRanges:::default_overlapsAny
)
setMethod("overlapsAny", c("hasGRanges", "hasGRanges"),
    IRanges:::default_overlapsAny
)

setMethod("subsetByOverlaps", c("hasGRanges", "GenomicRanges"),
    IRanges:::default_subsetByOverlaps
)
setMethod("subsetByOverlaps", c("GenomicRanges", "hasGRanges"),
    IRanges:::default_subsetByOverlaps
)
setMethod("subsetByOverlaps", c("hasGRanges", "hasGRanges"),
    IRanges:::default_subsetByOverlaps
)

