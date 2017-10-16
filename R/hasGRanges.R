setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("hasGRanges",
         representation(gr = "GRanges"))

setMethod("seqnames", signature(x = "hasGRanges"), function(x) {
    seqnames(x@gr)
})
setReplaceMethod("seqnames", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqnames(gr) <- value
    x@gr <- gr
    x
})

setMethod("seqlevels", signature(x = "hasGRanges"), function(x) {
    seqlevels(x@gr)
})
setReplaceMethod("seqlevels", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqlevels(gr) <- value
    x@gr <- gr
    x
})

setMethod("seqlengths", signature(x = "hasGRanges"), function(x) {
    seqlengths(x@gr)
})
setReplaceMethod("seqlengths", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqlengths(gr) <- value
    x@gr <- gr
    x
})

setMethod("granges", signature(x = "hasGRanges"),
          function(x) x@gr)

## FIXME: might want a granges replacement function

setMethod("start", "hasGRanges", function(x, ...) {
    start(x@gr, ...)
})
setReplaceMethod("start", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    start(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("end", "hasGRanges", function(x, ...) {
    end(x@gr, ...)
})
setReplaceMethod("end", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    end(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("width", "hasGRanges", function(x) {
    width(x@gr)
})
setReplaceMethod("width", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    width(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("strand", "hasGRanges", function(x) {
    strand(x@gr)
})
setReplaceMethod("strand", "hasGRanges", function(x, value) {
    gr <- granges(x)
    strand(gr) <- value
    x@gr <- gr
    x
})

setMethod("length", "hasGRanges", function(x) length(x@gr))

## setMethod("subsetByOverlaps",
##           signature(query = "hasGRanges", subject = "GenomicRanges"),
##           function(query, subject, maxgap = 0L, minoverlap = 1L,
##                    type = c("any", "start", "end", "within", "equal"),
##                    ignore.strand = FALSE, ...) {
##               ov <- findOverlaps(query = granges(query), subject = subject,
##                                  maxgap = maxgap, minoverlap = minoverlap,
##                                  type = match.arg(type), select = "arbitrary",
##                                  ignore.strand = ignore.strand, ... )
##               query[!is.na(ov)]
##           })

## setMethod("subsetByOverlaps",
##           signature(query = "hasGRanges", subject = "hasGRanges"),
##           function(query, subject, maxgap = 0L, minoverlap = 1L,
##                    type = c("any", "start", "end", "within", "equal"),
##                    ignore.strand = FALSE, ...) {
##               ov <- findOverlaps(query = granges(query), subject = granges(subject),
##                                  maxgap = maxgap, minoverlap = minoverlap,
##                                  type = match.arg(type), select = "arbitrary",
##                                  ignore.strand = ignore.strand, ... )
##               query[!is.na(ov)]
##           })

## setMethod("subsetByOverlaps",
##           signature(query = "GenomicRanges", subject = "hasGRanges"),
##           function(query, subject, maxgap = 0L, minoverlap = 1L,
##                    type = c("any", "start", "end", "within", "equal"),
##                    ignore.strand = FALSE, ...) {
##               ov <- findOverlaps(query = query, subject = granges(subject),
##                                  maxgap = maxgap, minoverlap = minoverlap,
##                                  type = match.arg(type), select = "arbitrary",
##                                  ignore.strand = ignore.strand, ... )
##               query[!is.na(ov)]
##           })

setMethod("findOverlaps",
          signature(query = "hasGRanges", subject = "GenomicRanges"),
          function (query, subject, maxgap = -1L, minoverlap = 0L,
                    type = c("any", "start", "end", "within", "equal"),
                    select = c("all", "first", "last", "arbitrary"),
                    ignore.strand = FALSE, ...) {
              findOverlaps(query = granges(query), subject = subject,
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
              findOverlaps(query = granges(query), subject = granges(subject),
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
              findOverlaps(query = query, subject = granges(subject),
                           maxgap = maxgap, minoverlap = minoverlap,
                           type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
          })

setMethod("[", "hasGRanges", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x
})

