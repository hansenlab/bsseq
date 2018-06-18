# Exported generics ------------------------------------------------------------

setGeneric(
    "findCytosines",
    function(x, context, seqlevels,
             strand = c("*", "+", "-")) standardGeneric("findCytosines"),
    signature = "x")

# Exported methods -------------------------------------------------------------

setMethod(
    "findCytosines",
    "BSgenome",
    function(x, context, seqlevels = seqlevels(x), strand = c("*", "+", "-")) {
        strand <- match.arg(strand)
        # NOTE: vmatchPattern,BSgenome-method returns a GRanges instance and
        #       automatically checks both forward and reverse strands.
        gr <- vmatchPattern(
            pattern = context,
            subject = x,
            exclude = setdiff(seqlevels(x), seqlevels),
            fixed = "subject")
        if (strand == "+") {
            gr <- gr[strand(gr) == "+"]
        } else if (strand == "-") {
            gr <- gr[strand(gr) == "-"]
        }
        # NOTE: Want just the position of the cytosine.
        resize(gr, width = 1L, fix = "start")
    }
)

# TODO: Replace with vmatchPDict() when/if this is available, where the pdict
#       argument is a DNAStringSet including the original context and its
#       reverse complement.
setMethod(
    "findCytosines",
    "DNAStringSet",
    function(x, context, seqlevels = seqlevels(x), strand = c("*", "+", "-")) {
        context <- DNAString(context)
        x <- x[seqlevels]
        strand <- match.arg(strand)
        if (strand %in% c("*", "+")) {
            fwd_gr <- as(
                vmatchPattern(
                    pattern = context,
                    subject = x,
                    fixed = "subject"),
                "GRanges")
            strand(fwd_gr) <- "+"
        } else {
            fwd_gr <- GRanges()
        }
        if (strand %in% c("*", "-")) {
            rev_gr <- as(
                vmatchPattern(
                    pattern = reverseComplement(context),
                    subject = x,
                    fixed = "subject"),
                "GRanges")
            strand(rev_gr) <- "-"
        } else {
            rev_gr <- GRanges()
        }
        gr <- c(fwd_gr, rev_gr)
        # NOTE: Want just the position of the cytosine.
        resize(gr, width = 1L, fix = "start")
    }
)

# TODOs ------------------------------------------------------------------------

# TODO: Default value of seqlevels isn't working.
