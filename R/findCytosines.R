# TODO: Need to update NAMESPACE
# TODO: Search both forward and reverse strands
# TODO: Default value of seqlevels as missing or NULL? Check what other
#       functions with this arg use.
findCytosines <- function(BSgenome, context = c("CG", "CA", "CC", "CT"),
                          seqlevels, verbose = getOption("verbose")) {
    stopifnot(is(BSgenome, "BSgenome"))
    context <- match.arg(context)
    context <- DNAString(context)
    seqinfo <- seqinfo(BSgenome)
    if (missing(seqlevels)) {
        seqlevels <- seqlevels(seqinfo)
    } else {
        seqinfo <- seqinfo[seqlevels]
    }
    bsparams <- new(
        "BSParams",
        X = BSgenome,
        FUN = matchPattern,
        exclude = setdiff(seqlevels(BSgenome), seqlevels))
    current_verbose <- getOption("verbose")
    options(verbose = verbose)
    on.exit(options(verbose = current_verbose), add = TRUE)
    fwd_strand <- IRangesList(
        endoapply(bsapply(bsparams, pattern = context), as, "IRanges"))
    n <- sum(as.numeric(lengths(fwd_strand)))
    fwd_strand <- .FWGRanges(
        seqnames = Rle(
            factor(names(fwd_strand), levels = seqlevels(seqinfo)),
            lengths(fwd_strand)),
        ranges = .FWIRanges(
            start = unlist(start(fwd_strand), use.names = FALSE),
            width = 1L),
        strand = Rle(strand("+"), n),
        seqinfo = seqinfo,
        elementMetadata = S4Vectors:::make_zero_col_DataFrame(n))
    rev_strand <- IRangesList(
        endoapply(
            bsapply(bsparams, pattern = reverseComplement(context)),
            as,
            "IRanges"))
    n <- sum(as.numeric(lengths(rev_strand)))
    rev_strand <- .FWGRanges(
        seqnames = Rle(
            factor(names(rev_strand), levels = seqlevels(seqinfo)),
            lengths(rev_strand)),
        ranges = .FWIRanges(
            start = unlist(start(rev_strand), use.names = FALSE),
            width = 1L),
        strand = Rle(strand("-"), n),
        seqinfo = seqinfo,
        elementMetadata = S4Vectors:::make_zero_col_DataFrame(n))
    sort(c(fwd_strand, rev_strand))
}

# TODO: findCytosines from a FASTA file (e.g., for lambda phage)
# lambda <- import("http://www.ebi.ac.uk/ena/data/view/J02459&display=fasta&download=fasta&filename=J02459.fasta")

#
# hg19_si <- suppressWarnings(keepStandardChromosomes(
#     seqinfo(BSgenome.Hsapiens.UCSC.hg19)))
# params <- new("BSParams",
#               X = BSgenome.Hsapiens.UCSC.hg19,
#               FUN = matchPattern,
#               exclude = setdiff(seqnames(BSgenome.Hsapiens.UCSC.hg19),
#                                 seqnames(hg19_si)))
# options(verbose = TRUE)
# hg19_cpgs <- IRangesList(
#     endoapply(bsapply(params, pattern = "CG"), as, "IRanges"))
# options(verbose = FALSE)
# hg19_cpgs <- GRanges(seqnames = Rle(names(hg19_cpgs), lengths(hg19_cpgs)),
#                      ranges = unlist(hg19_cpgs, use.names = FALSE),
#                      strand = "*",
#                      seqinfo = hg19_si)
# width(hg19_cpgs) <- 1L
