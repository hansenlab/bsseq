# Exported functions -----------------------------------------------------------

findLoci <- function(pattern, subject, include = seqlevels(subject),
                     strand = c("*", "+", "-"), fixed = "subject",
                     resize = TRUE) {
    pattern <- DNAString(pattern)
    strand <- match.arg(strand)

    if (is.character(subject) || is(subject, "RTLFile")) {
        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("rtracklayer package required for importing 'subject'.")
        }
        subject <- try(rtracklayer::import(subject), silent = TRUE)
        if (!is(subject, "DNAStringSet")) {
            stop("Unable to import 'subject' as a DNAStringSet.")
        }
    }

    if (is(subject, "BSgenome")) {
        # NOTE: vmatchPattern,BSgenome-method returns a GRanges instance and
        #       automatically checks both forward and reverse strands.
        gr <- vmatchPattern(
            pattern = pattern,
            subject = subject,
            # NOTE: This must be a regular expression;
            #       see https://github.com/Bioconductor/BSgenome/issues/1.
            exclude = paste0("^", setdiff(seqnames(subject), include), "$"),
            fixed = fixed)
        if (identical(strand, "+")) {
            gr <- gr[strand(gr) == "+"]
        } else if (identical(strand, "-")) {
            gr <- gr[strand(gr) == "-"]
        }
    } else if (is(subject, "DNAStringSet")) {
        subject <- subject[include]
        if (strand %in% c("*", "+")) {
            fwd_gr <- as(
                vmatchPattern(
                    pattern = pattern,
                    subject = subject,
                    fixed = fixed),
                "GRanges")
            strand(fwd_gr) <- "+"
        } else {
            fwd_gr <- GRanges()
        }
        if (strand %in% c("*", "-")) {
            rev_gr <- as(
                vmatchPattern(
                    pattern = reverseComplement(pattern),
                    subject = subject,
                    fixed = fixed),
                "GRanges")
            strand(rev_gr) <- "-"
        } else {
            rev_gr <- GRanges()
        }
        gr <- c(fwd_gr, rev_gr)
    } else {
        stop("Cannot handle 'subject' of class ", class(subject), ".")
    }

    if (resize) {
        gr <- resize(gr, width = 1L, fix = "start")
    }
    gr
}

# TODOs ------------------------------------------------------------------------

# TODO: Default value of seqlevels may not be working.
# TODO: What happens if subject is a non-filepath character (e.g., "CATGCG") or
#       a DNAString?
# TODO: Allow passing of seqinfo and/or autogenerate?
