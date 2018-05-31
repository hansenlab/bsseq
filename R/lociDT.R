# Internal functions -----------------------------------------------------------

# TODO: Think through logic of 'sort' argument; where is it needed, should it
#       be part of the constructor?
.constructSeqinfoFromLociDT <- function(loci_dt, sort = TRUE) {
    unique_seqnames <- loci_dt[, as.character(unique(seqnames))]
    si <- Seqinfo(seqnames = unique_seqnames)
    if (sort) return(sortSeqlevels(si))
    si
}

.lociDTAsGRanges <- function(loci_dt, seqinfo = NULL) {
    GRanges(
        seqnames = loci_dt[["seqnames"]],
        ranges = IRanges(loci_dt[["start"]], width = 1L),
        strand = loci_dt[["strand"]],
        seqinfo = seqinfo)
}

.grAsLociDT <- function(gr) {
    # NOTE: Don't use as.data.table(gr) because it will modify gr by
    #       reference (https://github.com/Rdatatable/data.table/issues/2278)
    data.table(
        seqnames = as.factor(seqnames(gr)),
        start = start(gr),
        strand = as.factor(strand(gr)),
        key = c("seqnames", "strand", "start"))
}

.strandCollapseLociDT <- function(loci_dt, by_ref = FALSE) {
    if (loci_dt[, all(strand == "*")]) return(loci_dt)
    if (!by_ref) {
        loci_dt <- copy(loci_dt)
    }
    # Shift loci on negative strand by 1 to the left
    loci_dt[strand == "-", start := start - 1L]
    # Unstrand all loci
    loci_dt[, strand := strand("*")]
    # data.table:::funique(loci_dt)
    unique(loci_dt, by = c("seqnames", "strand", "start"))
}

.readBismarkAsLociDT <- function(file, rmZeroCov = FALSE,
                                 strandCollapse = FALSE, verbose = FALSE) {
    if (!rmZeroCov) {
        dt <- .readBismarkAsDT(
            file = file,
            col_spec = "GRanges",
            check = FALSE,
            verbose = verbose)
        if (is.null(dt[["strand"]])) {
            dt[, strand := strand("*")]
        }
        if (strandCollapse) {
            return(.strandCollapseLociDT(dt))
        } else {
            return(dt)
        }
    } else {
        # Require M and U to remove loci with zero coverage
        dt <- .readBismarkAsBSseqDT(
            file = file,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            check = FALSE,
            verbose = verbose)
        # Drop 'M' and 'U'
        dt[, c("M", "U") := .(NULL, NULL)]
    }
}

# TODO: Allow the passing of query_seqinfo and subject_seqinfo if they are
#       already available.
.findOverlaps_lociDT <- function(query, subject, maxgap = -1L, minoverlap = 0L,
                                 type = c("any", "start", "end", "within",
                                          "equal"),
                                 select = c("all", "first", "last",
                                            "arbitrary"),
                                 ignore.strand = FALSE) {
    type <- match.arg(type)
    select <- match.arg(select)

    # NOTE: Modding GenomicRanges:::findOverlaps_GNCList()
    if (!(is(query, "data.table") && is(subject, "data.table"))) {
        stop("'query' and 'subject' must be data.table objects")
    }
    if (!isTRUEorFALSE(ignore.strand)) {
        stop("'ignore.strand' must be TRUE or FALSE")
    }
    query_seqinfo <- .constructSeqinfoFromLociDT(query, FALSE)
    subject_seqinfo <- .constructSeqinfoFromLociDT(subject, FALSE)
    si <- merge(query_seqinfo, subject_seqinfo)
    q_seqlevels <- seqlevels(query_seqinfo)
    s_seqlevels <- seqlevels(subject_seqinfo)
    common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
    NG <- length(common_seqlevels)
    q_group_idx <- match(common_seqlevels, q_seqlevels)  # of length NG
    s_group_idx <- match(common_seqlevels, s_seqlevels)  # of length NG

    ## Extract 'q_groups' and 's_groups' (both of length NG).
    # q_groups <- .extract_groups_from_GenomicRanges(query)[q_group_idx]
    q_groups <- splitAsList(
        x = seq.int(0, nrow(query) - 1L),
        f = factor(query[["seqnames"]],
                   seqlevels(query_seqinfo)))[q_group_idx]
    # s_groups <- .extract_groups_from_GenomicRanges(subject)[s_group_idx]
    s_groups <- splitAsList(
        x = seq.int(0, nrow(subject) - 1L),
        f = factor(subject[["seqnames"]],
                   seqlevels(subject_seqinfo)))[s_group_idx]

    ## Extract 'nclists' and 'nclist_is_q' (both of length NG).
    ## We'll do "on-the-fly preprocessing".
    nclists <- vector(mode = "list", length = NG)
    nclist_is_q <- rep.int(NA, NG)

    ## Extract 'circle_length' (of length NG).
    circle_length <- GenomicRanges:::.get_circle_length(si)[q_group_idx]

    ## Extract 'q_space' and 's_space'.
    if (ignore.strand) {
        q_space <- s_space <- NULL
    } else {
        q_space <- as.integer(query[["strand"]]) - 3L
        s_space <- as.integer(subject[["strand"]]) - 3L
    }

    # NOTE: Modding IRanges:::NCList_find_overlaps_in_groups
    if (!(is(query, "data.table") && is(subject, "data.table"))) {
        stop("'q' and 's' must be data.table object")
    }
    if (!is(q_groups, "CompressedIntegerList")) {
        stop("'q_groups' must be a CompressedIntegerList object")
    }
    if (!is(s_groups, "CompressedIntegerList")) {
        stop("'s_groups' must be a CompressedIntegerList object")
    }

    if (!isSingleNumber(maxgap)) {
        stop("'maxgap' must be a single integer")
    }
    if (!is.integer(maxgap)) {
        maxgap <- as.integer(maxgap)
    }

    if (!isSingleNumber(minoverlap)) {
        stop("'minoverlap' must be a single integer")
    }
    if (!is.integer(minoverlap)) {
        minoverlap <- as.integer(minoverlap)
    }

    circle.length <- circle_length
    q_circle_len <- circle.length
    q_circle_len[which(nclist_is_q)] <- NA_integer_
    # NOTE: No circular ranges
    # q <- .shift_ranges_in_groups_to_first_circle(q, q_groups, q_circle_len)
    s_circle_len <- circle.length
    s_circle_len[which(!nclist_is_q)] <- NA_integer_
    # NOTE: No circular ranges
    # s <- .shift_ranges_in_groups_to_first_circle(s, s_groups, s_circle_len)

    .Call2("NCList_find_overlaps_in_groups",
           query[["start"]], query[["start"]], q_space, q_groups,
           subject[["start"]], subject[["start"]], s_space, s_groups,
           nclists, nclist_is_q,
           maxgap, minoverlap, type, select, circle.length,
           PACKAGE = "IRanges")
}

.overlapsAny_lociDT <- function(query, subject, maxgap = -1L, minoverlap = 0L,
                                type = c("any", "start", "end", "within",
                                         "equal"), ...) {
    if (is.integer(query)) {
        stop("Not yet implemented")
    }
    type <- match.arg(type)
    if (missing(subject)) {
        stop("Not yet implemented")
    }
    else {
        ahit <- .findOverlaps_lociDT(query, subject, maxgap = maxgap,
                                     minoverlap = minoverlap, type = type,
                                     select = "arbitrary", ...)
    }
    !is.na(ahit)
}

.subsetByOverlaps_lociDT <- function(x, ranges, maxgap = -1L, minoverlap = 0L,
                                     type = c("any", "start", "end", "within",
                                              "equal"),
                                     invert = FALSE, ...) {
    ov_any <- .overlapsAny_lociDT(x, ranges, maxgap = maxgap,
                                  minoverlap = minoverlap,
                                  type = match.arg(type), ...)
    if (invert)
        ov_any <- !ov_any
    x[ov_any]
}
