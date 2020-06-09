# Internal classes -------------------------------------------------------------

# The FWGRanges class is a container for the fixed-width genomic locations and
# their associated annotations.
# NOTE: The intention is to make this a fully-fledged GenomicRanges subclass
#       that is part of the GenomicRRanges package. For now, this class is
#       internal to bsseq and only a subset of methods are properly implemented.
.FWGRanges <- setClass(
    "FWGRanges",
    contains = "GenomicRanges",
    representation(
        seqnames = "Rle",
        ranges = "FWIRanges",
        strand = "Rle",
        seqinfo = "Seqinfo"
    ),
    prototype(
        seqnames = Rle(factor()),
        strand = Rle(strand())
    )
)

# Internal methods -------------------------------------------------------------

# NOTE: Combine the new "parallel slots" with those of the parent class. Make
#       sure to put the new parallel slots **first**. See R/Vector-class.R file
#       in the S4Vectors package for what slots should or should not be
#       considered "parallel".
setMethod(
    "parallel_slot_names",
    "FWGRanges",
    function(x) {
        c(GenomicRanges:::extraColumnSlotNames(x), "seqnames", "ranges",
          "strand", callNextMethod())
    }
)

setMethod("seqnames", "FWGRanges", function(x) x@seqnames)

setMethod("strand", "FWGRanges", function(x) x@strand)

setMethod("seqinfo", "FWGRanges", function(x) x@seqinfo)

setMethod(
    "ranges",
    "FWGRanges",
    function(x, use.names = TRUE, use.mcols = FALSE) {
        if (!isTRUEorFALSE(use.names)) stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols)) stop("'use.mcols' must be TRUE or FALSE")
        # TODO: This is a work-around the requirement in the GenomicRanges
        #       validity method that ranges(GenomicRanges) returns a IRanges or
        #       IPos instance.
        # ans <- x@ranges
        ans <- as(x@ranges, "IRanges")
        if (!use.names) names(ans) <- NULL
        if (use.mcols) mcols(ans) <- mcols(x)
        ans
    }
)

# TODO: These wouldn't be necessary if ranges(FWGRanges) returned a FWIRanges
#       instead of an IRanges.
setMethod("start", "FWGRanges", function(x) start(x@ranges))
setMethod("width", "FWGRanges", function(x) width(x@ranges))
setMethod("end", "FWGRanges", function(x) {
    if (x@ranges@width == 1L) {
        start(x)
    } else {
        width(x) - 1L + start(x)
    }
})

# TODO: Follow shift,GRanges-method
setMethod("shift", "FWGRanges", function(x, shift = 0L, use.names = TRUE) {
    stopifnot(use.names)
    new_ranges <- shift(x@ranges, shift, use.names)
    x@ranges <- new_ranges
    validObject(x)
    x
})

# NOTE: Needed for seqinfo<-,FWGRanges-method
setMethod("update", "FWGRanges", function(object, ...) {
    BiocGenerics:::replaceSlots(object, ...)
})

# NOTE: This is pieced together in order to get something up and running. It
#       should, however, be possible to this 'properly' via inheritance. Most
#       of this hack is to work around the fact that ranges(FWGRanges)
#       currently has to return a IRanges instance instead of a FWIRanges
#       instance.
.findOverlaps_FWGRanges <- function(query, subject, maxgap = -1L,
                                    minoverlap = 0L,
                                    type = c("any", "start", "end", "within",
                                             "equal"),
                                    select = c("all", "first", "last",
                                               "arbitrary"),
                                    ignore.strand = FALSE) {
    # findOverlaps,GenomicRanges,GenomicRanges-method ----------------------
    type <- match.arg(type)
    select <- match.arg(select)

    # GenomicRanges:::findOverlaps_GNCList() -------------------------------
    if (!(is(query, "FWGRanges") && is(subject, "FWGRanges"))) {
        stop("'query' and 'subject' must be FWGRanges objects")
    }
    if (!isTRUEorFALSE(ignore.strand)) {
        stop("'ignore.strand' must be TRUE or FALSE")
    }

    si <- merge(seqinfo(query), seqinfo(subject))
    q_seqlevels <- seqlevels(query)
    s_seqlevels <- seqlevels(subject)
    common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
    NG <- length(common_seqlevels)
    q_group_idx <- match(common_seqlevels, q_seqlevels)
    s_group_idx <- match(common_seqlevels, s_seqlevels)

    # Extract 'q_groups' and 's_groups' (both of length NG).
    q_groups <- GenomicRanges:::.extract_groups_from_GenomicRanges(
        query)[q_group_idx]
    s_groups <- GenomicRanges:::.extract_groups_from_GenomicRanges(
        subject)[s_group_idx]

    # Extract 'nclists' and 'nclist_is_q' (both of length NG).
    if (is(subject, "GNCList")) {
        nclists <- subject@nclists[s_group_idx]
        nclist_is_q <- rep.int(FALSE, NG)
    } else if (is(query, "GNCList")) {
        nclists <- query@nclists[q_group_idx]
        nclist_is_q <- rep.int(TRUE, NG)
    } else {
        # We'll do "on-the-fly preprocessing".
        nclists <- vector(mode = "list", length = NG)
        nclist_is_q <- rep.int(NA, NG)
    }

    # Extract 'circle_length' (of length NG).
    circle_length <- GenomicRanges:::.get_circle_length(si)[q_group_idx]

    ##Extract 'q_space' and 's_space'.
    if (ignore.strand) {
        q_space <- s_space <- NULL
    } else {
        q_space <- as.integer(strand(query)) - 3L
        s_space <- as.integer(strand(subject)) - 3L
    }

    # IRanges:::find_overlaps_in_groups_NCList() ---------------------------
    q <- query@ranges
    s <- subject@ranges
    circle.length <- circle_length
    if (!(is(q, "IntegerRanges") && is(s, "IntegerRanges"))) {
        stop("'q' and 's' must be IntegerRanges object")
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

    q_circle_len <- circle.length
    q_circle_len[which(nclist_is_q)] <- NA_integer_
    q <- IRanges:::.shift_ranges_in_groups_to_first_circle(
        x = q,
        x_groups = q_groups,
        circle.length = q_circle_len)
    s_circle_len <- circle.length
    s_circle_len[which(!nclist_is_q)] <- NA_integer_
    s <- IRanges:::.shift_ranges_in_groups_to_first_circle(
        x = s,
        x_groups = s_groups,
        circle.length = s_circle_len)
    .Call2("C_find_overlaps_in_groups_NCList",
           start(q), end(q), q_space, q_groups,
           start(s), end(s), s_space, s_groups,
           nclists, nclist_is_q,
           maxgap, minoverlap, type, select, circle.length,
           PACKAGE = "IRanges")
}

setMethod("findOverlaps", c("FWGRanges", "FWGRanges"), .findOverlaps_FWGRanges)

# TODO: A dedicated duplicated method would be faster than
#       duplicated,GenomicRanges, because it could use end(x) instead of
#       width(x) in call to duplicatedIntegerQuads().

# Internal functions -----------------------------------------------------------

# TODO: Warn if x contains mix of + and * or - and * loci?
.strandCollapse <- function(x) {
    stopifnot(is(x, "GenomicRanges"))
    if (all(strand(x) == "*")) return(x)
    # Shift loci on negative strand by 1 to the left
    x[strand(x) == "-"] <- IRanges::shift(x[strand(x) == "-"], -1L)
    # Unstrand all loci
    x <- unstrand(x)
    # Extract unique loci
    unique(x)
}

# TODO: Document that the default 'sort = TRUE' applies sort(sortSeqlevels())
#       to the output. This is the default behaviour because it results in
#       the smallest returned object (albeit at the small cost of a sort).
.readBismarkAsFWGRanges <- function(file, rmZeroCov = FALSE,
                                    strandCollapse = FALSE, sort = TRUE,
                                    nThread = 1L, verbose = FALSE) {
    # Argument checks ----------------------------------------------------------

    stopifnot(isTRUEorFALSE(rmZeroCov))
    stopifnot(isTRUEorFALSE(strandCollapse))
    stopifnot(isTRUEorFALSE(sort))

    # Quieten R CMD check about 'no visible binding for global variable'
    M <- U <- NULL

    # Read file to construct data.table of valid loci --------------------------
    if (rmZeroCov) {
        dt <- .readBismarkAsDT(
            file = file,
            col_spec = "BSseq",
            check = TRUE,
            verbose = verbose)
        if (strandCollapse && !is.null(dt[["strand"]]) &&
            !dt[, all(strand == "*")]) {
            # Shift loci on negative strand by 1 to the left and then remove
            # strand since no longer valid.
            dt[strand == "-", start := start - 1L][, strand := NULL]
            # Aggregate counts at loci with the same 'seqnames' and 'start'.
            dt <- dt[,
                     list(M = sum(M), U = sum(U)), by = c("seqnames", "start")]
        }
        # Identify loci with non-zero coverage then drop 'M' and 'U' as no
        # longer required.
        dt <- dt[(M + U) > 0][, c("M", "U") := list(NULL, NULL)]
    } else {
        dt <- .readBismarkAsDT(
            file = file,
            col_spec = "GRanges",
            check = FALSE,
            nThread = nThread,
            verbose = verbose)
        if (strandCollapse && !is.null(dt[["strand"]]) &&
            !dt[, all(strand == "*")]) {
            # Shift loci on negative strand by 1 to the left and then remove
            # strand since no longer valid.
            dt[strand == "-", start := start - 1L][, strand := NULL]
            dt <- data.table:::funique(dt)
        }
    }

    # Construct FWGRanges from 'dt' --------------------------------------------

    # NOTE: Sorting results in a smaller FWGRanges object because the
    #       'seqnames' and 'strand' slots are more compressible in their Rle
    #       representation.
    if (sort) {
        if (is.null(dt[["strand"]])) {
            setkey(dt, seqnames, start)
        } else {
            setkey(dt, seqnames, strand, start)
        }
    }
    seqnames <- Rle(dt[["seqnames"]])
    dt[, seqnames := NULL]
    seqinfo <- Seqinfo(seqnames = levels(seqnames))
    ranges <- .FWIRanges(start = dt[["start"]], width = 1L)
    dt[, start := NULL]
    mcols <- make_zero_col_DFrame(length(ranges))
    if (is.null(dt[["strand"]])) {
        strand <- strand(Rle("*", length(seqnames)))
    } else {
        strand <- Rle(dt[["strand"]])
        dt[, strand := NULL]
    }
    fwgranges <- .FWGRanges(
        seqnames = seqnames,
        ranges = ranges,
        strand = strand,
        seqinfo = seqinfo,
        elementMetadata = mcols)
    # NOTE: Final sort is to re-order with respect to sorted seqlevels.
    if (sort) {
        fwgranges <- sort(sortSeqlevels(fwgranges))
    }
    fwgranges
}

# TODO: Document that this applies sort(sortSeqlevels()) to the output. It is
#       deliberate that there is no option to override this behaviour.
.contructFWGRangesFromBismarkFiles <- function(files,
                                               rmZeroCov,
                                               strandCollapse,
                                               verbose,
                                               nThread,
                                               BPPARAM) {
    subverbose <- max(as.integer(verbose) - 1L, 0L)

    # TODO: Instead of using the 'largest' file, use the largest
    #       'cytosine report' file, which will have all loci in the
    #       reference genome; provided all samples were aligned to the same
    #       reference genome, this means it contains all loci.
    # TODO: Initialise using the 'largest' file (i.e. largest number of lines)?
    #       Would like to do this without reading the data into memory.
    #       Some benchmarks can be found at
    #       https://gist.github.com/peterhurford/0d62f49fd43b6cf078168c043412f70a
    #       My initial tests using /users/phickey/GTExScripts/FlowSortingProject/hdf5/extdata/methylation/nonCG/5248_BA9_neg_CHG_report.txt (32 GB) give:
    #       wc -l:                                       77.000s
    #       R.utils::readLines():                      1165.299s
    #       nrow(fread(..., select = 1, nThread = 1)):  582.721s (359s re-run)
    #       nrow(fread(..., select = 1, nThread = 10)):  82.029s
    #       nrow(fread(..., select = 1, nThread = 40)):  81.408s
    #       file.size():                                  0.000s
    #       Of course, fread() only works directly with non-[b]gzipped files.
    #       And subsequent runs of fread() benefit from some cacheing effect
    #       that I don't fully understand except to know that subsequent runs
    #       are 'artificially' faster.
    #       And using file.size() will be innaccurate if files are a mix of
    #       compressed and uncompressed files.
    # Initalise `loci_dt` using the first file.
    if (verbose) {
        message("[.contructFWGRangesFromBismarkFiles] Extracting loci from ",
                "'", files[1L], "'")
    }
    loci_from_first_file <- .readBismarkAsFWGRanges(
        file = files[[1L]],
        rmZeroCov = rmZeroCov,
        strandCollapse = strandCollapse,
        nThread = nThread,
        verbose = subverbose)
    # Identify loci not found in first file.
    # TODO: Pre-process loci as a GNCList?
    # Set number of tasks to ensure the progress bar gives frequent updates.
    # NOTE: The progress bar increments once per task
    #       (https://github.com/Bioconductor/BiocParallel/issues/54).
    #       Although it is somewhat of a bad idea to overrides a user-specified
    #       bptasks(BPPARAM), the value of bptasks(BPPARAM) doesn't affect
    #       performance in this instance, and so we opt for a useful progress
    #       bar. Only SnowParam (and MulticoreParam by inheritance) have a
    #       bptasks<-() method.
    if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
        bptasks(BPPARAM) <- length(files) - 1L
    }
    list_of_loci_from_other_files_not_in_first_file <- bplapply(
        files[-1L], function(file, loci_from_first_file) {
            # Read this file.
            loci_from_this_file <- .readBismarkAsFWGRanges(
                file = file,
                rmZeroCov = rmZeroCov,
                strandCollapse = strandCollapse,
                verbose = subverbose)
            subsetByOverlaps(
                x = loci_from_this_file,
                ranges = loci_from_first_file,
                type = "equal",
                invert = TRUE)
        }, loci_from_first_file = loci_from_first_file,
        BPPARAM = BPPARAM)
    # Identify unique FWGRanges.
    loci_non_found_in_first_file <- unique(
        do.call(c, list_of_loci_from_other_files_not_in_first_file))
    loci <- c(loci_from_first_file, loci_non_found_in_first_file)

    # Sort the loci
    sort(sortSeqlevels(loci))
}

# TODOs ------------------------------------------------------------------------

# TODO: Document internal classes, methods, and functions for my own sanity.
#       Also, some may be useful to users of bsseq (although I still won't
#       export these for the time being). Ultimately, I'd like to promote
#       FWGRanges and FWIRanges to 'official' GenomicRanges and IntegerRanges
#       subclasses.
