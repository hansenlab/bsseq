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

# NOTE: Combine the new parallel slots with those of the parent class. Make
#       sure to put the new parallel slots *first*.
setMethod(
    "parallelSlotNames",
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

    # IRanges:::NCList_find_overlaps_in_groups() ---------------------------
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
    .Call2("NCList_find_overlaps_in_groups",
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

.strandCollapseFWGRanges <- function(fwgranges) {
    if (all(strand(fwgranges) == "*")) return(fwgranges)
    # Shift loci on negative strand by 1 to the left
    fwgranges[strand(fwgranges) == "-"] <- shift(
        x = fwgranges[strand(fwgranges) == "-"],
        shift = -1L)
    # Unstrand all loci
    fwgranges <- unstrand(fwgranges)
    # Extract unique loci
    unique(fwgranges)
}

# TODO: Document that the default 'sort = TRUE' applies sort(sortSeqlevels())
#       to the output. This is the default behaviour because it results in
#       the smallest returned object (albeit at the small cost of a sort).
.readBismarkAsFWGRanges <- function(file, rmZeroCov = FALSE,
                                    strandCollapse = FALSE, sort = TRUE,
                                    verbose = FALSE) {
    if (!rmZeroCov) {
        dt <- .readBismarkAsDT(
            file = file,
            col_spec = "GRanges",
            check = FALSE,
            verbose = verbose)
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
        mcols <- S4Vectors:::make_zero_col_DataFrame(length(ranges))
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
        if (strandCollapse) {
            fwgranges <- .strandCollapseFWGRanges(fwgranges)
        }
    } else {
        # Require M and U to remove loci with zero coverage
        dt <- .readBismarkAsBSseqDT(
            file = file,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            check = FALSE,
            verbose = verbose)
        # Drop 'M' and 'U' as no longer required.
        dt[, c("M", "U") := .(NULL, NULL)]
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
        mcols <- S4Vectors:::make_zero_col_DataFrame(length(ranges))
        if (is.null(dt[["strand"]])) {
            strand <- strand(Rle("*", nrow(dt)))
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
    }
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
                                               BPPARAM) {
    subverbose <- max(as.integer(verbose) - 1L, 0L)

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
    fwgranges_from_first_file <- .readBismarkAsFWGRanges(
        file = files[[1L]],
        rmZeroCov = rmZeroCov,
        strandCollapse = strandCollapse,
        verbose = subverbose)
    # Identify loci not found in first file.
    # TODO: Pre-process fwgranges as a GNCList?
    list_of_fwgranges_from_other_files <- bplapply(
        files[-1L], function(file, fwgranges_from_first_file) {
            # TODO: This message won't appear in main process, so probably remove.
            if (verbose) {
                message("[.contructFWGRangesFromBismarkFiles] Extracting loci ",
                        "from '", file, "'")
            }
            # Read this file.
            fwgranges_from_this_file <- .readBismarkAsFWGRanges(
                file = file,
                rmZeroCov = rmZeroCov,
                strandCollapse = strandCollapse,
                verbose = subverbose)
            subsetByOverlaps(
                x = fwgranges_from_this_file,
                ranges = fwgranges_from_first_file,
                type = "equal",
                invert = TRUE)
        }, fwgranges_from_first_file = fwgranges_from_first_file,
        BPPARAM = BPPARAM)
    # Identify unique FWGRanges.
    fwgranges_non_found_in_first_file <- unique(
        do.call(c, list_of_fwgranges_from_other_files))
    fwgranges <- c(fwgranges_from_first_file, fwgranges_non_found_in_first_file)

    # Sort the loci
    sort(sortSeqlevels(fwgranges))
}

.constructCountsFromBismarkFileAndFWGRanges <- function(b, files,
                                                        strandCollapse,
                                                        fwgranges, grid,
                                                        M_sink, Cov_sink,
                                                        M_sink_lock,
                                                        Cov_sink_lock,
                                                        verbose, BPPARAM) {
    # Read in data for b-th file -----------------------------------------------
    file <- files[b]
    if (verbose) {
        message("[.constructCountsFromBismarkFileAndFWGRanges] Extracting ",
                "counts from ", "'", file, "'")
    }
    bsseq_dt <- .readBismarkAsBSseqDT(
        file = file,
        # NOTE: Don't remove loci with zero coverage (it's unnecessary).
        rmZeroCov = FALSE,
        strandCollapse = strandCollapse,
        check = TRUE,
        verbose = verbose)
    seqnames <- Rle(bsseq_dt[["seqnames"]])
    bsseq_dt[, seqnames := NULL]
    seqinfo <- Seqinfo(seqnames = levels(seqnames))
    ranges <- .FWIRanges(start = bsseq_dt[["start"]], width = 1L)
    bsseq_dt[, start := NULL]
    mcols <- S4Vectors:::make_zero_col_DataFrame(length(ranges))
    if (is.null(bsseq_dt[["strand"]])) {
        strand <- strand(Rle("*", length(seqnames)))
    } else {
        strand <- Rle(bsseq_dt[["strand"]])
        bsseq_dt[, strand := NULL]
    }
    this_sample_fwgranges <- .FWGRanges(
        seqnames = seqnames,
        ranges = ranges,
        strand = strand,
        seqinfo = seqinfo,
        elementMetadata = mcols)
    M <- matrix(rep(0L, length(fwgranges)), ncol = 1)
    Cov <- matrix(rep(0L, length(fwgranges)), ncol = 1)
    ol <- findOverlaps(this_sample_fwgranges, fwgranges)
    M[subjectHits(ol)] <- bsseq_dt[["M"]][queryHits(ol)]
    Cov[subjectHits(ol)] <- bsseq_dt[, Cov := M + U][["Cov"]][queryHits(ol)]

    # Return M and U or write them to the RealizationSink objects --------------

    if (is.null(M_sink)) {
        return(list(M = M, Cov = Cov))
    }
    # Write to M_sink and Cov_sink while respecting the IPC locks.
    ipclock(M_sink_lock)
    write_block_to_sink(M, M_sink, grid[[b]])
    ipcunlock(M_sink_lock)
    ipclock(Cov_sink_lock)
    write_block_to_sink(Cov, Cov_sink, grid[[b]])
    ipcunlock(Cov_sink_lock)
    NULL
}

.constructCountsFromBismarkFilesAndFWGRanges <- function(files,
                                                         fwgranges,
                                                         strandCollapse,
                                                         verbose,
                                                         BPPARAM,
                                                         BACKEND) {
    # Set up ArrayGrid so that each block contains data for a single sample.
    ans_nrow <- length(fwgranges)
    ans_ncol <- length(files)
    grid <- RegularArrayGrid(c(ans_nrow, ans_ncol), c(ans_nrow, 1L))
    # Construct RealizationSink objects.
    if (is.null(BACKEND)) {
        M_sink <- NULL
        Cov_sink <- NULL
        M_sink_lock <- NULL
        Cov_sink_lock <- NULL
    } else {
        M_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(Cov_sink), add = TRUE)
        M_sink_lock <- ipcid()
        on.exit(ipcremove(M_sink_lock), add = TRUE)
        Cov_sink_lock <- ipcid()
        on.exit(ipcremove(Cov_sink_lock), add = TRUE)
    }
    # Set number of tasks to ensure the progress bar gives frequent updates.
    # NOTE: The progress bar increments once per task
    #       (https://github.com/Bioconductor/BiocParallel/issues/54).
    #       Although it is somewhat of a bad idea to overrides a user-specified
    #       bptasks(BPPARAM), the value of bptasks(BPPARAM) doesn't affect
    #       performance in this instance, and so we opt for a useful progress
    #       bar. Only SnowParam (and MulticoreParam by inheritance) have a
    #       bptasks<-() method.
    # TODO: Check that setting number of tasks doesn't affect things (e.g.,
    #       the cost of transfering loci_dt to the workers may be substantial).
    if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
        bptasks(BPPARAM) <- length(grid)
    }
    counts <- bptry(bplapply(
        X = seq_along(grid),
        FUN = .constructCountsFromBismarkFileAndFWGRanges,
        files = files,
        strandCollapse = strandCollapse,
        fwgranges = fwgranges,
        grid = grid,
        M_sink = M_sink,
        Cov_sink = Cov_sink,
        M_sink_lock = M_sink_lock,
        Cov_sink_lock = Cov_sink_lock,
        verbose = verbose,
        BPPARAM = BPPARAM))
    if (!all(bpok(counts))) {
        # TODO: This isn't yet properly impelemented
        # TODO: Feels like stop() rather than warning() should be used, but
        #       stop() doesn't allow for the return of partial results;
        #       see https://support.bioconductor.org/p/109374/
        warning("read.bismark() encountered errors: ",
                sum(!bpok(counts)), " of ", length(counts),
                " files failed.\n",
                "read.bismark() has returned partial results, including ",
                "errors, for debugging purposes.\n",
                "It may be possible to re-run just these failed files.\n",
                "See help(\"read.bismark\")",
                call. = FALSE)
        # NOTE: Return intermediate results as well as all derived variables.
        return(list(counts = counts,
                    M_sink = M_sink,
                    Cov_sink = Cov_sink,
                    BACKEND = BACKEND))
    }
    # Construct M and Cov from results of
    # .constructCountsFromBismarkFileAndLociDT().
    if (is.null(BACKEND)) {
        # Returning matrix objects.
        M <- do.call(c, lapply(counts, "[[", "M"))
        attr(M, "dim") <- c(ans_nrow, ans_ncol)
        Cov <- do.call(c, lapply(counts, "[[", "Cov"))
        attr(Cov, "dim") <- c(ans_nrow, ans_ncol)
    } else {
        # Returning DelayedMatrix objects.
        M <- as(M_sink, "DelayedArray")
        Cov <- as(Cov_sink, "DelayedArray")
    }

    return(list(M = M, Cov = Cov))
}
