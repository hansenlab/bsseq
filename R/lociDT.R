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
        dt
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

.contructLociDTAndSeqinfoFromBismarkFiles <- function(files,
                                                      rmZeroCov,
                                                      strandCollapse,
                                                      seqinfo,
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
        message("[.contructLociDTFromBismarkFiles] Extracting loci from ",
                "'", files[1L], "'")
    }
    loci_dt <- .readBismarkAsLociDT(
        file = files[[1L]],
        rmZeroCov = rmZeroCov,
        strandCollapse = strandCollapse,
        verbose = subverbose)

    # Identify loci not found in first file.
    loci_from_other_files_dt <- bplapply(files[-1L], function(file, loci_dt) {
        # TODO: This message won't appear in main process, so probably remove.
        if (verbose) {
            message("[.contructLociDTFromBismarkFiles] Extracting loci from ",
                    "'", file, "'")
        }
        # Read this file.
        loci_from_this_file_dt <- .readBismarkAsLociDT(
            file = file,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            verbose = subverbose)
        .subsetByOverlaps_lociDT(loci_from_this_file_dt, loci_dt,
                                 type = "equal", invert = TRUE)
    }, loci_dt = loci_dt, BPPARAM = BPPARAM)

    # Take union of FWGRanges.
    loci_dt <- funion(loci_dt, Reduce(funion, loci_from_other_files_dt))

    # Construct seqinfo if none supplied
    if (is.null(seqinfo)) {
        # TODO: Document that sort = TRUE is used.
        seqinfo <- .constructSeqinfoFromLociDT(loci_dt, sort = TRUE)
    }
    # Sort the loci using the "natural order" of ordering the elements of a
    # GenomicRanges object (see ?`sort,GenomicRanges-method`)
    loci_dt[, seqnames := factor(seqnames, levels = seqlevels(seqinfo))]
    setkey(loci_dt, seqnames, strand, start)
    list(loci_dt = loci_dt, seqinfo = seqinfo)
}


.constructCountsFromBismarkFileAndLociDT <- function(b, files, strandCollapse,
                                                     loci_dt, seqinfo, grid,
                                                     M_sink, Cov_sink,
                                                     M_sink_lock, Cov_sink_lock,
                                                     verbose, BPPARAM) {
    # Read in data for b-th file -----------------------------------------------
    file <- files[b]
    if (verbose) {
        message("[.constructCounts] Extracting counts from ", "'", file,
                "'")
    }
    bsseq_dt <- .readBismarkAsBSseqDT(
        file = file,
        # NOTE: Don't remove loci with zero coverage (it's unnecessary).
        rmZeroCov = FALSE,
        strandCollapse = strandCollapse,
        check = TRUE,
        verbose = verbose)
    # Sort bsseq_dt to use same order as loci_dt.
    bsseq_dt[, seqnames := factor(seqnames, seqlevels(seqinfo))]
    setkeyv(bsseq_dt, key(loci_dt))

    # Combine data for b-th file with loci_dt and construct counts -------------
    # TODO: (UP TO HERE) Need to use GenomicRanges' strand matching behaviour
    # NOTE: Can't use use nomatch = 0, so have to do this in two steps.
    counts <- bsseq_dt[loci_dt, c("M", "U"), nomatch = NA]
    counts[is.na(M), M := 0L]
    counts[is.na(U), U := 0L]
    counts[, c("Cov", "U") := .(M + U, NULL)]

    # Return counts or write them to the RealizationSink objects ---------------

    if (is.null(M_sink)) {
        return(counts)
    }
    # Write to M_sink and Cov_sink while respecting the IPC locks.
    ipclock(M_sink_lock)
    write_block_to_sink(as.matrix(counts[["M"]]), M_sink, grid[[b]])
    ipcunlock(M_sink_lock)
    ipclock(Cov_sink_lock)
    write_block_to_sink(as.matrix(counts[["Cov"]]), Cov_sink, grid[[b]])
    ipcunlock(Cov_sink_lock)
    NULL
}

.constructCountsFromBismarkFilesAndLociDT <- function(files,
                                                      loci_dt,
                                                      seqinfo = seqinfo,
                                                      strandCollapse,
                                                      verbose,
                                                      BPPARAM,
                                                      BACKEND) {
    # Set up ArrayGrid so that each block contains data for a single sample.
    ans_nrow <- nrow(loci_dt)
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
        FUN = .constructCountsFromBismarkFileAndLociDT,
        files = files,
        strandCollapse = strandCollapse,
        loci_dt = loci_dt,
        seqinfo = seqinfo,
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
