# Internal functions -----------------------------------------------------------

# TODO: Test against some of Bismark's other output files. May need a stricter
#       and more accurate heuristic to avoid 'passing' bad files.
# TODO: Test for 'bismark_methylation_extractor' and 'bedGraph' formats
#       (https://github.com/FelixKrueger/Bismark/tree/master/Docs#output-1)
#       and error out (they're not supported by .readBismaskAsDT()).
.guessBismarkFileType <- function(files, n_max = 10L) {
    vapply(files, function(file) {
        x <- read.delim(
            file = file,
            header = FALSE,
            nrows = n_max,
            stringsAsFactors = FALSE)
        n_fields <- ncol(x)
        if (isTRUE(all(n_fields == 6L))) {
            return("cov")
        }
        if (isTRUE(all(n_fields == 7L))) {
            return("cytosineReport")
        }
        # Couldn't guess the file type.
        NA_character_
    }, character(1L))
}

# NOTE: In brief benchmarking, readr::read_csv() is ~1.3-1.6x faster than
#       utils::read.delim() when reading a gzipped file, albeit it with ~1.6-2x
#       more total memory allocated. Therefore, there may be times users prefer
#       to trade off faster speed for lower memory usage.
# TODO: (long term) Formalise benchmarks of utils::read.table(),
#       data.table::fread(), and readr::read_tsv() as a document in the bsseq
#       package so that we can readily re-visit these as needed.
.readBismarkAsDT <- function(file,
                             col_spec = c("all", "BSseq", "GRanges"),
                             check = FALSE,
                             showProgress = FALSE,
                             nThread = 1L,
                             verbose = FALSE) {
    # Argument checks ----------------------------------------------------------

    file <- file_path_as_absolute(file)
    file_type <- .guessBismarkFileType(file)
    col_spec <- match.arg(col_spec)
    stopifnot(isTRUEorFALSE(check))
    stopifnot(isTRUEorFALSE(showProgress))
    # NOTE: 'verbose' controls the verbosity of .readBismarkAsDT() and
    #       'fread_verbose' is derived from that to control the verbosity of
    #       data.table::fread().
    fread_verbose <- as.logical(max(verbose - 1L, 0L))

    # Construct the column spec ------------------------------------------------

    if (file_type == "cov") {
        header <- FALSE
        if (col_spec == "BSseq") {
            colClasses <- c("factor", "integer", "NULL", "NULL", "integer",
                            "integer")
            drop <- c(3L, 4L)
            col.names <- c("seqnames", "start", "M", "U")
        } else if (col_spec == "GRanges") {
            colClasses <- c("factor", "integer", "NULL", "NULL", "NULL", "NULL")
            drop <- c(3L, 4L, 5L, 6L)
            col.names <- c("seqnames", "start")
        } else if (col_spec == "all") {
            colClasses <- c("factor", "integer", "integer", "numeric",
                            "integer", "integer")
            drop <- integer(0L)
            col.names <- c("seqnames", "start", "end", "beta", "M", "U")
        }
    } else if (file_type == "cytosineReport") {
        header <- FALSE
        if (col_spec == "BSseq") {
            colClasses <- c("factor", "integer", "factor", "integer",
                            "integer", "NULL", "NULL")
            drop <- c(6L, 7L)
            col.names <- c("seqnames", "start", "strand", "M", "U")
        } else if (col_spec == "GRanges") {
            colClasses <- c("factor", "integer", "factor", "NULL", "NULL",
                            "NULL", "NULL")
            drop <- c(4L, 5L, 6L, 7L)
            col.names <- c("seqnames", "start", "strand")
        } else if (col_spec == "all") {
            colClasses <- c("factor", "integer", "factor", "integer",
                            "integer", "factor", "factor")
            drop <- integer(0L)
            col.names <- c("seqnames", "start", "strand", "M", "U",
                           "dinucleotide_context", "trinucleotide_context")
        }
    }

    # Read the file ------------------------------------------------------------

    if (verbose) {
        message("[.readBismarkAsDT] Parsing '", file, "'")
    }
    ptime1 <- proc.time()
    if (verbose) {
        message("[.readBismarkAsDT] Reading file ...")
    }
    x <- fread(
        # TODO: This should be file = file, but need to wait for data.table
        #       v.1.11.9 for fix.
        input = file,
        sep = "\t",
        header = header,
        verbose = fread_verbose,
        drop = drop,
        colClasses = colClasses,
        col.names = col.names,
        quote = "",
        showProgress = showProgress,
        nThread = nThread)

    # Construct the result -----------------------------------------------------

    if (!is.null(x[["strand"]])) {
        x[, strand := strand(strand)]
    }
    # NOTE: Quieten R CMD check about 'no visible binding for global variable'.
    M <- U <- NULL
    if (check && all(c("M", "U") %in% colnames(x))) {
        if (verbose) {
            message("[.readBismarkAsDT] Checking validity of counts in file.")
        }
        valid <- x[, isTRUE(all(M >= 0L & U >= 0L))]
        if (!valid) {
            stop("[.readBismarkAsDT] Invalid counts detected.\n",
                 "'M' and 'U' columns should be non-negative integers.")
        } else {
            if (verbose) {
                message("[.readBismarkAsDT] All counts in file are valid.")
            }
        }
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        message("Done in ", round(stime, 1), " secs")
    }
    x
}

.constructCountsFromSingleFile <- function(b, files, loci, strandCollapse,
                                           grid, M_sink, Cov_sink, sink_lock,
                                           nThread, verbose) {

    # Argument checks ----------------------------------------------------------

    subverbose <- max(as.integer(verbose) - 1L, 0L)
    # NOTE: Quieten R CMD check about 'no visible binding for global variable'.
    M <- U <- NULL

    # Read b-th file to construct data.table of valid loci and their counts ----

    file <- files[b]
    if (verbose) {
        message("[.constructCountsFromBismarkFileAndFWGRanges] Extracting ",
                "counts from ", "'", file, "'")
    }
    dt <- .readBismarkAsDT(
        file = file,
        col_spec = "BSseq",
        check = TRUE,
        nThread = nThread,
        verbose = subverbose)
    # NOTE: Can only collapse by strand if the data are stranded!
    if (strandCollapse && !is.null(dt[["strand"]])) {
        # Shift loci on negative strand by 1 to the left and then remove
        # strand since no longer valid.
        dt[strand == "-", start := start - 1L][, strand := NULL]
        # Aggregate counts at loci with the same 'seqnames' and 'start'.
        dt <- dt[, list(M = sum(M), U = sum(U)), by = c("seqnames", "start")]
    }

    # Construct FWGRanges ------------------------------------------------------

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
    loci_from_this_sample <- .FWGRanges(
        seqnames = seqnames,
        ranges = ranges,
        strand = strand,
        seqinfo = seqinfo,
        elementMetadata = mcols)

    # Construct 'M' and 'Cov' matrices -----------------------------------------

    ol <- findOverlaps(loci_from_this_sample, loci, type = "equal")
    M <- matrix(rep(0L, length(loci)), ncol = 1)
    Cov <- matrix(rep(0L, length(loci)), ncol = 1)
    M[subjectHits(ol)] <- dt[queryHits(ol), ][["M"]]
    Cov[subjectHits(ol)] <- dt[queryHits(ol), list(Cov = (M + U))][["Cov"]]

    # Return 'M' and 'Cov' or write them to the RealizationSink objects --------

    if (is.null(M_sink)) {
        return(list(M = M, Cov = Cov))
    }
    # Write to M_sink and Cov_sink while respecting the IPC lock.
    viewport <- grid[[b]]
    ipclock(sink_lock)
    write_block(M_sink, viewport = viewport, block = M)
    write_block(Cov_sink, viewport = viewport, block = Cov)
    ipcunlock(sink_lock)
    NULL
}

.constructCounts <- function(files, loci, strandCollapse, BPPARAM,
                             BACKEND, dir, chunkdim, level, nThread,
                             verbose = FALSE) {

    # Argument checks ----------------------------------------------------------

    subverbose <- max(as.integer(verbose) - 1L, 0L)

    # Construct ArrayGrid ------------------------------------------------------

    # NOTE: Each block contains data for a single sample.
    ans_nrow <- length(loci)
    ans_ncol <- length(files)
    ans_dim <- c(ans_nrow, ans_ncol)
    grid <- RegularArrayGrid(
        refdim = ans_dim,
        spacings = c(ans_nrow, 1L))

    # Construct RealizationSink objects ----------------------------------------

    if (is.null(BACKEND)) {
        M_sink <- NULL
        Cov_sink <- NULL
        sink_lock <- NULL
    } else if (BACKEND == "HDF5Array") {
        h5_path <- file.path(dir, "assays.h5")
        M_sink <- HDF5RealizationSink(
            dim = ans_dim,
            dimnames = NULL,
            type = "integer",
            filepath = h5_path,
            name = "M",
            chunkdim = chunkdim,
            level = level)
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- HDF5RealizationSink(
            dim = ans_dim,
            dimnames = NULL,
            type = "integer",
            filepath = h5_path,
            name = "Cov",
            chunkdim = chunkdim,
            level = level)
        on.exit(close(Cov_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
    } else {
        # TODO: This branch should probably never be entered because we
        #       (implicitly) only support in-memory or HDF5Array backends.
        #       However, we retain it for now (e.g., fstArray backend would use
        #       this until a dedicated branch was implemented).
        M_sink <- DelayedArray::AutoRealizationSink(
            dim = ans_dim,
            type = "integer")
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- DelayedArray::AutoRealizationSink(
            dim = ans_dim,
            type = "integer")
        on.exit(close(Cov_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
    }

    # Apply .constructCountsFromSingleFile() to each file ----------------------

    # Set number of tasks to ensure the progress bar gives frequent updates.
    # NOTE: The progress bar increments once per task
    #       (https://github.com/Bioconductor/BiocParallel/issues/54).
    #       Although it is somewhat of a bad idea to overrides a user-specified
    #       bptasks(BPPARAM), the value of bptasks(BPPARAM) doesn't affect
    #       performance in this instance, and so we opt for a useful progress
    #       bar. Only SnowParam (and MulticoreParam by inheritance) have a
    #       bptasks<-() method.
    if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
        bptasks(BPPARAM) <- length(grid)
    }
    counts <- bptry(bplapply(
        X = seq_along(grid),
        FUN = .constructCountsFromSingleFile,
        files = files,
        loci = loci,
        strandCollapse = strandCollapse,
        grid = grid,
        M_sink = M_sink,
        Cov_sink = Cov_sink,
        sink_lock = sink_lock,
        verbose = subverbose,
        nThread = nThread,
        BPPARAM = BPPARAM))
    if (!all(bpok(counts))) {
        stop(".constructCounts() encountered errors for these files:\n  ",
             paste(files[!bpok], collapse = "\n  "))
    }

    # Construct 'M' and 'Cov' --------------------------------------------------

    if (is.null(BACKEND)) {
        # Returning matrix objects.
        M <- do.call(c, lapply(counts, "[[", "M"))
        attr(M, "dim") <- ans_dim
        Cov <- do.call(c, lapply(counts, "[[", "Cov"))
        attr(Cov, "dim") <- ans_dim
    } else {
        # Returning DelayedMatrix objects.
        M <- as(M_sink, "DelayedArray")
        Cov <- as(Cov_sink, "DelayedArray")
    }

    return(list(M = M, Cov = Cov))
}

# Exported functions -----------------------------------------------------------

# TODO: If you have N cores available, are you better off using
#       bpworkers() = N in the BPPARAM or nThread = N and use
#       data.table::fread()? Or something in between?
# TODO: Allow user to determine `nThread`?
# TODO: Document that you can't use strandCollapse with files that lack strand
#       (e.g., .cov files).
# TODO: (long term) Report a run-time warning/message if strandCollapse = TRUE
#       is used in conjunction with files without strand information.
# TODO: Try a bpiterate()-based approach to getting the set of loci from files.
read.bismark <- function(files,
                         loci = NULL,
                         colData = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE,
                         BPPARAM = bpparam(),
                         BACKEND = NULL,
                         dir = tempfile("BSseq"),
                         replace = FALSE,
                         chunkdim = NULL,
                         level = NULL,
                         nThread = 1L,
                         verbose = getOption("verbose")) {
    # Argument checks ----------------------------------------------------------

    # Check 'files' is valid.
    if (anyDuplicated(files)) {
        stop("'files' cannot have duplicate entries.")
    }
    file_exists <- file.exists(files)
    if (!isTRUE(all(file_exists))) {
        stop("These files cannot be found:\n  ",
             paste(files[!file_exists], collapse = "\n  "))
    }
    guessed_file_types <- .guessBismarkFileType(files)
    if (anyNA(guessed_file_types)) {
        stop("Could not guess Bismark file type for these files:\n  ",
             "  ", paste(files[!is.na(guessed_file_types)], collapse = "\n  "))
    }
    # Check 'loci' is valid.
    if (!is.null(loci)) {
        if (!is(loci, "GenomicRanges")) {
            stop("'loci' must be a GenomicRanges instance if not NULL.")
        }
        if (any(width(loci) != 1L)) {
            stop("All elements of 'loci' must have width equal to 1.")
        }
    }
    # Check 'colData' is valid.
    if (is.null(colData)) {
        colData <- DataFrame(row.names = files)
    }
    if (nrow(colData) != length(files)) {
        stop("Supplied 'colData' must have nrow(colData) == length(files).")
    }
    # Check 'rmZeroCov' and 'strandCollapse' are valid.
    stopifnot(isTRUEorFALSE(rmZeroCov))
    stopifnot(isTRUEorFALSE(strandCollapse))
    # Register 'BACKEND' and return to current value on exit.
    # TODO: Is this strictly necessary?
    current_BACKEND <- getAutoRealizationBackend()
    on.exit(setAutoRealizationBackend(current_BACKEND), add = TRUE)
    setAutoRealizationBackend(BACKEND)
    # Check compatability of 'BPPARAM' with 'BACKEND'.
    if (!.areBackendsInMemory(BACKEND)) {
        if (!.isSingleMachineBackend(BPPARAM)) {
            stop("The parallelisation strategy must use a single machine ",
                 "when using an on-disk realization backend.\n",
                 "See help(\"read.bismark\") for details.",
                 call. = FALSE)
        }
    } else {
        if (!is.null(BACKEND)) {
            # NOTE: Currently do not support any in-memory realization
            #       backends. If the realization backend is NULL then an
            #       ordinary matrix is returned rather than a matrix-backed
            #       DelayedMatrix.
            stop("The '", BACKEND, "' realization backend is not supported.",
                 "\n  See help(\"read.bismark\") for details.",
                 call. = FALSE)
        }
    }
    # If using HDF5Array as BACKEND, check remaining options are sensible.
    if (identical(BACKEND, "HDF5Array")) {
        # NOTE: Most of this copied from
        #       HDF5Array::saveHDF5SummarizedExperiment().
        if (!isSingleString(dir)) {
            stop(wmsg("'dir' must be a single string specifying the path to ",
                      "the directory where to save the BSseq object (the ",
                      "directory will be created)."))
        }
        if (!isTRUEorFALSE(replace)) {
            stop("'replace' must be TRUE or FALSE")
        }
        if (!dir.exists(dir)) {
            HDF5Array::create_dir(dir)
        } else {
            HDF5Array::replace_dir(dir, replace)
        }
    }
    # Set verbosity used by internal functions.
    subverbose <- as.logical(max(verbose - 1L, 0L))

    # Construct FWGRanges with all valid loci ----------------------------------

    # NOTE: "Valid" loci are those that remain after collapsing by strand (if
    #       strandCollapse == TRUE) and then removing loci with zero coverage
    #       in all samples (if rmZeroCov == TRUE).
    if (is.null(loci)) {
        ptime1 <- proc.time()
        if (verbose) {
            message(
                "[read.bismark] Parsing files and constructing valid loci ...")
        }
        loci <- .contructFWGRangesFromBismarkFiles(
            files = files,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            verbose = subverbose,
            nThread = nThread,
            BPPARAM = BPPARAM)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            message("Done in ", round(stime, 1), " secs")
        }
    } else {
        if (verbose) {
            message("[read.bismark] Using 'loci' as candidate loci.")
        }
        if (strandCollapse) {
            if (verbose) {
                message("[read.bismark] Collapsing strand of 'loci' ...")
            }
            ptime1 <- proc.time()
            loci <- .strandCollapse(loci)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                message("Done in ", round(stime, 1), " secs")
            }
        }
        if (rmZeroCov) {
            if (verbose) {
                message("[read.bismark] Parsing files to identify elements of ",
                        "'loci' with non-zero coverage ...")
            }
            ptime1 <- proc.time()
            # Construct loci with non-zero coverage in files.
            loci_from_files <- .contructFWGRangesFromBismarkFiles(
                files = files,
                rmZeroCov = rmZeroCov,
                strandCollapse = strandCollapse,
                verbose = subverbose,
                nThread = nThread,
                BPPARAM = BPPARAM)
            # Retain those elements of 'loci' with non-zero coverage.
            loci <- subsetByOverlaps(loci, loci_from_files, type = "equal")
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                message("Done in ", round(stime, 1), " secs")
            }
        }
    }

    # Construct 'M' and 'Cov' matrices -----------------------------------------
    ptime1 <- proc.time()
    if (verbose) {
        message("[read.bismark] Parsing files and constructing 'M' and 'Cov' ",
                "matrices ...")
    }
    counts <- .constructCounts(
        files = files,
        loci = loci,
        strandCollapse = strandCollapse,
        BPPARAM = BPPARAM,
        BACKEND = BACKEND,
        dir = dir,
        chunkdim = chunkdim,
        level = level,
        nThread = nThread,
        verbose = subverbose)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        message("Done in ", round(stime, 1), " secs")
    }

    # Construct BSseq object, saving it if it is HDF5-backed -------------------

    if (verbose) {
        message("[read.bismark] Constructing BSseq object ... ")
    }
    se <- SummarizedExperiment(
        assays = counts,
        # NOTE: For now, have to use GRanges instance as rowRanges but
        #       ultimately would like to use a FWGRanges instance.
        rowRanges = as(loci, "GRanges"),
        colData = colData)
    # TODO: Is there a way to use the internal constructor with `check = FALSE`?
    #       Don't need to check M and Cov because this has already happened
    #       when files were parsed.
    # .BSseq(se, trans = function(x) NULL, parameters = list())
    bsseq <- new2("BSseq", se, check = FALSE)
    if (!is.null(BACKEND) && BACKEND == "HDF5Array") {
        # NOTE: Save BSseq object; mimicing
        #       HDF5Array::saveHDF5SummarizedExperiment().
        x <- bsseq
        x@assays <- HDF5Array::shorten_assay2h5_links(x@assays)
        saveRDS(x, file = file.path(dir, "se.rds"))
    }
    bsseq
}

# TODOs ------------------------------------------------------------------------

# TODO: Document internal functions for my own sanity. Also, some may be useful
#       to users of bsseq (although I still won't export these for the time
#       being).
# TODO: (long term) Current implementation requires that user can load at least
#       one sample's worth of data into memory per worker. Could instead read
#       chunks of data, write to sink, load next chunk, etc.
# TODO: Document that if 'loci' is NULL and any 'files' (especially the first
#       file) are .cov files, then any loci present in the .cov files will have
#       their strand set to *. If you are mixing and matching .cov and
#       .cytosineReport files and don't want this behaviour (i.e. you want to
#       retain strand) then you'll need to construct your own 'gr' and pass
#       this to the function. Add unit tests for this behaviour.
