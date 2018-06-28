# Internal functions -----------------------------------------------------------

# TODO: Test against some of Bismark's other output files. May need a stricter
#       and more accurate heuristic to avoid 'passing' bad files.
# TODO: Test for 'bismark_methylation_extractor' and 'bedGraph' formats
#       (https://github.com/FelixKrueger/Bismark/tree/master/Docs#output-1)
.guessBismarkFileType <- function(files, n_max = 10L) {
    guessed_file_types <- setNames(vector("character", length(files)), files)
    for (file in files) {
        # NOTE: Not using readr::count_fields() because it is very slow on
        #       large, compressed files. This is because it reads the entire
        #       file into memory as a raw vector regardless of the value of
        #       `n_max` (see https://github.com/tidyverse/readr/issues/610).
        #       But good old utils::read.delim() doesn't have this limitation!
        x <- read.delim(
            file = file,
            header = FALSE,
            nrows = n_max,
            stringsAsFactors = FALSE)
        n_fields <- ncol(x)
        if (isTRUE(all(n_fields == 6L))) {
            guessed_file_types[file] <- "cov"
        } else if (isTRUE(all(n_fields == 7L))) {
            guessed_file_types[file] <- "cytosineReport"
        } else {
            stop("Could not guess Bismark file type for '", file, "'")
        }
    }
    guessed_file_types
}

# TODO: (long term) Choose between utils::read.delim(), readr::read_tsv(), and
#       data.table::fread() based on 'file'. If plain text, use fread(). If a
#       compressed file, use readr::read_tsv() if available, otherwise
#       utils::read_delim(). Longer term, combine data.table::fread() with
#       shell commands (where available) to pass compressed files. Will need to
#       be careful of the interaction between BPPARAM and fread()'s nThread.
#       Once implemented, move readr to Suggests. Finally, allow user to
#       specify which function to use.
# NOTE: In brief  benchmarking, readr::read_csv() is ~1.3-1.6x faster than
#       utils::read.delim() when reading a gzipped file, albeit it with ~1.6-2x
#       more total memory allocated. Therefore, there may be times users prefer
#       to trade off faster speed for lower memory usage.
# TODO: (long term) Formalise these benchmarks as a document in the bsseq
#       package so that we can readily re-visit these as needed.
.readBismarkAsDT <- function(file,
                             col_spec = c("all", "BSseq", "GRanges"),
                             check = FALSE,
                             verbose = FALSE,
                             ...) {
    # Quieten R CMD check about 'no visible binding for global variable' -------
    M <- U <- NULL

    # Construct the column spec ------------------------------------------------

    col_spec <- match.arg(col_spec)
    file_type <- .guessBismarkFileType(file)
    # TODO: Test for 'bismark_methylation_extractor' and 'bedGraph' formats,
    #       and error out (they're not supported).
    stopifnot(isTRUEorFALSE(check))
    if (file_type == "cov") {
        col_names <- c("seqnames", "start", "end", "beta", "M", "U")
        if (col_spec == "BSseq") {
            cols <- cols_only(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                end = col_skip(),
                beta = col_skip(),
                M = col_integer(),
                U = col_integer())
        } else if (col_spec == "GRanges") {
            cols <- cols(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                end = col_skip(),
                beta = col_skip(),
                M = col_skip(),
                U = col_skip())
        } else if (col_spec == "all") {
            cols <- cols(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                end = col_integer(),
                beta = col_double(),
                M = col_integer(),
                U = col_integer())
        }
    } else if (file_type == "cytosineReport") {
        col_names = c("seqnames", "start", "strand", "M", "U",
                      "dinucleotide_context", "trinucleotide_context")
        if (col_spec == "BSseq") {
            cols <- cols_only(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "GRanges") {
            cols <- cols_only(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_skip(),
                U = col_skip(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "all") {
            cols <- cols_only(
                seqnames = col_factor(levels = NULL),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_character(),
                trinucleotide_context = col_character())
        }
    }

    # Read the file ------------------------------------------------------------

    if (verbose) {
        message("[.readBismarkAsDT] Reading file '", file, "'")
    }
    ptime1 <- proc.time()
    x <- read_tsv(
        file = file,
        col_names = col_names,
        col_types = cols,
        na = character(),
        quoted_na = FALSE,
        progress = verbose,
        ...)

    # Construct the result -----------------------------------------------------

    x <- setDT(x)
    if (check && all(c("M", "U") %in% colnames(x))) {
        if (verbose) {
            message("[.readBismarkAsDT] Checking validity of counts in file.")
        }
        # NOTE: .checkMandCov() only accepts matrix input
        # M <- as.matrix(x[["M"]])
        # Cov <- as.matrix(x[["M"]] + x[["U"]])
        # msg <- .checkMandCov(M, Cov)
        # if (!is.null(msg)) {
        #     stop(msg)
        # } else {
        #     if (verbose) {
        #         message("[.readBismarkAsDT] Counts in file are valid!")
        #     }
        # }
        # TODO: Write .checkMandU()? Should be written to avoid copying, use a
        #       low amount of memory, stop as soon as an error is detected,
        #       and be check the exact same stuff as .checkMandCov().
        # TODO: Benchmark the below implementation of .checkMandU().
        valid <- x[, isTRUE(all(M >= 0L & U >= 0L))]
        if (!valid) {
            stop("[.readBismarkAsDT] Invalid counts detected!\n",
                 "M and U columns should be non-negative integers.")
        } else {
            if (verbose) {
                message("[.readBismarkAsDT] All counts in file are valid.")
            }
        }
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    x
}


.constructCountsFromSingleFile <- function(b, files, strandCollapse, loci,
                                           grid, M_sink, Cov_sink, sink_lock,
                                           verbose, BPPARAM) {

    # Quieten R CMD check about 'no visible binding for global variable' -------
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
        verbose = verbose)
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
    mcols <- S4Vectors:::make_zero_col_DataFrame(length(ranges))
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
    ipclock(sink_lock)
    write_block(x = M_sink, viewport = grid[[b]], block = M)
    write_block(x = Cov_sink, viewport = grid[[b]], block = Cov)
    ipcunlock(sink_lock)
    NULL
}

.constructCounts <- function(files, loci, strandCollapse, verbose, BPPARAM,
                             BACKEND, ...) {
    # Set up ArrayGrid so that each block contains data for a single sample.
    ans_nrow <- length(loci)
    ans_ncol <- length(files)
    grid <- RegularArrayGrid(c(ans_nrow, ans_ncol), c(ans_nrow, 1L))
    # Construct RealizationSink objects.
    if (is.null(BACKEND)) {
        M_sink <- NULL
        Cov_sink <- NULL
        sink_lock <- NULL
    } else if (BACKEND == "HDF5Array") {
        M_sink <- HDF5RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            # NOTE: Never allow dimnames.
            dimnames = NULL,
            type = "integer",
            name = "M",
            # TODO: Can chunkdim be specified if data are written to
            #       column-by-column?
            # chunkdim = NULL,
            # level = NULL,
            ...)
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- HDF5Array::HDF5RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            # NOTE: Never allow dimnames.
            dimnames = NULL,
            type = "integer",
            name = "Cov",
            # TODO: Can chunkdim be specified if data are written to
            #       column-by-column?
            # chunkdim = NULL,
            # level = NULL,
            ...)
        on.exit(close(Cov_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
    } else {
        # TODO: This branch should probably never be entered because we
        #       (implicitly) only support in-memory or HDF5Array backends.
        #       However, we retain it for now (e.g., fstArray backend would use
        #       this until a dedicated branch was implemented).
        M_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(Cov_sink), add = TRUE)
        sink_lock <- ipcid()
        on.exit(ipcremove(sink_lock), add = TRUE)
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
    #       the cost of transfering loci to the workers may be substantial).
    if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
        bptasks(BPPARAM) <- length(grid)
    }
    counts <- bptry(bplapply(
        X = seq_along(grid),
        FUN = .constructCountsFromSingleFile,
        files = files,
        strandCollapse = strandCollapse,
        loci = loci,
        grid = grid,
        M_sink = M_sink,
        Cov_sink = Cov_sink,
        sink_lock = sink_lock,
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

# Exported functions -----------------------------------------------------------

# TODO: Support passing a colData so that metadata is automatically added to
#       samples?
# TODO: Document that `...` are used to pass filepath, chunkdim, level, etc. to
#       HDF5RealizationSink().
# TODO: (long term) Formalise `...` by something analogous to the
#       BiocParallelParam class (RealizationSinkParam) i.e. something that
#       encapsulates the various arguments available when constructing a
#       RealizationSink.
# TODO: sampleNames = basename(files) isn't guaranteed to be unique.
# TODO: Document that you can't use strandCollapse with files that lack strand
#       (e.g., .cov files).
# TODO: (long term) Report a run-time warning/message if strandCollapse = TRUE
#       is used in conjunction with files without strand information.
read.bismark <- function(files,
                         sampleNames = basename(files),
                         loci = NULL,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE,
                         verbose = TRUE,
                         BPPARAM = bpparam(),
                         BACKEND = getRealizationBackend(),
                         ...,
                         fileType = c("cov", "oldBedGraph", "cytosineReport"),
                         mc.cores = 1) {
    # Argument checks ----------------------------------------------------------

    # Check for deprecated arguments and issue warning(s) if found.
    if (!missing(fileType)) {
        warning(
            "'filetype' is deprecated.\n",
            "Replaced with automatic detection of file types.\n",
            "See help(\"read.bismark\") for details.",
            call. = FALSE,
            immediate. = TRUE)
    }
    if (!missing(mc.cores)) {
        # TODO: What if user has provided a BPPARAM?
        warning(
            "'mc.cores' is deprecated.\n",
            "Replaced with 'BPPARAM = MulticoreParam(workers = mc.cores)'",
            ".\nSee help(\"BSmooth\").",
            call. = FALSE,
            immediate. = TRUE)
        BPPARAM <- MulticoreParam(workers = mc.cores)
    }
    # Check 'files' is valid.
    # TODO: Allow duplicate files? Useful in testing/debugging but generally a
    #       bad idea.
    if (anyDuplicated(files)) stop("'files' cannot have duplicate entires")
    file_exists <- file.exists(files)
    if (!isTRUE(all(file_exists))) {
        stop("These files cannot be found:\n",
             "  ", paste(files[!file_exists], collapse = "\n  "))
    }
    # Check 'sampleNames' is valid.
    # TODO: Check if this is still required
    # NOTE: SummarizedExperiment validator makes calls to identical() that will
    #       fail due to how sampleNames are propagated sometimes with and
    #       without names. To make life simpler, we simply strip the names.
    sampleNames <- unname(sampleNames)
    if (length(sampleNames) != length(files)) {
        stop("'sampleNames' must have the same length as 'files'.")
    }
    if (anyDuplicated(sampleNames)) {
        stop("'sampleNames' cannot have duplicate entires")
    }
    # Check 'rmZeroCov' and 'strandCollapse' are valid.
    stopifnot(isTRUEorFALSE(rmZeroCov))
    stopifnot(isTRUEorFALSE(strandCollapse))
    # Check 'loci' is valid.
    if (!is.null(loci)) {
        if (!is(loci, "GenomicRanges")) {
            stop("'loci' must be a GenomicRanges instance if not NULL.")
        }
        if (any(width(loci) != 1L)) {
            stop("All elements of 'loci' must have width equal to 1.")
        }
    }
    # Set verbosity used by internal functions
    subverbose <- as.logical(max(verbose - 1L, 0L))
    # Register 'BACKEND' and return to current value on exit
    current_BACKEND <- getRealizationBackend()
    on.exit(setRealizationBackend(current_BACKEND), add = TRUE)
    setRealizationBackend(BACKEND)
    # Check compatability of 'BPPARAM' with 'BACKEND'.
    if (!.areBackendsInMemory(BACKEND)) {
        if (!.isSingleMachineBackend(BPPARAM)) {
            stop("The parallelisation strategy must use a single machine ",
                 "when using an on-disk realization backend.\n",
                 "See help(\"BSmooth\") for details.",
                 call. = FALSE)
        }
    } else {
        if (!is.null(BACKEND)) {
            # NOTE: Currently do not support any in-memory realization
            #       backends. If the realization backend is NULL then an
            #       ordinary matrix is returned rather than a matrix-backed
            #       DelayedMatrix.
            stop("The '", BACKEND, "' realization backend is not supported.\n",
                 "See help(\"BSmooth\") for details.",
                 call. = FALSE)
        }
    }

    # Construct FWGRanges with all valid loci ----------------------------------

    # NOTE: "Valid loci" are those that remain after collapsing by strand (if
    #       strandCollapse == TRUE) and then removing loci with zero coverage
    #       in all samples (if rmZeroCov == TRUE).
    if (is.null(loci)) {
        ptime1 <- proc.time()
        if (verbose) {
            message("[read.bismark] Parsing files to construct valid loci ...")
        }
        loci <- .contructFWGRangesFromBismarkFiles(
            files = files,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            verbose = subverbose,
            BPPARAM = BPPARAM)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))
        }
    } else {
        ptime1 <- proc.time()
        if (verbose) message("[read.bismark] Using 'loci' as candidate loci")
        if (strandCollapse) {
            if (verbose) {
                message("[read.bismark] Collapsing strand of 'loci' ...")
            }
            ptime1 <- proc.time()
            loci <- .strandCollapse(loci)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("done in %.1f secs\n", stime))
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
                BPPARAM = BPPARAM)
            # Retain those elements of 'loci' with non-zero coverage.
            loci <- subsetByOverlaps(loci, loci_from_files, type = "equal")
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("done in %.1f secs\n", stime))
            }
        }
    }

    # Construct 'M' and 'Cov' matrices -----------------------------------------
    ptime1 <- proc.time()
    if (verbose) {
        message("[read.bismark] Parsing files to construct 'M' and 'Cov' ",
                "matrices ...")
    }
    counts <- .constructCounts(
        files = files,
        loci = loci,
        strandCollapse = strandCollapse,
        verbose = subverbose,
        BPPARAM = BPPARAM,
        BACKEND = BACKEND,
        ...)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }

    # Construct BSseq object ---------------------------------------------------

    if (verbose) {
        cat(sprintf("[read.bismark] Constructing BSseq object ... "))
    }
    se <- SummarizedExperiment(
        assays = counts,
        # NOTE: For now, have to use GRanges instance as rowRanges but
        #       ultimately would like to use a FWGRanges instance.
        rowRanges = as(loci, "GRanges"),
        colData = DataFrame(row.names = sampleNames))
    # TODO: Is there a way to use the internal constructor with `check = FALSE`?
    #       Don't need to check M and Cov because this has already happened
    #       when files were parsed.
    # .BSseq(se, trans = function(x) NULL, parameters = list())
    new2("BSseq", se, check = FALSE)
}

# TODOs ------------------------------------------------------------------------

# TODO: Consolidate use of message()/cat()/etc.
# TODO: Add function like minfi::read.metharray.sheet()?
# TODO: Should BACKEND really be an argument of read.bismark(); see related
#       issue on minfi repo https://github.com/hansenlab/minfi/issues/140
# TODO: May receive warning "In read_tokens_(data, tokenizer, col_specs, col_names,  ... : length of NULL cannot be changed". This is fixed in devel version of
#       readr (https://github.com/tidyverse/readr/issues/833)
# TODO: Think about naming scheme for functions. Try to have the function that
#       is bpapply()-ed have a similar name to its parent.
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
# TODO: Use verbose = getOption("verbose") as default.
