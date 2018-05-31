# Internal functions -----------------------------------------------------------

# TODO: Test against some of Bismark's other output files. May need a stricter
#       and more accurate heuristic to avoid 'passing' bad files.
# NOTE: Not using readr::count_fields() because it is very slow on large,
#       compressed files. This is because it reads the entire file into memory
#       as a raw vector regardless of the value of `n_max` (see
#       https://github.com/tidyverse/readr/issues/610). But good old
#       utils::read.delim() doesn't have this limitation!
.guessBismarkFileType <- function(files, n_max = 10L) {
    guessed_file_types <- setNames(vector("character", length(files)), files)
    for (file in files) {
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

# TODO: (longterm, see "Alternatively ..." for a better idea)
#       .readBismarkAsDT2(): exact same as .readBismarkAsDT() but
#       uses utils::read.delim() instead of readr::read_tsv(). In brief
#       benchmarking, readr::read_csv() is ~1.3-1.6x faster than
#       utils::read.delim() when reading a gzipped file, albeit it with ~1.6-2x
#       more total memory allocated. Therefore, there may be times users prefer
#       to trade off faster speed for lower memory usage. When written, move
#       readr to Suggests. Alternatively, re-write .readBismarkAsDT() using
#       data.table::fread() for uncompressed files and utils::read.delim() for
#       compressed files. This removes the dependency on readr, albeit it with
#       slightly slower reading of compressed files. Could even then use
#       data.table::fread() coupled with shell commands (where available) to
#       pass compressed files. Ultimately, we want to use data.table beyond
#       data.table::fread() whereas readr is only used for file input.
# NOTE: This returns the file as a data.table. However, to do this it uses
#       readr::read_tsv() + data.table::setDT() instead of data.table::fread()!
#       Although the latter is faster, this uses the former because Bismark
#       files are commonly compressed and readr::read_tsv() supports reading
#       directly from compressed files whereas data.table::fread() does not.
.readBismarkAsDT <- function(file,
                             col_spec = c("all", "BSseq", "GRanges"),
                             check = FALSE,
                             verbose = FALSE,
                             ...) {
    col_spec <- match.arg(col_spec)
    file_type <- .guessBismarkFileType(file)
    stopifnot(S4Vectors:::isTRUEorFALSE(check))
    if (file_type == "cov") {
        col_names <- c("seqnames", "start", "end", "beta", "M", "U")
        if (col_spec == "BSseq") {
            cols <- cols_only(
                seqnames = col_character(),
                start = col_integer(),
                end = col_skip(),
                beta = col_skip(),
                M = col_integer(),
                U = col_integer())
        } else if (col_spec == "GRanges") {
            cols <- cols(
                seqnames = col_character(),
                start = col_integer(),
                end = col_skip(),
                beta = col_skip(),
                M = col_skip(),
                U = col_skip())
        } else if (col_spec == "all") {
            cols <- cols(
                seqnames = col_character(),
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
                seqnames = col_character(),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "GRanges") {
            cols <- cols_only(
                seqnames = col_character(),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_skip(),
                U = col_skip(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "all") {
            cols <- cols_only(
                seqnames = col_character(),
                start = col_integer(),
                strand = col_factor(levels(strand())),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_character(),
                trinucleotide_context = col_character())
        }
    }
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
    # Take union of loci.
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
                                                     loci_dt, grid, M_sink,
                                                     Cov_sink, M_sink_lock,
                                                     Cov_sink_lock, verbose,
                                                     BPPARAM) {
    # Read in data for b-th file -----------------------------------------------
    file <- files[b]
    if (verbose) {
        message("[.constructCounts] Extracting counts from ", "'", file,
                "'")
    }
    bsseq_dt <- .readBismarkAsBSseqDT(
        file = file,
        # NOTE: No need to remove zero coverage loci because we combine with
        #       loci_dt, which by construction only contains loci with non-zero
        #       coverage.
        rmZeroCov = FALSE,
        strandCollapse = strandCollapse,
        check = TRUE,
        verbose = verbose)
    # Sort bsseq_dt to use same order as loci_dt.
    bsseq_dt[, c("seqnames", "strand") :=
                 .(factor(seqnames, levels(loci_dt[["seqnames"]])),
                   strand(strand))]
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

# Exported functions -----------------------------------------------------------

# TODO: Support BPREDO?
# TODO: Support passing a colData so that metadata is automatically added to
#       samples?
# TODO: Probably pass down verbose as subverbose to functions that take verbose
#       argument.
read.bismark <- function(files,
                         sampleNames = basename(files),
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE,
                         fileType = c("cov", "oldBedGraph", "cytosineReport"),
                         seqinfo = NULL,
                         gr = NULL,
                         mc.cores = 1,
                         verbose = TRUE,
                         BPPARAM = bpparam(),
                         BACKEND = getRealizationBackend()) {
    # Argument checks ----------------------------------------------------------

    # TODO: Register BACKEND and return to current value on exit!

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
    if (!missing(verbose)) {
        warning(
            "'verbose' is deprecated.\n",
            "Replaced by setting 'bpprogressbar(BPPARAM) <- TRUE'.\n",
            "See help(\"BSmooth\") for details.",
            call. = FALSE,
            immediate. = TRUE)
        if (verbose) bpprogressbar(BPPARAM) <- TRUE
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
    # Check 'gr' and 'seqinfo' are compatible
    if (!is.null(gr)) {
        if (is.null(seqinfo)) {
            # TODO: stop() if both are provided?
            warning("'seqinfo' is ignored when 'gr' is supplied.")
        }
        seqinfo <- seqinfo(gr)
    }

    # Construct data.table and GRanges with all valid loci ---------------------

    # NOTE: "Valid loci" are those that remain after collapsing by strand (if
    #       strandCollapse == TRUE) and then removing loci with zero coverage
    #       in all samples (if rmZeroCov == TRUE).
    ptime1 <- proc.time()
    if (is.null(gr)) {
        if (verbose) {
            message("[read.bismark] Reading files to construct GRanges with ",
                    "valid loci ...")
        }
        loci_dt_and_seqinfo <- .contructLociDTAndSeqinfoFromBismarkFiles(
            files = files,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            seqinfo = seqinfo,
            verbose = subverbose,
            BPPARAM = BPPARAM)
        loci_dt <- loci_dt_and_seqinfo[["loci_dt"]]
        seqinfo <- loci_dt_and_seqinfo[["seqinfo"]]
        gr <- .lociDTAsGRanges(loci_dt, seqinfo)
    } else {
        if (verbose) message("[read.bismark] Using 'gr' as GRanges with loci")
        loci_dt <- .grAsLociDT(gr)
        if (strandCollapse) {
            if (verbose) {
                message("[read.bismark] Collapsing strand of loci in 'gr' ...")
            }
            # ptime1 <- proc.time()
            loci_dt <- .strandCollapseLociDT(loci_dt, has_counts = FALSE)
            # ptime2 <- proc.time()
            # stime <- (ptime2 - ptime1)[3]
            # if (verbose) {
            #     cat(sprintf("done in %.1f secs\n", stime))
            # }
        }
        if (rmZeroCov) {
            # TODO: Test this branch
            # NOTE: Have to parse files in order to identify loci with zero
            #       coverage in all samples.
            if (verbose) {
                message("[read.bismark] Identifying loci in 'gr' with zero ",
                        "coverage in all samples ...")
            }
            # ptime1 <- proc.time()
            loci_from_files_dt <- .contructLociDTFromBismarkFiles(
                files = files,
                rmZeroCov = rmZeroCov,
                strandCollapse = strandCollapse,
                seqinfo = seqinfo,
                verbose = subverbose)
            # Retain the intersection of loci_dt and loci_from_files_dt
            # TODO: Need to use GenomicRanges' strand matching behaviour
            loci_dt <- loci_from_files_dt[loci_dt]
            # ptime2 <- proc.time()
            # stime <- (ptime2 - ptime1)[3]
            # if (verbose) {
            #     cat(sprintf("done in %.1f secs\n", stime))
            # }
        }
        if (strandCollapse || rmZeroCov) {
            # NOTE: Have to update 'gr' if loci have been collapsed by strand
            #       or loci have been filtered to remove loci with zero
            #       coverage.
            if (verbose) {
                message("[read.bismark] Filtering 'gr' to retain valid loci ",
                        "...")
            }
            # ptime1 <- proc.time()
            gr <- .lociDTAsGRanges(loci_dt, seqinfo)
            # ptime2 <- proc.time()
            # stime <- (ptime2 - ptime1)[3]
            # if (verbose) {
            #     cat(sprintf("done in %.1f secs\n", stime))
            # }
        }
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }

    # Construct 'M' and 'Cov' matrices -----------------------------------------
    ptime1 <- proc.time()
    if (verbose) {
        message("[read.bismark] Reading files to construct 'M' and 'Cov' ",
                "matrices ...")
    }
    counts <- .constructCountsFromBismarkFilesAndLociDT(
        files = files,
        loci_dt = loci_dt,
        strandCollapse = strandCollapse,
        verbose = subverbose,
        BPPARAM = BPPARAM,
        BACKEND = BACKEND)
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
        rowRanges = gr,
        colData = DataFrame(row.names = sampleNames))
    # TODO: Is there a way to use the internal constructor with `check = FALSE`?
    # .BSseq(se, trans = function(x) NULL, parameters = list())
    new2("BSseq", se, check = FALSE)
}

# TODOs ------------------------------------------------------------------------

# TODO: The documentation needs a complete overhaul and checked against Bismark
#       docs (e.g., the description of the .cov file is wrong)
# TODO: Consolidate use of message()/cat()/etc.
# TODO: Add function like minfi::read.metharray.sheet()?
# TODO: Should BACKEND really be an argument of read.bismark(); see related
#       issue on minfi repo https://github.com/hansenlab/minfi/issues/140
# TODO: May receive warning "In read_tokens_(data, tokenizer, col_specs, col_names,  ... : length of NULL cannot be changed". This is fixed in devel version of
#       readr (https://github.com/tidyverse/readr/issues/833)
# TODO: Decide whether to preserve verbose argument of several functions.
# TODO: Document that if seqinfo is supplied then only loci from those
#       seqlevels will be retained. Similarly, if gr is supplied then only loci
#       in the gr will be retained.
# TODO: Think about naming scheme for functions. Try to have the function that
#       is bpapply()-ed have a similar name to its parent.
# TODO: Document internal functions for my own sanity. Also, some may be useful
#       to users of bsseq (although I still won't export these for the time
#       being).
# TODO: Allow user to specify HDF5 file and have both M and Cov written to that
#       file.
# TODO: Helper function to obtain CpX loci from BSgenome as GRanges to pass as
#       'gr' argument in read.bismark().
# TODO: (long term) Current implementation requires that user can load at least
#       one sample's worth of data into memory per worker. Could instead read
#       chunks of data, write to sink, load next chunk, etc.
# TODO: Add big note to documentation that .cov file does not contain strand
#       information, which means strandCollapse can't be used (unless 'gr' is
#       supplied or one of the other files is a cytosineReport file).
# TODO: Document that if 'gr' is NULL and any 'files' (especially the first
#       file) are .cov files, then any loci present in the .cov files will have
#       their strand set to *. If you are mixing and matching .cov and
#       .cytosineReport files and don't want this behaviour (i.e. you want to
#       retain strand) then you'll need to construct your own 'gr' and pass
#       this to the function. Add unit tests for this behaviour.
