# Internal functions -----------------------------------------------------------

.guessBismarkFileType <- function(files, n_max = 10L) {
    guessed_file_types <- setNames(vector("character", length(files)), files)
    for (file in files) {
        n_fields <- count_fields(file, tokenizer_tsv(), n_max = n_max)
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

.readBismark <- function(file,
                         col_spec = c("BSseq", "GRanges", "all"),
                         verbose = FALSE) {
    col_spec <- match.arg(col_spec)
    file_type <- .guessBismarkFileType(file)
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
                strand = col_character(),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "GRanges") {
            cols <- cols_only(
                seqnames = col_character(),
                start = col_integer(),
                strand = col_character(),
                M = col_skip(),
                U = col_skip(),
                dinucleotide_context = col_skip(),
                trinucleotide_context = col_skip())
        } else if (col_spec == "all") {
            cols <- cols_only(
                seqnames = col_character(),
                start = col_integer(),
                strand = col_character(),
                M = col_integer(),
                U = col_integer(),
                dinucleotide_context = col_character(),
                trinucleotide_context = col_character())
        }
    }
    if (verbose) {
        message("[.readBismark] Reading file '", file, "'")
    }
    ptime1 <- proc.time()
    x <- read_tsv(
        file = file,
        col_names = col_names,
        col_types = cols,
        progress = FALSE)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    x
}

.readBismarkAsBSseqDT <- function(file, rmZeroCov, strandCollapse, verbose) {
    x <- .readBismark(file, "BSseq", verbose)
    setDT(x)
    # Data is unstranded if none is provided
    # TODO: What's the data.table way to check if column exists and if create
    #       it if it doesn't exist?
    if (is.null(x[["strand"]])) {
        x[, strand := "*"]
    }
    if (strandCollapse) {
        x <- .strandCollapseDT(x, has_counts = TRUE)
    }
    if (rmZeroCov) {
        return(x[(M + U) > 0])
    }
    setkey(x, seqnames, strand, start)
    x
}

.constructSeqinfoFromLociDT <- function(loci_dt) {
    seqnames <- loci_dt[, as.character(unique(seqnames))]
    sortSeqlevels(Seqinfo(seqnames = seqnames))
}

.contructLociDTFromBismarkFiles <- function(files,
                                            rmZeroCov,
                                            strandCollapse,
                                            seqinfo,
                                            verbose) {
    subverbose <- max(as.integer(verbose) - 1L, 0L)

    # Initalise `loci_dt` using the first file.
    if (verbose) {
        message("[.contructLociDTFromBismarkFiles] Extracting loci from ",
                "'", files[1L], "'")
    }
    loci_dt <- .readBismarkAsBSseqDT(
        file = files[1L],
        rmZeroCov = rmZeroCov,
        strandCollapse = strandCollapse,
        verbose = subverbose)
    # Drop 'M' and 'U'
    loci_dt[, c("M", "U") := .(NULL, NULL)]

    # TODO: Do in batches in parallel (n_batches = mc.cores)
    # Loop over remaining files
    for (file in files[-1L]) {
        if (verbose) {
            message("[.contructLociDTFromBismarkFiles] Extracting loci from ",
                    "'", file, "'")
        }
        # Process this file
        loci_from_this_file_dt <- .readBismarkAsBSseqDT(
            file = file,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            verbose = subverbose)
        # Drop 'M' and 'U'
        loci_from_this_file_dt[, c("M", "U") := .(NULL, NULL)]

        # Update `loci_dt`
        loci_dt <- merge(loci_dt, loci_from_this_file_dt, all = TRUE)
    }

    # Construct seqinfo if none supplied
    if (is.null(seqinfo)) {
        seqinfo <- .constructSeqinfoFromLociDT(loci_dt)
    }
    # Sort the loci using the "natural order" of ordering the elements of a
    # GenomicRanges object (see ?`sort,GenomicRanges-method`)
    loci_dt[,
            c("seqnames", "strand") :=
                .(factor(seqnames, levels = seqlevels(seqinfo)),
                  strand(strand))]
    setkey(loci_dt, seqnames, strand, start)
    loci_dt
}

.constructCountsFromBismarkFilesAndLociDT <- function(files,
                                                      loci_dt,
                                                      strandCollapse,
                                                      mc.cores,
                                                      BACKEND,
                                                      verbose = FALSE) {
    ans_nrow <- nrow(loci_dt)
    ans_ncol <- length(files)
    if (is.null(BACKEND)) {
        M <- matrix(NA_integer_, nrow = ans_nrow, ncol = ans_ncol)
        Cov <- matrix(NA_integer_, nrow = ans_nrow, ncol = ans_ncol)
    } else {
        # Set up intermediate RealizationSink objects of appropriate dimensions
        # and type. These are ultimately coerced to the output DelayedMatrix
        # objects, `M` and `U`.
        # NOTE: Need to register the supplied BACKEND while being sure to
        #       re-register any existing backend upon exiting the function.
        current_BACKEND <- getRealizationBackend()
        on.exit(setRealizationBackend(current_BACKEND))
        setRealizationBackend(BACKEND)
        # NOTE: Don't do `U_sink <- M_sink` or else these will reference the
        #       same object and clobber each other when written to!
        M_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(M_sink), add = TRUE)
        Cov_sink <- DelayedArray:::RealizationSink(
            dim = c(ans_nrow, ans_ncol),
            type = "integer")
        on.exit(close(Cov_sink), add = TRUE)
        grid <- RegularArrayGrid(dim(M_sink), c(ans_nrow, 1L))
    }

    # TODO: Load as many files as safe according to block.size (if using HDF5,
    #       otherwise load them all up)
    for (k in seq_along(files)) {
        file <- files[k]
        if (verbose) {
            message("[.constructCounts] Extracting counts from ", "'", file,
                    "'")
        }
        bsseq_dt <- .readBismarkAsBSseqDT(
            file = file,
            rmZeroCov = FALSE,
            strandCollapse = strandCollapse,
            verbose = verbose)
        # TODO: Need to use GenomicRanges' strand matching behaviour
        counts <- bsseq_dt[loci_dt, c("M", "U"), nomatch = NA]
        counts[is.na(M), M := 0L]
        counts[is.na(U), U := 0L]
        counts[, c("Cov", "U") := .(M + U, NULL)]

        if (is.null(BACKEND)) {
            M[, k] <- counts[["M"]]
            Cov[, k] <- counts[["Cov"]]
        } else {
            viewport <- grid[[k]]
            # TODO: Does as.matrix() cause a copy?
            write_block_to_sink(
                block = as.matrix(counts[["M"]]),
                sink = M_sink,
                viewport = viewport)
            write_block_to_sink(
                block = as.matrix(counts[["Cov"]]),
                sink = Cov_sink,
                viewport = viewport)
        }
    }
    if (!is.null(BACKEND)) {
        M <- as(M_sink, "DelayedArray")
        Cov <- as(Cov_sink, "DelayedArray")
    }
    list(M = M, Cov = Cov)
}

.lociDTAsGRanges <- function(loci_dt, seqinfo = NULL) {
    GRanges(
        seqnames = loci_dt[["seqnames"]],
        ranges = IRanges(loci_dt[["start"]], width = 1L),
        strand = loci_dt[["strand"]],
        seqinfo = seqinfo)
}

.grAsLociDT <- function(gr) {
    data.table(
        seqnames = as.factor(seqnames(gr)),
        start = start(gr),
        strand = as.factor(strand(gr)),
        key = c("seqnames", "strand", "start"))
}

.strandCollapseDT <- function(x, has_counts = TRUE) {
    # Shift loci on negative strand by 1 to the left
    x[strand == "-", start := start - 1L]
    # Unstrand all loci
    x[, strand := "*"]
    # Set key
    setkey(x, seqnames, strand, start)
    # TODO: Is by
    if (has_counts) return(x[, .(M = sum(M), U = sum(U)), by = key(x)])
    unique(x)
}

# Exported functions -----------------------------------------------------------

read.bismark <- function(files,
                         sampleNames = basename(files),
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE,
                         fileType = c("cov", "oldBedGraph", "cytosineReport"),
                         seqinfo = NULL,
                         gr = NULL,
                         mc.cores = 1,
                         verbose = TRUE,
                         BACKEND = NULL) {
    # Argument checking
    if (anyDuplicated(files)) stop("'files' cannot have duplicate entires")
    file_exists <- file.exists(files)
    if (!isTRUE(all(file_exists))) {
        stop("These files cannot be found:\n",
             "  ", paste(files[!file_exists], collapse = "\n  "))
    }
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
    # TODO: What's the proper way to deprecate a function argument?
    if (!missing(fileType)) {
        warning("'fileType' is deprecated and ignored.")
    }
    if (!is.null(gr) && !is.null(seqinfo)) {
        warning("'seqinfo' is ignored when 'gr' is supplied.")
    }
    subverbose <- max(as.integer(verbose) - 1L, 0L)

    # Construct a data.table and a GRanges with all "valid loci".
    # NOTE: "Valid loci" are those that remain after collapsing by strand (if
    #       strandCollapse == TRUE) and then removing loci with zero coverage
    #       in all samples (if rmZeroCov == TRUE).
    if (is.null(gr)) {
        if (verbose) {
            message("[read.bismark] Constructing GRanges with valid loci ...")
        }
        ptime1 <- proc.time()
        loci_dt <- .contructLociDTFromBismarkFiles(
            files = files,
            rmZeroCov = rmZeroCov,
            strandCollapse = strandCollapse,
            seqinfo = seqinfo,
            verbose = subverbose)
        if (is.null(seqinfo)) {
            seqinfo <- .constructSeqinfoFromLociDT(loci_dt)
        }
        gr <- .lociDTAsGRanges(loci_dt, seqinfo)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))
        }
    } else {
        if (verbose) message("[read.bismark] Using 'gr' as GRanges with loci")
        # NOTE: Don't use as.data.table() because it will modify gr by
        #       reference (https://github.com/Rdatatable/data.table/issues/2278)
        ptime1 <- proc.time()
        loci_dt <- .grAsLociDT(gr)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))
        }
        if (strandCollapse) {
            if (verbose) {
                message("[read.bismark] Collapsing strand of loci in 'gr' ...")
            }
            ptime1 <- proc.time()
            loci_dt <- .strandCollapseDT(loci_dt, has_counts = FALSE)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("done in %.1f secs\n", stime))
            }
        }
        if (rmZeroCov) {
            # NOTE: Have to parse files in order to identify loci with zero
            #       coverage in all samples.
            if (verbose) {
                message("[read.bismark] Identifying loci in 'gr' with zero ",
                        "coverage in all samples ...")
            }
            ptime1 <- proc.time()
            loci_from_files_dt <- .contructLociDTFromBismarkFiles(
                files = files,
                rmZeroCov = rmZeroCov,
                strandCollapse = strandCollapse,
                seqinfo = seqinfo,
                verbose = subverbose)
            # Retain the intersection of loci_dt and loci_from_files_dt
            # TODO: Need to use GenomicRanges' strand matching behaviour
            loci_dt <- loci_from_files_dt[loci_dt]
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("done in %.1f secs\n", stime))
            }
        }
        if (strandCollapse || rmZeroCov) {
            # NOTE: Have to update 'gr' if loci have been collapsed by strand
            #       or loci have been filtered to remove loci with zero
            #       coverage.
            if (verbose) {
                message("[read.bismark] Filtering 'gr' to retain valid loci ",
                        "...")
            }
            ptime1 <- proc.time()
            gr <- .lociDTAsGRanges(loci_dt, seqinfo)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("done in %.1f secs\n", stime))
            }
        }
    }

    # Construct 'M' and 'Cov' matrices
    if (verbose) {
        message("[read.bismark] Constructing 'M' and 'Cov' matrices ...")
    }
    ptime1 <- proc.time()
    counts <- .constructCountsFromBismarkFilesAndLociDT(
        files = files,
        loci_dt = loci_dt,
        strandCollapse = strandCollapse,
        mc.cores = mc.cores,
        BACKEND = BACKEND,
        verbose = subverbose)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }

    # Construct BSseq object
    if (verbose) {
        cat(sprintf("[read.bismark] Constructing BSseq object ... "))
    }
    ptime1 <- proc.time()
    # TODO: Use new_BSseq(), an internal, fast/minimal BSseq constructor
    #       (this function needs to be written).
    BSseq <- BSseq(
        M = counts[["M"]],
        Cov = counts[["Cov"]],
        gr = gr,
        sampleNames = sampleNames)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    BSseq
}

# TODOs ------------------------------------------------------------------------

# TODO: The documentation needs a complete overhaul and checked against Bismark
#       docs (e.g., the description of the .cov file is wrong)
# TODO: Consolidate use of message()/cat()/etc.
# TODO: mc.cores isn't used anywhere; can it be?
# TODO: Add function like minfi::read.metharray.sheet()?
# TODO: Should BACKEND really be an argument of read.bismark(); see related
#       issue on minfi repo https://github.com/hansenlab/minfi/issues/140
# TODO: May receive warning "In read_tokens_(data, tokenizer, col_specs, col_names,  ... : length of NULL cannot be changed". This is fixed in devel version of
#       readr (https://github.com/tidyverse/readr/issues/833)
