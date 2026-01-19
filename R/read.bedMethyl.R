# read.bedMethyl files into a extended BSseq object

#Internal function to read a single bedMethyl file and return a FWGRanges object
.readbedMethylAsFWGRanges <- function(file, rmZeroCov = TRUE,
                                     strandCollapse = TRUE, sort = TRUE,
                                     nThread = 1L, verbose = FALSE) {
  # Check that the provided arugments are boolean (logical) values
  # This ensures the function parameters are correctly set.
  stopifnot(isTRUEorFALSE(rmZeroCov)) #Note: bedMethyl files doesn't contain zero coverage loci and can be deprecated, but we can also keep it as a safeguard
  stopifnot(isTRUEorFALSE(strandCollapse))
  stopifnot(isTRUEorFALSE(sort))
  # Initialize all relevant modification-related columns as NULL
  M <- H <- U <- DE <- DI <- NO <- mod <- NULL
   #Note: This part can be removed if we deprecate rmZeroCov
   if (rmZeroCov) {
    # Read the bedMethyl data as a data.table object
    dt <- .readbedMethylAsDT(file = file,
                            col_spec = "BSseq",
                            check = TRUE,
                            nThread = nThread,
                            verbose = verbose)
    if (strandCollapse && !is.null(dt[["strand"]]) &&
        !dt[, all(strand == "*")]) {
      dt[strand == "-", `:=`(start, start - 1L)][, `:=`(strand, NULL)]
      dt <- unique(dt[, list(M = sum(M), H = sum(H), U= sum(U), D = sum(DE + DI + NO)),
                      by = c("seqnames", "start")])
    }
    dt <- dt[(M + H + U) > 0][, `:=`(c("M", "H", "U"),
                                     list(NULL, NULL, NULL))]
   #Strand collapse and sort the input
  } else {
    # Read the bedMethyl data as a data.table object
    dt <- .readbedMethylAsDT(file = file,
                            col_spec = "BSseq",
                            check = TRUE,
                            nThread = nThread,
                            verbose = verbose)
    # Perform strand collapsing and aggregation for the input data
    # If the strand information exists and is valid, adjust and merge data across strands
    if (strandCollapse && !is.null(dt[["strand"]]) &&!dt[, all(strand == "*")]) {
      # Adjust reverse ('-') strand coordinates and discard strand column.
      dt[strand == "-", `:=`(start, start - 1L)][, `:=`(strand, NULL)]
      # Aggregate modification data by sequence name (seqnames) and start position.
      dt <- unique(dt[, list(M = sum(M), H = sum(H), U = sum(U), D = sum(DE + DI + NO)),
                      by = c("seqnames", "start")])
    } else {
        # Perform aggregation while preserving strand information
      dt <- unique(dt[, list(M = sum(M), H = sum(H), U = sum(U), D = sum(DE + DI + NO)),
               by = c("seqnames", "start", "strand")])
    }
  }
  # Sort the data if sorting is requested
  # Group by sequence name and start position (consider strand if available)
  if (sort) {
    if (is.null(dt[["strand"]])) {
      setkey(dt, seqnames, start)
    } else {
      setkey(dt, seqnames, strand, start)
    }
  }
  #Prepare input for FWGRanges constructor
  seqnames <- Rle(dt[["seqnames"]])
  dt[, `:=`(seqnames, NULL)]
  seqinfo <- Seqinfo(seqnames = levels(seqnames))
  ranges <- .FWIRanges(start = dt[["start"]], width = 1L)
  dt[, `:=`(start, NULL)]
  mcols <- make_zero_col_DFrame(length(ranges))
  # Assign strand information, either as strand-agnostic ("*") or strand-specific.
  if (is.null(dt[["strand"]])) {
    strand <- strand(Rle("*", length(seqnames)))
  } else {
    strand <- Rle(dt[["strand"]])
    dt[, `:=`(strand, NULL)]
  }
  # Construct the FWGRanges object with relevant metadata
  fwgranges <- .FWGRanges(
    seqnames = seqnames,
    ranges = ranges,
    strand = strand,
    seqinfo = seqinfo,
    elementMetadata = mcols)
  if (sort) {
    fwgranges <- sort(sortSeqlevels(fwgranges))
  }
  fwgranges
}

# Function to read multiple bedMethyl files and construct a combined FWGRanges object
.constructFWGRangesFrombedMethylFiles <- function(files,
                                                 rmZeroCov,
                                                 strandCollapse,
                                                 verbose,
                                                 nThread,
                                                 BPPARAM) {
  subverbose <- max(as.integer(verbose) - 1L, 0L)
  #Write message with files being processed
  if (verbose) {
    message("[.constructFWGRangesFrombedMethylFiles] Extracting loci from: ")
    message(writeLines(c(files)), sep = "\n")
  }
  # Read and process loci from the first file
  # This serves as the base loci set for comparison with other files
  loci_from_first_file <- .readbedMethylAsFWGRanges(
    file = files[[1L]],
    rmZeroCov = rmZeroCov, #Note: This part can be removed if we deprecate rmZeroCov
    strandCollapse = strandCollapse,
    nThread = nThread,
    verbose = subverbose)
  # Configure parallel tasks for the remaining files
  if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
    bptasks(BPPARAM) <- length(files) - 1L
  }
  # Process all other files and identify loci not found in the first file
  list_of_loci_from_other_files_not_in_first_file <- bplapply(files[-1L],
    function(file, loci_from_first_file) {
    loci_from_this_file <- .readbedMethylAsFWGRanges(
      file = file,
      rmZeroCov = rmZeroCov, #Note: This part can be removed if we deprecate rmZeroCov
      strandCollapse = strandCollapse,
      verbose = subverbose)
    # Keep only loci that don't overlap with those in the first file
    subsetByOverlaps(
      x = loci_from_this_file,
      ranges = loci_from_first_file,
      type = "equal",
      invert = TRUE)
    }, loci_from_first_file = loci_from_first_file,
    BPPARAM = BPPARAM)
  # Combine loci from all files, ensuring unique loci across files
  loci_non_found_in_first_file <- unique(
    do.call(c, list_of_loci_from_other_files_not_in_first_file))
  loci <- c(loci_from_first_file, loci_non_found_in_first_file)
  # Sort all loci for consistent downstream processing
  sort(sortSeqlevels(loci))
}

## Internal function to read a bedMethyl file into a data.table
.readbedMethylAsDT <- function(file,
                              col_spec = c("all", "BSseq", "GRanges"),
                              check = FALSE,
                              showProgress = FALSE,
                              nThread = 1L,
                              verbose = verbose) {
  # Convert file path to an absolute path for robust file access
  file <- file_path_as_absolute(file)
  col_spec <- match.arg(col_spec)
  # Check that the provided arugments are boolean (logical) values
  stopifnot(isTRUEorFALSE(check))
  stopifnot(isTRUEorFALSE(showProgress))
  fread_verbose <- as.logical(max(verbose - 1L, 0L))

  header <- FALSE # you can include a header in bedMethyl file, so maybe add a check
  # Define the expected columns and their data types in the bedMethyl file:
  #  - "factor", "integer": Necessary columns (e.g., seqnames, start, modification data)
  #  - "NULL": Columns to be ignored during reading
  colClasses <- c("factor", "integer", "NULL", "factor", "NULL",
                  "character", "NULL", "NULL", "NULL", "NULL",
                  "NULL", "integer", "integer", "integer", "integer",
                  "NULL", "integer", "integer")
  # Specify the indices of unused columns to be dropped during file reading
  drop <- c(3L, 5L, 7L, 8L, 9L, 10L, 11L, 16L)
  # Specify the column names for the output data.table
  col.names <- c("seqnames", "start", "mod", "strand", "M", "U", "H", "DE", "DI", "NO")
  if (verbose) {
    message("[.readbedMethylAsDT] Parsing '", file, "'")
  }
  ptime1 <- proc.time()
  if (verbose) {
    message("[.readbedMethylAsDT] Reading file ...")
  }
  # Read the file into a data.table using fread (optimized for large files)
  x <- fread(input = file,
             sep = "\t",
             header = header,
             verbose = fread_verbose,
             drop = drop,
             colClasses = colClasses,
             col.names = col.names,
             quote = "",
             showProgress = showProgress,
             nThread = nThread)

  # Process the strand column:
  # If strand information exists and is "." (unknown), change it to "*" (strand-agnostic).
  if (!is.null(x[["strand"]])) {
    x[strand == ".", strand := "*"]
    x[, `:=`(strand, strand(strand))]
  }
  # Set the key for data.table operations (organize by the "mod" column)
  setkey(x, mod)
  # Filter rows to retain only modifications of interest (5mC and 5hmC)
  # This step and .check_bedMethyl ensures that only relevant modifications are included in downstream analyses.
  x<-x[.("m")] # remove all other mods except 5mC (they are still included in the count)
  if (verbose) {
    message("Reading in 5mC and 5hmC")
  }
  # Initialize all relevant modification-related columns as NULL
  M <- H <- U <- DE <- DI <- NO <- mod <- . <- NULL
  # Validate that all values in the count columns are non-negative integers
  if (check && all(c("M", "H", "U") %in% colnames(x))) {
    if (verbose) {
      message("[.readbedMethylAsDT] Checking validity of counts in file.")
    }
    valid <- x[, isTRUE(all(M >= 0L & U >= 0L & H >= 0L & DE >= 0L & DI >= 0L & NO >= 0L))]
    if (!valid) {
      stop("[.readbedMethylAsDT] Invalid counts detected.\n",
           "'M', 'H' and 'U' columns should be non-negative integers.")
    }
    else {
      if (verbose) {
        message("[.readbedMethylAsDT] All counts in file are valid.")
      }
    }
  }
  # Measure the time taken to process the file (and print it, if verbose)
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose) {
    message("Done in ", round(stime, 1), " secs")
  }
  # Return the processed data.table
  x
}


# Internal function to extract methylation counts from a single bedMethyl file
.constructCountsFromSinglebedMethylFile <- function(b, files, loci, strandCollapse,
                                                   grid, M_sink, H_sink, C_sink, D_sink,
                                                   Cov_sink, sink_lock,
                                                   nThread, verbose) {

  subverbose <- max(as.integer(verbose) - 1L, 0L)
  # Initialize all relevant modification-related columns as NULL
  M <- H <- U <- DE <- DI <- NO <- mod <- NULL

  # Get the file corresponding to the given index `b` in the files vector
  file <- files[b]
  if (verbose) {
    message("[.constructCountsFrombedMethylFileAndFWGRanges] Extracting ",
            "counts from ", "'", file, "'")
  }
  # Read the bedMethyl file into a data.table
  dt <- .readbedMethylAsDT(file = file, col_spec = "BSseq", check = TRUE,
                          nThread = nThread, verbose = subverbose)
  # Perform strand collapsing and aggregation for the input data
  if (strandCollapse && !is.null(dt[["strand"]]) &&
      !dt[, all(strand == "*")]) {
    dt[strand == "-", `:=`(start, start - 1L)][, `:=`(strand, NULL)]
    dt <- unique(dt[, list(M = sum(M), H = sum(H), U = sum(U), D = sum(DE + DI + NO)),
                    by = c("seqnames", "start")])
  } else {
    # Perform aggregation while preserving strand information
    dt <- unique(dt[, list(M = sum(M), H = sum(H), U = sum(U), D = sum(DE + DI + NO)),
                    by = c("seqnames", "start", "strand")])
  }
  # Convert the data.table into an FWGRanges object
  seqnames <- Rle(dt[["seqnames"]])
  dt[, `:=`(seqnames, NULL)]
  seqinfo <- Seqinfo(seqnames = levels(seqnames))
  ranges <- .FWIRanges(start = dt[["start"]], width = 1L)
  dt[, `:=`(start, NULL)]
  mcols <- make_zero_col_DFrame(length(ranges))
  if (is.null(dt[["strand"]])) {
    strand <- strand(Rle("*", length(seqnames)))
  }
  else {
    strand <- Rle(dt[["strand"]])
    dt[, `:=`(strand, NULL)]
  }
  # Construct an FWGRanges object for loci in this sample
  loci_from_this_sample <- .FWGRanges(
    seqnames = seqnames,
    ranges = ranges,
    strand = strand,
    seqinfo = seqinfo,
    elementMetadata = mcols)

  # Find overlaps between the sample-specific loci and the given loci
  ol <- findOverlaps(loci_from_this_sample, loci, type = "equal")

  # Initialize matrices for methylation metrics, each row corresponds to loci
  M <- matrix(rep(0L, length(loci)), ncol = 1)  # Methylated reads
  H <- matrix(rep(0L, length(loci)), ncol = 1)  # Hydroxymethylated reads
  U <- matrix(rep(0L, length(loci)), ncol = 1)  # Unmethylated reads
  D <- matrix(rep(0L, length(loci)), ncol = 1)  # Non-CpG coverage
  Cov <- matrix(rep(0L, length(loci)), ncol = 1)  # CpG coverage (M + H + U)

  # Fill the matrices with counts from the data.table for matching loci
  M[subjectHits(ol)] <- dt[queryHits(ol), ][["M"]]
  H[subjectHits(ol)] <- dt[queryHits(ol), ][["H"]]
  U[subjectHits(ol)] <- dt[queryHits(ol), ][["U"]]
  D[subjectHits(ol)] <- dt[queryHits(ol), ][["D"]]
  Cov[subjectHits(ol)] <- dt[queryHits(ol), list(Cov = (M + H + U))][["Cov"]]

  # If no sink is provided, return matrices as list
  if (is.null(M_sink)) {
    return(list(M = M, H = H, U = U, D = D, Cov = Cov))
  }
  # Write the matrices into sinks using locks (thread-safe parallel writes)
  viewport <- grid[[b]]
  ipclock(sink_lock)
  write_block(M_sink, viewport = viewport, block = M)
  write_block(H_sink, viewport = viewport, block = H)
  write_block(C_sink, viewport = viewport, block = U)
  write_block(D_sink, viewport = viewport, block = D)
  write_block(Cov_sink, viewport = viewport, block = Cov)
  ipcunlock(sink_lock)
  NULL
}

# Internal function to process multiple bedMethyl files and construct methylation counts
.constructbedMethylCounts <- function(files, loci, strandCollapse, BPPARAM,
                            BACKEND, dir, chunkdim, level, nThread,
                            verbose = FALSE) {

  subverbose <- max(as.integer(verbose) - 1L, 0L)
  # Dimensions of the output matrices
  ans_nrow <- length(loci) # Number of loci
  ans_ncol <- length(files) # Number of files (samples)
  ans_dim <- c(ans_nrow, ans_ncol)
  # Create a grid for dividing the output matrix into chunks for parallel processing
  grid <- RegularArrayGrid(refdim = ans_dim, spacings = c(ans_nrow,1L))
  # Initialize sink for outputs based on the storage backend
  if (is.null(BACKEND)) {
    M_sink <- NULL
    H_sink <- NULL
    C_sink <- NULL
    D_sink <- NULL
    Cov_sink <- NULL
    sink_lock <- NULL
  }
  else if (BACKEND == "HDF5Array") {
    # Set up HDF5 storage for efficient on-disk realization
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
    H_sink <- HDF5RealizationSink(
      dim = ans_dim,
      dimnames = NULL,
      type = "integer",
      filepath = h5_path,
      name = "H",
      chunkdim = chunkdim,
      level = level)
    on.exit(close(H_sink), add = TRUE)
    C_sink <- HDF5RealizationSink(
      dim = ans_dim,
      dimnames = NULL,
      type = "integer",
      filepath = h5_path,
      name = "U",
      chunkdim = chunkdim,
      level = level)
    on.exit(close(C_sink), add = TRUE)
    D_sink <- HDF5RealizationSink(
      dim = ans_dim,
      dimnames = NULL,
      type = "integer",
      filepath = h5_path,
      name = "D",
      chunkdim = chunkdim,
      level = level)
    on.exit(close(D_sink), add = TRUE)
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
    # Set up DelayedArray sink for memory-based data realization
    M_sink <- DelayedArray::AutoRealizationSink(
      dim = ans_dim,
      type = "integer")
    on.exit(close(M_sink), add = TRUE)
    H_sink <- DelayedArray::AutoRealizationSink(
      dim = ans_dim,
      type = "integer")
    on.exit(close(H_sink), add = TRUE)
    C_sink <- DelayedArray::AutoRealizationSink(
      dim = ans_dim,
      type = "integer")
    on.exit(close(C_sink), add = TRUE)
    D_sink <- DelayedArray::AutoRealizationSink(
      dim = ans_dim,
      type = "integer")
    on.exit(close(D_sink), add = TRUE)
    Cov_sink <- DelayedArray::AutoRealizationSink(
      dim = ans_dim,
      type = "integer")
    on.exit(close(Cov_sink), add = TRUE)
    sink_lock <- ipcid()
    on.exit(ipcremove(sink_lock), add = TRUE)
  }
  if (is(BPPARAM, "SnowParam") && bpprogressbar(BPPARAM)) {
    bptasks(BPPARAM) <- length(grid)
  }
  # Process each chunk of the grid in parallel using bplapply
  counts <- bptry(bplapply(X = seq_along(grid),
                           FUN = .constructCountsFromSinglebedMethylFile,
                           files = files,
                           loci = loci,
                           strandCollapse = strandCollapse,
                           grid = grid,
                           M_sink = M_sink,
                           H_sink = H_sink,
                           C_sink = C_sink,
                           D_sink = D_sink,
                           Cov_sink = Cov_sink,
                           sink_lock = sink_lock,
                           verbose = subverbose,
                           nThread = nThread,
                           BPPARAM = BPPARAM))
  # Check for errors during parallel runs and handle them
  if (!all(bpok(counts))) {
    stop(".constructbedMethylCounts() encountered errors for these files:\n  ",
         paste(files[!bpok], collapse = "\n  "))
  }
  # If no storage backend is provided, gather results into in-memory matrices
  if (is.null(BACKEND)) {
    M <- do.call(c, lapply(counts, "[[", "M"))
    attr(M, "dim") <- ans_dim
    H <- do.call(c, lapply(counts, "[[", "H"))
    attr(H, "dim") <- ans_dim
    U <- do.call(c, lapply(counts, "[[", "U"))
    attr(U, "dim") <- ans_dim
    D <- do.call(c, lapply(counts, "[[", "D"))
    attr(D, "dim") <- ans_dim
    Cov <- do.call(c, lapply(counts, "[[", "Cov"))
    attr(Cov, "dim") <- ans_dim
  } else {
    # Activate on-disk storage backend (e.g., HDF5 or DelayedArray)
    # Convert the HDF5 or DelayedArray sinks (initialized earlier) into actual array-like objects
    M <- as(M_sink, "DelayedArray")
    H <- as(H_sink, "DelayedArray")
    U <- as(C_sink, "DelayedArray")
    D <- as(D_sink, "DelayedArray")
    Cov <- as(Cov_sink, "DelayedArray")
  }
  return(list(M = M, H = H, U = U, D = D, Cov = Cov))
}

##
read.bedMethyl <- function (files,
                            loci = NULL,
                            colData = NULL,
                            rmZeroCov = TRUE,
                            strandCollapse = TRUE,
                            BPPARAM = bpparam(),
                            BACKEND = NULL,
                            dir = tempfile("BSseq"),
                            replace = FALSE,
                            chunkdim = NULL,
                            level = NULL,
                            nThread = 1L,
                            output = c("BSseq", "MethylCounts"),
                            check_input = TRUE,
                            verbose = getOption("verbose")) {
  # Choose the desired output format: "BSseq" or "MethylCounts"
  output <- match.arg(output)
  # Ensure no duplicate file names are present (duplicates would cause ambiguity).
  if (anyDuplicated(files)) {
    stop("'files' cannot have duplicate entries.")
  }
  # Check that all specified files exist. Abort if any files are missing.
  file_exists <- file.exists(files)
  if (!isTRUE(all(file_exists))) {
    stop("These files cannot be found:\n  ", paste(files[!file_exists],
                                                   collapse = "\n  "))
  }
  # Optionally validate each bedMethyl file (if `check_input` is enabled).
  # This step calls `.check_bedMethyl` to ensure the file has the expected format and structure.
  if (check_input) {
   for (file in files) {
      message("Validating file: ", file)
      .check_bedMethyl(file = file, output = output)
   }
  }
  # Verify the validity of the `loci` argument (if provided)
  if (!is.null(loci)) {
    if (!is(loci, "GenomicRanges")) {
      stop("'loci' must be a GenomicRanges instance if not NULL.")
    }
    if (any(width(loci) != 1L)) {
      stop("All elements of 'loci' must have width equal to 1.")
    }
  }
  # Create default colData if none is provided
  if (is.null(colData)) {
    colData <- DataFrame(row.names = files)
  }
  # Ensure colData rows match the number of files
  if (nrow(colData) != length(files)) {
    stop("Supplied 'colData' must have nrow(colData) == length(files).")
  }
  # Validate logical arguments: `rmZeroCov` and `strandCollapse`.
  stopifnot(isTRUEorFALSE(rmZeroCov))
  stopifnot(isTRUEorFALSE(strandCollapse))
  # Handle realization backends (in-memory or on-disk).
  current_BACKEND <- getAutoRealizationBackend()
  on.exit(setAutoRealizationBackend(current_BACKEND), add = TRUE)
  setAutoRealizationBackend(BACKEND)

  # Validate backend compatibility (e.g., in-memory vs. on-disk). #Should be updated to read.bedMethyl specifics
  if (!.areBackendsInMemory(BACKEND)) {
    if (!.isSingleMachineBackend(BPPARAM)) {
      stop("The parallelisation strategy must use a single machine ",
           "when using an on-disk realization backend.\n",
           "See help(\"read.bismark\") for details.", call. = FALSE)
    }
  }
  else {
    if (!is.null(BACKEND)) {
      stop("The '", BACKEND, "' realization backend is not supported.",
           "\n  See help(\"read.bismark\") for details.",
           call. = FALSE)
    }
  }
  if (identical(BACKEND, "HDF5Array")) {
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
    }
    else {
      HDF5Array::replace_dir(dir, replace)
    }
  }
  subverbose <- as.logical(max(verbose - 1L, 0L))
  if (is.null(loci)) {
    ptime1 <- proc.time()
    if (verbose) {
      message("[read.bedMethyl] Parsing files and constructing valid loci ...")
    }
    loci <- .constructFWGRangesFrombedMethylFiles(
      files = files,
      rmZeroCov = rmZeroCov,
      strandCollapse = strandCollapse,
      verbose = verbose,
      nThread = nThread,
      BPPARAM = BPPARAM)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
      message("Done in ", round(stime, 1), " secs")
    }
  }
  else {
    if (verbose) {
      message("[read.bedMethyl] Using 'loci' as candidate loci.")
    }
    if (strandCollapse) {
      if (verbose) {
        message("[read.bedMethyl] Collapsing strand of 'loci' ...")
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
        message("[read.bedMethyl] Parsing files to identify elements of ",
                "'loci' with non-zero coverage ...")
      }
      ptime1 <- proc.time()
      loci_from_files <- .constructFWGRangesFrombedMethylFiles(
        files = files,
        rmZeroCov = rmZeroCov,
        strandCollapse = strandCollapse,
        verbose = subverbose,
        nThread = nThread,
        BPPARAM = BPPARAM)
      loci <- subsetByOverlaps(loci, loci_from_files, type = "equal")
      ptime2 <- proc.time()
      stime <- (ptime2 - ptime1)[3]
      if (verbose) {
        message("Done in ", round(stime, 1), " secs")
      }
    }
  }
  ptime1 <- proc.time()
  if (verbose) {
    message("[read.bedMethyl] Parsing files and constructing",
            "'M', 'H', 'U', 'D' and 'Cov' ", "matrices ...")
  }
  counts <- .constructbedMethylCounts(
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
  if (verbose) {
      message("[read.bedMethyl] Constructing ", if(output=="MethylCounts") "MethylCounts" else "BSseq", " object ... ")
  }
  if (output == "MethylCounts") {
      assays_mc <- list(M = counts$M, H = counts$H, U = counts$U, D = counts$D)
      se <- SummarizedExperiment(assays = assays_mc, rowRanges = shift(as(loci, "GRanges"),1), colData = colData)
      mc <- new2("MethylCounts", se, check = FALSE)
      return(mc)
  }
  else {
  se <- SummarizedExperiment(assays = counts, rowRanges = shift(as(loci,"GRanges"),1), colData = colData)
  bsseq <- new2("BSseq", se, check = FALSE)
  if (!is.null(BACKEND) && BACKEND == "HDF5Array") {
    x <- bsseq
    x@assays <- HDF5Array::shorten_assay2h5_links(x@assays)
    saveRDS(x, file = file.path(dir, "se.rds"))
  }
  bsseq
  }
}
