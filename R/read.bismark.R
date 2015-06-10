read.bismark <- function(files, sampleNames, rmZeroCov = FALSE, verbose = TRUE){
    ## Argument checking
    if (anyDuplicated(files)){
        stop("duplicate entries in 'files'")
    }
    if (length(sampleNames) != length(files) | anyDuplicated(sampleNames)){
        stop("argument 'sampleNames' has to have the same length as argument 'files', without duplicate entries")
    }
    ## Process each file
    idxes <- seq_along(files)
    names(idxes) <- sampleNames
    allOut <- lapply(idxes, function(ii){
        if (verbose) {
            cat(sprintf("[read.bismark] Reading file '%s' ... ", files[ii]))
        }
        ptime1 <- proc.time()
        raw <- read.bismarkFileRaw(thisfile = files[ii])
        M <- matrix(mcols(raw)[, "mCount"], ncol = 1)
        Cov <- M + mcols(raw)[, "uCount"]
        mcols(raw) <- NULL
        out <- BSseq(gr = raw, M = M, Cov = Cov,
                     sampleNames = sampleNames[ii], rmZeroCov = rmZeroCov)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))  
        }
        out  
    })
    if (verbose) {
        cat(sprintf("[read.bismark] Joining samples ... "))
    }
    ptime1 <- proc.time()
    allOut <- combineList(allOut)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    allOut
}

read.bismarkFileRaw <- function(thisfile, verbose = TRUE){
    ## Set up the 'what' argument for scan()
    columnHeaders <- c("chr", "start", "end", "mPerc", "mCount", "uCount")
    what0 <- replicate(length(columnHeaders), character(0))
    names(what0) <- columnHeaders
    int <- which(columnHeaders %in% c("start", "end", "mCount", "uCount"))
    what0[int] <- replicate(length(int), integer(0))
    null <- which(columnHeaders %in% "mPerc")
    what0[null] <- replicate(length(null), NULL)
    ## Read in the file
    if (grepl("\\.gz$", thisfile)) 
        con <- gzfile(thisfile)
    else 
        con <- file(thisfile, open = "r")
    out <- scan(file = con, what = what0, sep = "\t", quote = "", na.strings = "NA", quiet = TRUE)
    close(con)
    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = out[["chr"]], ranges = IRanges(start = out[["start"]], width = 1))
    out[["chr"]] <- out[["start"]] <- out[["end"]] <- NULL
    out <- out[!sapply(out, is.null)]
    df <- DataFrame(out)
    mcols(gr) <- df
    gr
}

read.bismarkCytosineRaw <- function(thisfile, keepContext = FALSE) {
    out <- fread(thisfile)
    if(length(out) != 6 && length(out) != 7)
        stop("unknown file format")
    if(length(out) == 6) {
        setnames(out, c("chr", "start", "end", "methPerc", "mCount", "uCount"))
        gr <- GRanges(seqnames = out[["chr"]],
                      ranges = IRanges(start = out[["start"]], width = 1),
                      mCount = out[["mCount"]], uCount = out[["uCount"]])
    }
    if(length(out) == 7) {
        setnames(out, c("chr", "start", "strand", "mCount", "uCount", "type", "context"))
        gr <- GRanges(seqnames = out[["chr"]], strand = out[["strand"]],
                      ranges = IRanges(start = out[["start"]], width = 1),
                      mCount = out[["mCount"]], uCount = out[["uCount"]])
        if(keepContext) {
            gr$type <- out[["type"]]
            gr$context <- out[["context"]]
        }                      
    }
    gr
}
