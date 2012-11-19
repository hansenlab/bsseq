read.bismark <- function(files, sampleNames, returnRaw = FALSE, rmZeroCov = FALSE, verbose = TRUE){
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
      cat(sprintf("Reading file '%s' ... ", files[ii]))
    }
    stime <- system.time({
      if (returnRaw) {
        out <- read.bismarkFileRaw(thisfile = files[ii])
      } else{
        raw <- read.bismarkFileRaw(thisfile = files[ii])
        M <- matrix(elementMetadata(raw)[, "mCount"], ncol = 1)
        Cov <- M + elementMetadata(raw)[, "uCount"]
        elementMetadata(raw) <- NULL
        out <- BSseq(gr = raw, M = M, Cov = Cov, sampleNames = sampleNames[ii], rmZeroCov = rmZeroCov)
      }
    })[3]
    if (verbose) {
      cat(sprintf("in %.1f secs\n", stime))  
    }
    out  
  })
  if (!returnRaw) {
    if (verbose) {
      cat(sprintf("Joining samples ... "))
    }
    stime <- system.time({
      allOut <- Reduce("combine", allOut)
    })[3]
    if (verbose) {
      cat(sprintf("in %.1f secs\n", stime))
    }
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
  dbl <- which(columnHeaders %in% "mPerc")
  what0[dbl] <- replicate(length(dbl), double(0))
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
  elementMetadata(gr) <- df
  gr
}
