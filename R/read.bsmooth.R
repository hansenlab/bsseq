read.umtab2 <- function(dirs, sampleNames = NULL, rmZeroCov = FALSE,
                        readCycle = FALSE, keepFilt = FALSE,
                        pattern = NULL, keepU, keepM, verbose = TRUE) {
    filesPerDir <- lapply(dirs, function(xx) sort(list.files(xx, pattern = pattern, full.names = TRUE)))
    if(!all(sapply(filesPerDir, function(xx) all(basename(xx) == basename(filesPerDir[[1]])))))
        warning("'dirs' does not contain the exact same file names")
    allChrs <- as.list(seq_along(filesPerDir[[1]]))
    for(ii in seq_along(filesPerDir[[1]])) {
        chrRead <- read.umtab2.chr(sapply(filesPerDir, function(xx) xx[ii]),
                                   keepM = keepM, keepU = keepU, sampleNames = sampleNames,
                                   readCycle = readCycle, keepFilt = keepFilt,
                                   verbose = verbose)
        if(rmZeroCov && length(chrRead) != 0){
            wh <- which(rowSums2(chrRead$U + chrRead$M) > 0)
            chrRead$M <- chrRead$M[wh,]
            chrRead$U <- chrRead$U[wh,]
            if(readCycle) {
                chrRead$Mcy <- chrRead$Mcy[wh,]
                chrRead$Ucy <- chrRead$Ucy[wh,]
            }
            chrRead$chr <- chrRead$chr[wh]
            chrRead$pos <- chrRead$pos[wh]
            if(keepFilt) {
                chrRead$filt_cycle <- chrRead$filt_cycle[wh,]
                chrRead$filt_allele <- chrRead$filt_allele[wh,]
                chrRead$filt_mapq <- chrRead$filt_mapq[wh,]
                chrRead$filt_baseq <- chrRead$filt_baseq[wh,]
            }
        }
        allChrs[[ii]] <- chrRead
        print(gc())
    }
    chr <- do.call(c, lapply(allChrs, function(xx) xx$chr))
    pos <- do.call(c, lapply(allChrs, function(xx) xx$pos))
    M <- do.call(rbind, lapply(allChrs, function(xx) xx$M))
    Cov <- do.call(rbind, lapply(allChrs, function(xx) xx$M)) +
        do.call(rbind, lapply(allChrs, function(xx) xx$U))
    if(readCycle) {
        Mcy <- do.call(rbind, lapply(allChrs, function(xx) xx$Mcy))
        Ucy <- do.call(rbind, lapply(allChrs, function(xx) xx$Ucy))
    } else {
        Mcy <- NULL
        Ucy <- NULL
    }
    if(keepFilt) {
        filt_cycle <- do.call(rbind, lapply(allChrs, function(xx) xx$filt_cycle))
        filt_allele <- do.call(rbind, lapply(allChrs, function(xx) xx$filt_allele))
        filt_mapq <- do.call(rbind, lapply(allChrs, function(xx) xx$filt_mapq))
        filt_baseq <- do.call(rbind, lapply(allChrs, function(xx) xx$filt_baseq))
    } else {
        filt_cycle <- filt_allele <- filt_mapq <- filt_baseq <- NULL
    }
    print(gc())
    csums <- lapply(allChrs, function(xx) xx$csums)
    names(csums) <- basename(filesPerDir[[1]])
    rm(allChrs)
    BSdata <- BSseq(chr = chr, pos = pos, M = M, Cov = Cov, rmZeroCov = rmZeroCov)
    print(gc())
    return(list(BSdata = BSdata, filt_cycle = filt_cycle, filt_allele = filt_allele,
                filt_mapq = filt_mapq, filt_baseq = filt_baseq,
                Mcy = Mcy, Ucy = Ucy, csums = csums))
}


read.umtab2.chr <- function(files, sampleNames = NULL,
                            keepM, keepU, readCycle = FALSE, keepFilt = FALSE,
                            verbose = TRUE) {
    columnHeaders <- strsplit(readLines(files[1], n = 1), "\t")[[1]]
    if("strand" %in% columnHeaders)
        stranded <- TRUE
    what0 <- replicate(length(columnHeaders), integer(0))
    names(what0) <- columnHeaders
    what0[["ref"]] <- character(0)
    if(stranded)
        what0[["strand"]] <- character(0)
    if(!readCycle) {
        what0 <- what0[! what0 %in% c("Mcy", "Ucy")]
    }
    if(!keepFilt) {
        what0 <- what0[!grepl("^filt", what0)]
    }
    if(verbose) cat("[read.umtab2.chr] reading", files[1], "\n")
    scanPars <- list(sep = "", quote = "", quiet = TRUE, skip = 1,
                     what = what0, na.strings = c("NA", "?"))
    intab <- do.call(scan, c(scanPars, file = files[1]))
    if(length(intab[[1]]) == 0)
        return(list())
    pos <- intab[["off"]]
    chr <- intab[["ref"]]
    nPos <- length(pos)
    nSamples <- length(files)
    M <- matrix(0L, nrow = nPos, ncol = nSamples)
    colnames(M)  <- sampleNames
    U <- M
    if(readCycle)
        Mcy <- Ucy <- M
    else
        Mcy <- Ucy <- NULL
    if(keepFilt)
        filt_cycle <- filt_allele <- filt_mapq <- filt_baseq <- M
    else
        filt_cycle <- filt_allele <- filt_mapq <- filt_baseq <- NULL
    allMnames <- grep("^M[[:digit:]][[:digit:]]*", names(intab), value = TRUE)
    allUnames <- grep("^U[[:digit:]][[:digit:]]*", names(intab), value = TRUE)
    csums <- matrix(0, nrow = length(c(allMnames, allUnames)), ncol = length(files))
    colnames(csums) <- sampleNames
    rownames(csums) <- c(allMnames, allUnames)
    if(missing(keepM)) keepM <- allMnames
    if(missing(keepU)) keepU <- allUnames
    stopifnot(all(c(keepM, keepU) %in% c(allMnames, allUnames)))
    M[,1] <- as.integer(Reduce("+", intab[keepM]))
    U[,1] <- as.integer(Reduce("+", intab[keepU]))
    csums[,1] <- as.integer(sapply(intab[c(allMnames, allUnames)], sum))
    if(readCycle) {
        Mcy[,1] <- intab[["Mcy"]]
        Ucy[,1] <- intab[["Ucy"]]
    }
    if(keepFilt) {
        filt_cycle[,1] <- intab[["filt_cycle"]]
        filt_allele[,1] <- intab[["filt_allele"]]
        filt_mapq[,1] <- intab[["filt_mapq"]]
        filt_baseq[,1] <- intab[["filt_baseq"]]
    }
    for(ii in seq_along(files[-1]) + 1) {
        if(verbose) cat("reading", files[ii], "\n")
        intab <- do.call(scan, c(scanPars, file = files[ii]))
        if(length(intab[[1]]) == 0)
            next
        stopifnot(all(pos == intab[["Off"]]) && all(chr == intab[["Chr"]]))
        M[,ii] <- as.integer(Reduce("+", intab[keepM]))
        U[,ii] <- as.integer(Reduce("+", intab[keepU]))
        csums[,ii] <- as.integer(sapply(intab[c(allMnames, allUnames)], sum))
        if(readCycle) {
            Mcy[,ii] <- intab[["Mcy"]]
            Ucy[,ii] <- intab[["Ucy"]]
        }
        if(keepFilt) {
            filt_cycle[,ii] <- intab[["filt_cycle"]]
            filt_allele[,ii] <- intab[["filt_allele"]]
            filt_mapq[,ii] <- intab[["filt_mapq"]]
            filt_baseq[,ii] <- intab[["filt_baseq"]]
        }
    }
    return(list(chr = chr, pos = pos, M = M, U = U,
                Mcy = Mcy, Ucy = Ucy,
                filt_cycle = filt_cycle, filt_allele = filt_allele,
                filt_mapq = filt_mapq, filt_baseq = filt_baseq,
                csums = csums))
}

read.umtab2.chr2 <- function(files, sampleNames = NULL,
                             keepM, keepU, verbose = TRUE) {
    columnHeaders <- strsplit(readLines(files[1], n = 1), "\t")[[1]]
    if("strand" %in% columnHeaders)
        stranded <- TRUE
    what0 <- replicate(length(columnHeaders), integer(0))
    names(what0) <- columnHeaders
    what0[["ref"]] <- character(0)
    if(stranded)
        what0[["strand"]] <- character(0)
    if(verbose) cat("[read.umtab2.chr] parsing all samples first\n")
    scanPars <- list(sep = "", quote = "", quiet = TRUE, skip = 1,
                     what = what0, na.strings = c("NA", "?"))
    ## First we get the locations of everything in all files
    getLocs <- function(file) {
        if(stranded) {
            intab <- scan(file, what = list("ref" = character(0),
                                "off" = integer(0),
                                "strand" = character(0)),
                          flush = TRUE, skip = 1, quote = "", sep = "")
            gr <- GRanges(seqnames = intab[["ref"]],
                          strand = ifelse(intab[["strand"]]== "W", 1L, -1L),
                          ranges = IRanges(start = intab[["off"]], width = 1))
        } else {
            intab <- scan(file, what = list("ref" = character(0),
                                "off" = integer(0)),
                          flush = TRUE, quote = "", sep = "")
            gr <- GRanges(seqnames = intab[["ref"]],
                          ranges = IRanges(start = intab[["pos"]], width = 1))
        }
        gr
    }
    grLoc <- as.list(files)
    grLoc[[1]] <- getLocs(files[1])
    grLoc <- Reduce(function(gr, file) {
        union(gr, getLocs(file))
    }, grLoc)
    nPos <- length(grLoc)
    nSamples <- length(files)

    M <- matrix(0L, nrow = nPos, ncol = nSamples)
    colnames(M)  <- sampleNames
    U <- Mcy <- Ucy <- M
    filt_cycle <- filt_allele <- filt_mapq <- filt_baseq <- M
    allMnames <- grep("^M[[:digit:]][[:digit:]]*", columnHeaders, value = TRUE)
    allUnames <- grep("^U[[:digit:]][[:digit:]]*", columnHeaders, value = TRUE)
    csums <- matrix(0, nrow = length(c(allMnames, allUnames)), ncol = length(files))
    colnames(csums) <- sampleNames
    rownames(csums) <- c(allMnames, allUnames)
    if(missing(keepM)) keepM <- allMnames
    if(missing(keepU)) keepU <- allUnames
    stopifnot(all(c(keepM, keepU) %in% c(allMnames, allUnames)))

    for(ii in seq_along(files)) {
        if(verbose) cat("[read.umtab2.chr] reading", files[ii], "\n")
        intab <- do.call(scan, c(scanPars, file = files[ii]))
        if(length(intab[[1]]) == 0)
            next
        if(stranded) {
            gr <- GRanges(seqnames = intab[["ref"]],
                          strand = ifelse(intab[["strand"]] == "W", 1L, -1L),
                          ranges = IRanges(start = intab[["off"]], width = 1))
        } else {
            gr <- GRanges(seqnames = intab[["ref"]],
                          ranges = IRanges(start = intab[["off"]], width = 1))
        }
        sh <- subjectHits(findOverlaps(gr, grLoc))
        stopifnot(length(sh) == length(intab$ref) & !anyDuplicated(sh))
        M[sh,ii] <- as.integer(Reduce("+", intab[keepM]))
        U[sh,ii] <- as.integer(Reduce("+", intab[keepU]))
        Mcy[sh,ii] <- intab[["Mcy"]]
        Ucy[sh,ii] <- intab[["Ucy"]]
        filt_cycle[sh,ii] <- intab[["filt_cycle"]]
        filt_allele[sh,ii] <- intab[["filt_allele"]]
        filt_mapq[sh,ii] <- intab[["filt_mapq"]]
        filt_baseq[sh,ii] <- intab[["filt_baseq"]]
        csums[,ii] <- as.integer(sapply(intab[c(allMnames, allUnames)], sum))
    }
    return(list(chr = as.character(seqnames(grLoc)), pos = start(grLoc), M = M, U = U,
                Mcy = Mcy, Ucy = Ucy,
                filt_cycle = filt_cycle, filt_allele = filt_allele,
                filt_mapq = filt_mapq, filt_baseq = filt_baseq,
                csums = csums))
}



read.umtab <- function(dirs, sampleNames = NULL, rmZeroCov = FALSE,
                       pattern = NULL,
                       keepU = c("U10", "U20", "U30", "U40"),
                       keepM = c("M10", "M20", "M30", "M40"), verbose = TRUE) {
    filesPerDir <- lapply(dirs, function(xx) sort(list.files(xx, pattern = pattern, full.names = TRUE)))
    if(!all(sapply(filesPerDir, function(xx) all(basename(xx) == basename(filesPerDir[[1]])))))
        warning("'dirs' does not contain the exact same file names")
    allChrs <- as.list(seq_along(filesPerDir[[1]]))
    for(ii in seq_along(filesPerDir[[1]])) {
        chrRead <- read.umtab.chr(sapply(filesPerDir, function(xx) xx[ii]),
                                  keepM = keepM, keepU = keepU, sampleNames = sampleNames,
                                  verbose = verbose)
        if(rmZeroCov && length(chrRead) != 0){
            wh <- which(rowSums2(chrRead$U + chrRead$M) > 0)
            chrRead$M <- chrRead$M[wh,]
            chrRead$U <- chrRead$U[wh,]
            chrRead$Mcy <- chrRead$Mcy[wh,]
            chrRead$Ucy <- chrRead$Ucy[wh,]
            chrRead$chr <- chrRead$chr[wh]
            chrRead$pos <- chrRead$pos[wh]
            chrRead$Map <- chrRead$Map[wh]
            chrRead$GC <- chrRead$GC[wh]
        }
        allChrs[[ii]] <- chrRead
    }
    BSdata <- BSseq(chr = do.call(c, lapply(allChrs, function(xx) xx$chr)),
                    pos = do.call(c, lapply(allChrs, function(xx) xx$pos)),
                    M = do.call(rbind, lapply(allChrs, function(xx) xx$M)),
                    Cov = do.call(rbind, lapply(allChrs, function(xx) xx$M)) +
                    do.call(rbind, lapply(allChrs, function(xx) xx$U)),
                    rmZeroCov = rmZeroCov)
    GC <- do.call(c, lapply(allChrs, function(xx) xx$GC))
    Map <- do.call(c, lapply(allChrs, function(xx) xx$Map))
    Mcy <- do.call(rbind, lapply(allChrs, function(xx) xx$Mcy))
    Ucy <- do.call(rbind, lapply(allChrs, function(xx) xx$Ucy))
    csums <- lapply(allChrs, function(xx) xx$csums)
    names(csums) <- basename(filesPerDir[[1]])
    return(list(BSdata = BSdata, GC = GC, Map = Map,
                Mcy = Mcy, Ucy = Ucy, csums = csums))
}


read.umtab.chr <- function(files, sampleNames = NULL,
                           keepM = c("M10", "M20", "M30", "M40"),
                           keepU = c("U10", "U20", "U30", "U40"), verbose = TRUE) {
    columnHeaders <- c("Chr", "Off", "M0", "M10", "M20", "M30", "M40", "Mqual", "Mcy",
                       "U0", "U10", "U20", "U30", "U40", "Uqual", "Ucy",
                       "MapaFW", "MapaRC", "CGContFW", "CGContRC")
    what0 <- c(list(character(0)), replicate(19, integer(0)))
    names(what0) <- columnHeaders
    line1 <- readLines(files[1], n = 1)
    if(length(line1) == 0)
        return(list())
    line1 <- strsplit(line1, "\t")[[1]]
    stopifnot(all(line1 == columnHeaders))
    if(verbose) cat("[read.umtab.chr] reading", files[1], "\n")
    scanPars <- list(sep = "", quote = "", quiet = TRUE, skip = 1,
                  what = what0, na.strings = c("NA", "?"))
    intab <- do.call(scan, c(scanPars, file = files[1]))
    if(length(intab[[1]]) == 0)
        return(list())
    GC <- as.integer(intab[["CGContFW"]] + intab[["CGContRC"]])
    Map <- as.integer(intab[["MapaFW"]] + intab[["MapaRC"]])
    pos <- intab[["Off"]]
    chr <- intab[["Chr"]]
    nPos <- length(pos)
    nSamples <- length(files)
    M <- matrix(0L, nrow = nPos, ncol = nSamples)
    U <- matrix(0L, nrow = nPos, ncol = nSamples)
    Mcy <- matrix(0L, nrow = nPos, ncol = nSamples)
    Ucy <- matrix(0L, nrow = nPos, ncol = nSamples)
    csums <- matrix(0, nrow = 10, ncol = length(files))
    colnames(M) <- colnames(U) <- colnames(csums) <- sampleNames
    allMnames <- c("M0", "M10", "M20", "M30", "M40")
    allUnames <- c("U0", "U10", "U20", "U30", "U40")
    stopifnot(all(c(keepM, keepU) %in% c(allMnames, allUnames)))
    rownames(csums) <- c(allMnames, allUnames)
    M[,1] <- as.integer(Reduce("+", intab[keepM]))
    U[,1] <- as.integer(Reduce("+", intab[keepU]))
    Mcy[,1] <- intab[["Mcy"]]
    Ucy[,1] <- intab[["Ucy"]]
    csums[,1] <- as.integer(sapply(intab[c(allMnames, allUnames)], sum))
    for(ii in seq_along(files[-1]) + 1) {
        if(verbose) cat("[read.umtab2] reading", files[ii], "\n")
        intab <- do.call(scan, c(scanPars, file = files[ii]))
        if(length(intab[[1]]) == 0)
            next
        stopifnot(all(pos == intab[["Off"]]) && all(chr == intab[["Chr"]]))
        M[,ii] <- as.integer(Reduce("+", intab[keepM]))
        U[,ii] <- as.integer(Reduce("+", intab[keepU]))
        Mcy[,ii] <- intab[["Mcy"]]
        Ucy[,ii] <- intab[["Ucy"]]
        csums[,ii] <- as.integer(sapply(intab[c(allMnames, allUnames)], sum))
    }
    return(list(chr = chr, pos = pos, Map = Map, GC = GC, M = M, U = U,
                Mcy = Mcy, Ucy = Ucy, csums = csums))
}

read.bsmoothDirRaw <- function(dir, seqnames = NULL, keepCycle = FALSE, keepFilt = FALSE,
                               header = TRUE, verbose = TRUE) {
    dir <- normalizePath(dir)
    inpattern <- "\\.cpg\\.tsv(|\\.gz)$"
    if(length(dir) != 1 || !file.info(dir)$isdir)
        stop("argument 'dir' needs to be a single directory")
    allChrFiles <- list.files(dir, pattern = inpattern, full.names = TRUE)
    if(!is.null(seqnames))
        allChrFiles <- allChrFiles[sub(inpattern, "", basename(allChrFiles)) %in% seqnames]
    if(length(allChrFiles) == 0) {
        warning(sprintf("dir '%s' is empty or has no output from 'seqnames'", dir))
        return(NULL)
    }
    if(header) {
        columnHeaders <- sapply(allChrFiles, function(thisfile) {
            if(grepl("\\.gz$", thisfile))
                con <- gzfile(thisfile)
            else
                con <- file(thisfile, open = "r")
            out <- readLines(con, n = 1)
            close(con)
            out
        })
        columnHeaders <- strsplit(columnHeaders, "\t")
        if(!all(sapply(columnHeaders, function(xx) all.equal(columnHeaders[[1]], xx))))
            stop(sprintf("input files in dir '%s' does not have the same headers", dir))
        columnHeaders <- columnHeaders[[1]]
    } else
        columnHeaders <- c("ref", "off", "strand", "Mstr", "Mcy", "Ustr", "Ucy",
                           "filt_cycle", "filt_readlen", "filt_allele", "filt_mapq", "filt_baseq")
    what0 <- replicate(length(columnHeaders), character(0))
    names(what0) <- columnHeaders
    int <- c(which(columnHeaders %in% c("off", "Mcy", "Ucy")), grep("^filt", columnHeaders))
    what0[int] <- replicate(length(int), integer(0))
    if(!keepCycle)
        what0[c("Mcy", "Ucy")] <- replicate(2, NULL)
    if(!keepFilt)
        what0[grep("^filt", names(what0))] <- replicate(length(grep("^filt", names(what0))), NULL)
    outList <- lapply(allChrFiles, function(thisfile) {
        if(verbose)
            cat(sprintf("[read.bsmoothDirRaw] Reading '%s'\n", thisfile))
        if(grepl("\\.gz$", thisfile))
            con <- gzfile(thisfile)
        else
            con <- file(thisfile)
        out <- scan(con, skip = header, what = what0, sep = "\t",
                    quote = "", na.strings = "NA", quiet = TRUE)
        close(con)
        out
    })
    listNames <- names(outList[[1]])
    names(listNames) <- listNames
    out <- lapply(listNames, function(name) {
        do.call(c, lapply(outList, function(xx) xx[[name]]))
    })
    rm(outList)
    gr <- GRanges(seqnames = paste0("chr", out[["ref"]]),
                  ranges = IRanges(start = out[["off"]], width = 1))
    out[["ref"]] <- out[["off"]] <- NULL
    names(out)[names(out) == "strand"] <- "bstrand"
    out <- out[!sapply(out, is.null)]
    df <- DataFrame(out)
    mcols(gr) <- df
    gr
}


sampleRawToBSseq <- function(gr, qualityCutoff = 20, sampleName = NULL, rmZeroCov = FALSE) {
    numberQualsGreaterThan <- function(cvec) {
        onestring <- paste(cvec, collapse = "")
        greater <- (as.integer(charToRaw(onestring)) - 33L >= qualityCutoff)
        out <- tapply(greater, rep(1:length(cvec), times = nchar(cvec)), sum)
        out
    }
    strToCov <- function(vec) {
        Cov <- rep(0, length(vec))
        wh <- which(! vec %in% c("", "0"))
        if(length(wh) > 0)
            Cov[wh] <- numberQualsGreaterThan(vec[wh])
        Cov
    }
    M <- matrix(strToCov(mcols(gr)[, "Mstr"]), ncol = 1)
    Cov <- M + strToCov(mcols(gr)[, "Ustr"])
    mcols(gr) <- NULL
    BSseq(gr = gr, M = M, Cov = Cov, sampleNames = sampleName, rmZeroCov = rmZeroCov)
}

read.bsmooth <- function(dirs, sampleNames = NULL, seqnames = NULL, returnRaw = FALSE,
                         qualityCutoff = 20, rmZeroCov = FALSE, verbose = TRUE) {
    dirs <- normalizePath(dirs, mustWork = TRUE)
    if(!(all(file.info(dirs)$isdir)))
        stop("argument 'dirs' has to be directories")
    if(anyDuplicated(dirs))
        stop("duplicate entries in 'dirs'")
    if(is.null(sampleNames) && !anyDuplicated(basename(dirs)))
        sampleNames <- basename(dirs)
    if(is.null(sampleNames))
        sampleNames <- dirs
    if(!is.null(sampleNames) && (length(sampleNames) != length(dirs) || anyDuplicated(sampleNames)))
        stop("argument 'sampleNames' (if not NULL) has to have the same length as argument 'dirs', without duplicate entries")
    idxes <- seq_along(dirs)
    names(idxes) <- sampleNames
    allOut <- lapply(idxes, function(ii) {
        if(verbose) cat(sprintf("[read.bsmooth] Reading dir '%s' ... ", dirs[ii]))
        ptime1 <- proc.time()
        if(returnRaw) {
            out <- read.bsmoothDirRaw(dir = dirs[ii], seqnames = seqnames, keepCycle = TRUE,
                                      keepFilt = TRUE, header = TRUE, verbose = FALSE)
        } else {
            raw <- read.bsmoothDirRaw(dir = dirs[ii], seqnames = seqnames, keepCycle = FALSE,
                                      keepFilt = FALSE, header = TRUE, verbose = FALSE)
            out <- sampleRawToBSseq(raw, qualityCutoff = qualityCutoff, rmZeroCov, sampleName = sampleNames[ii])
        }
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) cat(sprintf("done in %.1f secs\n", stime))
        out
    })
    if(!returnRaw) {
        if(verbose) cat(sprintf("[read.bsmooth] Joining samples ... "))
        ptime1 <- proc.time()
        allOut <- combineList(allOut)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) cat(sprintf("done in %.1f secs\n", stime))
    }
    allOut
}

parsingPipeline <- function(dirs, qualityCutoff = 20, outDir, seqnames = NULL,
                            subdir = "ev_bt2_cpg_tab", timing = FALSE) {
    if(!all(file.exists(dirs)))
        stop("not all directories in 'dirs' exists.")
    cat("[parsingPipeline] Parsing all files.\n")
    oneDir <- function(dir) {
        cat("[parsingPipeline]  dir", basename(dir), ": ")
        base <- basename(dir)
        cat("parsing, ")
        ptime1 <- proc.time()
        raw <- read.bsmoothDirRaw(file.path(dir, subdir),
                                  keepCycle = TRUE, keepFilt = TRUE,
                                  verbose = FALSE)
        assign(paste0(base, ".raw"), raw)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(timing) {
            cat(sprintf("\ndone in %.1f secs\n", stime))
            print(gc())
        }
        cat("saving, ")
        save(list = paste0(base, ".raw"),
             file = file.path(outDir, paste0(base, ".raw.rda")))
        cat("converting, ")
        ptime1 <- proc.time()
        bsseq <- sampleRawToBSseq(raw, qualityCutoff = qualityCutoff,
                                  sampleName = base)
        seqlevels(bsseq)[seqlevels(bsseq) == "chrgi|9626243|ref|NC_001416.1|"] <- "chrLambda"
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(timing) {
            cat(sprintf("\ndone in %.1f secs\n", stime))
            print(gc())
        }
        cat("ordering, ")
        ptime1 <- proc.time()
        bsseq <- chrSelectBSseq(bsseq, order = TRUE,
                                seqnames = seqnames)
        assign(paste0(base, ".bsseq"), bsseq)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(timing) {
            cat(sprintf("\nIn %.1f secs\n", stime))
            print(gc())
        }
        cat("saving, ")
        save(list = paste0(base, ".bsseq"),
             file = file.path(outDir, paste0(base, ".bsseq.rda")))
        cat("done\n")
        NULL
    }
    for(dir in dirs)
        oneDir(dir)
    invisible(dirs)
}
