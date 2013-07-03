setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("BSseq", contains = "hasGRanges", 
         representation(M = "matrix",
                        Cov = "matrix",
                        coef = "matrixOrNULL",
                        se.coef = "matrixOrNULL",
                        trans = "function",
                        parameters = "list",
                        phenoData = "AnnotatedDataFrame"))


setClass("BSseqTstat", contains = "hasGRanges", 
         representation(stats = "matrix",
                        parameters = "list")
         )
         

##
## Simple accessor and replacement functions
##

setValidity("BSseq", function(object) {
    msg <- NULL
    nCpGs <- length(object@gr)
    nSamples <- ncol(object@M)
    if(nCpGs != nrow(object@M) ||
       nCpGs != nrow(object@Cov))
        msg <- c(msg, "length of the 'gr' slot is not equal to the nrows the 'M' or 'Cov' slots")
    if(nSamples != ncol(object@Cov))
        msg <- c(msg, "ncols of the 'M' slot is not equal to the ncols of the 'Cov' slot")
    if(nSamples != nrow(object@phenoData))
        msg <- c(msg, "the 'phenoData' slot does not match the dimensions of the 'M' slot")
    if(any(object@M < 0))
        msg <- c(msg, "the 'M' slot has negative entries")
    if(any(object@Cov < 0))
        msg <- c(msg, "the 'Cov' slot has negative entries")
    if(any(object@M > object@Cov))
        msg <- c(msg, "the 'M' slot has at least one entry bigger than the 'Cov' slot")
    if(!is.null(object@coef) && (nCpGs != nrow(object@coef) ||
                                 nSamples != ncol(object@coef)))
        msg <- c(msg, "nrows and/or ncols of the 'coef' slot does not have the same dimensions as the 'M' slot")
    if(!is.null(object@se.coef) && (nCpGs != nrow(object@se.coef) ||
                                    nSamples != ncol(object@se.coef)))
        msg <- c(msg, "nrows and/or ncols of the 'se.coef' slot does not have the same dimensions as the 'M' slot")
    if(!is.null(rownames(object@M)) ||
       !is.null(rownames(object@Cov)) ||
       !is.null(rownames(object@coef)) ||
       !is.null(rownames(object@se.coef))) warning("unnecessary rownames in object")
    ## FIXME: check samplenames
    if(any(colnames(object@M) != rownames(pData(object))) ||
       any(colnames(object@M) != rownames(object@Cov)))
        msg <- c(msg, "sample names are messed up: colnames of the M slot has to equal colnames of the Cov slot which has to equal rownames of the phenoData slot")
    if(is.null(msg)) TRUE else msg
})



setValidity("BSseqTstat", function(object) {
    msg <- NULL
    if(length(object@gr) != nrow(object@stats))
        msg <- c(msg, "length of 'gr' is different from the number of rows of 'stats'")
    if(is.null(msg)) TRUE else msg
})


setMethod("show", signature(object = "BSseq"),
          function(object) {
              cat("An object of type 'BSseq' with\n")
              cat(" ", nrow(object), "methylation loci\n")
              cat(" ", ncol(object), "samples\n")
              if(hasBeenSmoothed(object)) {
                  cat("has been smoothed with\n")
                  cat(" ", object@parameters$smoothText, "\n")
              } else {
                  cat("has not been smoothed\n")
              }
          })

setMethod("show", signature(object = "BSseqTstat"),
          function(object) {
              cat("An object of type 'BSseqTstat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
              cat("with parameters\n")
              cat(" ", object@parameters$tstatText, "\n")
          })


##
## eSet stuff
##

setMethod("dim", "BSseq", function(x) {
    dim(x@M)
})
setMethod("nrow", "BSseq", function(x) {
    nrow(x@M)
})
setMethod("ncol", "BSseq", function(x) {
    ncol(x@M)
})

setMethod("phenoData", "BSseq", function(object) {
    object@phenoData
})
setReplaceMethod("phenoData", signature = signature(
                              object = "BSseq",
                              value = "AnnotatedDataFrame"),
                 function(object, value) {
                     object@phenoData <- value
                     object
                 })

setMethod("pData", "BSseq", function(object) {
    pData(phenoData(object))
})
setReplaceMethod("pData", signature = signature(
                              object = "BSseq",
                              value = "data.frame"),
                 function(object, value) {
                     pd <- object@phenoData
                     pData(pd) <- value
                     object@phenoData <- pd
                     object
                 })

setMethod("sampleNames", "BSseq", function(object) {
    sampleNames(phenoData(object))
})

setReplaceMethod("sampleNames",
                 signature(object = "BSseq", value = "ANY"),
                 function(object, value) {
                     sampleNames(phenoData(object)) <- value
                     colnames(object@M) <- value
                     colnames(object@Cov) <- value
                     if(!is.null(object@coef))
                         colnames(object@coef) <- value
                     if(!is.null(object@se.coef))
                         colnames(object@se.coef) <- value
                     object
                 })

hasBeenSmoothed <- function(BSseq) {
    !is.null(BSseq@coef)
}

##
## Subsetting and combining
##

setMethod("[", "BSseq", function(x, i, j, ...) {
    if(missing(drop))
        drop <- FALSE
    if(missing(i) && missing(j))
        stop("need [i,j] for subsetting")
    if(!missing(j))
        x@phenoData <- phenoData(x)[j,, ..., drop = drop]
    else
        x@phenoData <- phenoData(x)
    if(missing(i)) {
        x@M <- getBSseq(x, "M")[, j, drop = FALSE]
        x@Cov <- getBSseq(x, "Cov")[, j, drop = FALSE]
        x@coef <- getBSseq(x, "coef")[, j, drop = FALSE]
        x@se.coef <- getBSseq(x, "se.coef")[, j, drop = FALSE]
        x@gr <- granges(x)
    } else {
        x@M <- getBSseq(x, "M")[i, j, drop = FALSE]
        x@Cov <- getBSseq(x, "Cov")[i, j, drop = FALSE]
        x@coef <- getBSseq(x, "coef")[i, j, drop = FALSE]
        x@se.coef <- getBSseq(x, "se.coef")[i, j, drop = FALSE]
        x@gr <- granges(x)[i]
    }
    x
})

setMethod("[", "BSseqTstat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x@stats <- x@stats[i,, drop = FALSE]
    x
})

setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    ## All of this assumes that we are place x and y "next" to each other,
    ##  ie. we are not combining the same set of samples sequenced at different times
    if (class(x) != class(y))
        stop(paste("objects must be the same class, but are ",
                   class(x), ", ", class(y), sep=""))
    if(hasBeenSmoothed(x) && hasBeenSmoothed(y) && !all.equal(x@trans, y@trans))
        stop("'x' and 'y' need to be smoothed on the same scale")
    phenoData <- combine(phenoData(x), phenoData(y))
    if(identical(granges(x), granges(y))) {
        gr <- granges(x)
        M <- cbind(getBSseq(x, "M"), getBSseq(y, "M"))
        Cov <- cbind(getBSseq(x, "Cov"), getBSseq(y, "Cov"))
        if(!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            coef <- cbind(getBSseq(x, "coef"), getBSseq(y, "coef"))
            se.coef <- cbind(getBSseq(x, "se.coef"), getBSseq(y, "se.coef"))
            trans <- getBSseq(x, "trans")
        }
    } else {
        gr <- reduce(c(granges(x), granges(y)), min.gapwidth = 0L)
        mm.x <- as.matrix(findOverlaps(gr, granges(x)))
        mm.y <- as.matrix(findOverlaps(gr, granges(y)))
        sampleNames <- c(sampleNames(x), sampleNames(y))
        ## FIXME: there is no check that the two sampleNames are disjoint.
        M <- Cov <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
        colnames(M) <- colnames(Cov) <- sampleNames
        M[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "M")[mm.x[,2],]
        M[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "M")[mm.y[,2],]
        Cov[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "Cov")[mm.x[,2],]
        Cov[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "Cov")[mm.y[,2],]
        if(!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            trans <- x@trans
            coef <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
            colnames(coef) <- sampleNames(phenoData)
            if(hasBeenSmoothed(x))
                coef[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "coef")[mm.x[,2],]
            if(hasBeenSmoothed(y))
                coef[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "coef")[mm.y[,2],]
            if(is.null(getBSseq(x, "se.coef")) && is.null(getBSseq(x, "se.coef")))
                se.coef <- NULL
            else {
                se.coef <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
                colnames(se.coef) <- sampleNames(phenoData)
                if(!is.null(getBSseq(x, "se.coef")))
                    se.coef[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "se.coef")[mm.x[,2],]
                if(!is.null(getBSseq(y, "se.coef")))
                    se.coef[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "se.coef")[mm.y[,2],]
            }
        }
    }
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          phenoData = phenoData, trans = trans, rmZeroCov = FALSE)
})

combineList <- function(x, ...) {
    if(class(x) == "BSseq")
        x <- list(x, ...)
    stopifnot(all(sapply(x, class) == "BSseq"))
    gr <- getBSseq(x[[1]], "gr")
    trans <- getBSseq(x[[1]], "trans")
    sameTrans <- sapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "trans"))
    })
    if(!all(sameTrans))
        stop("all elements of '...' in combineList needs to have the same trans")
    sameGr <- sapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "gr"))
    })
    if(all(sameGr)) {
        M <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "M")))
        Cov <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "Cov")))
    } else {
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))), min.gapwidth = 0L))
        sampleNames <- do.call(c, lapply(x, sampleNames))
        M <- matrix(0, ncol = length(sampleNames), nrow = length(gr))
        colnames(M) <- sampleNames
        Cov <- M
        idxes <- split(seq(along = sampleNames), rep(seq(along = x), times = sapply(x, ncol)))
        for(ii in seq(along = idxes)) {
            ov <- findOverlaps(gr, granges(x[[ii]]))
            idx <- idxes[[ii]]
            M[queryHits(ov), idx] <- getBSseq(x[[ii]], "M")[subjectHits(ov),]
            Cov[queryHits(ov), idx] <- getBSseq(x[[ii]], "Cov")[subjectHits(ov),]
        }
    }
    if(any(!sapply(x, hasBeenSmoothed)) || !(all(sameGr))) {
        coef <- NULL
        se.coef <- NULL
        trans <- NULL
    } else {
        coef <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "coef")))
        se.coef <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "se.coef")))
    }
    phenoData <- Reduce(combine, lapply(x, phenoData))
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          phenoData = phenoData, trans = trans, rmZeroCov = FALSE)
}

getBSseq <- function(BSseq, type = c("Cov", "M", "gr", "coef", "se.coef", "trans", "parameters")) {
    type <- match.arg(type)
    switch(type,
           "Cov" = {
               return(BSseq@Cov)
           },
           "M" = {
               return(BSseq@M)
           },
           "gr" = {
               return(BSseq@gr)
           },
           "coef" = {
               return(BSseq@coef)
           },
           "se.coef" = {
               return(BSseq@se.coef)
           },
           "trans" = {
               return(BSseq@trans)
           },
           "parameters" = {
               return(BSseq@parameters)
           })
}



BSseq <- function(M = NULL, Cov = NULL, coef = NULL, se.coef = NULL,
                  trans = NULL, parameters = NULL, phenoData = NULL,
                  gr = NULL, pos = NULL, chr = NULL, sampleNames = NULL,
                  rmZeroCov = FALSE) {
    if(is.null(gr)) {
        if(is.null(pos) || is.null(chr))
            stop("Need pos and chr")
        gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
    }
    if(!is(gr, "GRanges"))
        stop("'gr' needs to be a GRanges")
    if(any(width(gr)) != 1)
        stop("'gr' needs to have widths of 1")
    if(is.null(M) || is.null(Cov))
        stop("Need M and Cov")
    if(!is.matrix(M))
        stop("'M' needs to be a matrix")
    if(!is.matrix(Cov))
        stop("'Cov' needs to be a matrix")
    if(length(gr) != nrow(M) ||
       length(gr) != nrow(Cov) ||
       ncol(Cov) != ncol(M))
        stop("'gr', 'M' and 'Cov' need to have similar dimensions")
    if(!is.null(rownames(M)))
        rownames(M) <- NULL
    if(!is.null(rownames(Cov)))
        rownames(Cov) <- NULL
    if(!is.null(names(gr)))
        names(gr) <- NULL
    ## deal with sampleNames
    if(is.null(sampleNames) && !is.null(phenoData) && !is.null(sampleNames(phenoData)))
        sampleNames <- sampleNames(phenoData)
    if(is.null(sampleNames) && !is.null(colnames(M)))
        sampleNames <- colnames(M)
    if(is.null(sampleNames) && !is.null(colnames(Cov)))
        sampleNames <- colnames(Cov)
    if(is.null(sampleNames))
        sampleNames <- paste("V", 1:ncol(M), sep = "")
    if(length(unique(sampleNames)) != ncol(M))
        stop("sampleNames need to be unique and of the right length.")
    ## check that 0 <= M <= Cov and remove positions with Cov = 0
    if(any(M < 0) || any(M > Cov) || any(is.na(M)) || any(is.na(Cov)) ||
       any(is.infinite(Cov)))
        stop("'M' and 'Cov' may not contain NA or infinite values and 0 <= M <= Cov")
    if(rmZeroCov) {
        wh <- which(rowSums(Cov) == 0)
        if(length(wh) > 0) {
            gr <- gr[-wh]
            M <- M[-wh,,drop = FALSE]
            Cov <- Cov[-wh,,drop = FALSE]
        }
    }
    grR <- reduce(gr, min.gapwidth = 0L)
    if(!identical(grR, gr)) {
        ## Now we either need to re-order or collapse or both
        mm <- as.matrix(findOverlaps(grR, gr))
        mm <- mm[order(mm[,1]),]
        if(length(grR) == length(gr)) {
            ## only re-ordering is necessary 
            gr <- grR
            M <- M[mm[,2],,drop = FALSE]
            Cov <- Cov[mm[,2],,drop = FALSE]
            if(!is.null(coef))
                coef <- coef[mm[,2],,drop = FALSE]
            if(!is.null(se.coef))
                se.coef <- se.coef[mm[,2],, drop = FALSE]
        } else {
            warning("multiple positions, collapsing BSseq object\n")
            if(!is.null(coef) || !is.null(se.coef))
                stop("Cannot collapse when 'coef' or 'se.coef' are present")
            gr <- grR
            sp <- split(mm[,2], mm[,1])[as.character(1:length(grR))]
            names(sp) <- NULL
            M <- do.call(rbind, lapply(sp, function(ii) colSums(M[ii,, drop = FALSE])))
            Cov <- do.call(rbind, lapply(sp, function(ii) colSums(Cov[ii,, drop = FALSE])))
        }
    }
    if(is.null(colnames(M)) || any(sampleNames != colnames(M)))
        colnames(M) <- sampleNames
    if(is.null(colnames(Cov)) || any(sampleNames != colnames(Cov)))
        colnames(Cov) <- sampleNames
    if(is.null(phenoData))
        phenoData <- annotatedDataFrameFrom(M, byrow = FALSE)
    BSseq <- new("BSseq", gr = gr, M = M, Cov = Cov, phenoData = phenoData)
    if(!is.null(coef)) {
        if(!is.matrix(coef) ||
           nrow(coef) != nrow(BSseq) ||
           ncol(coef) != ncol(BSseq))
            stop("'coef' does not have the right dimensions")
        if(is.null(colnames(coef)) || any(sampleNames != colnames(coef)))
            colnames(coef) <- sampleNames
        if(!is.null(rownames(coef)))
            rownames(coef) <- NULL
        BSseq@coef <- coef
    }
    if(!is.null(se.coef)) {
        if(!is.matrix(se.coef) ||
           nrow(se.coef) != nrow(BSseq) ||
           ncol(se.coef) != ncol(BSseq))
            stop("'se.coef' does not have the right dimensions")
        if(is.null(colnames(se.coef)) || any(sampleNames != colnames(se.coef)))
            colnames(se.coef) <- sampleNames
        if(!is.null(rownames(se.coef)))
            rownames(se.coef) <- NULL
        BSseq@se.coef <- se.coef
    }
    if(is.function(trans))
        BSseq@trans <- trans
    if(is.list(parameters))
        BSseq@parameters <- parameters
    BSseq
}


BSseqTstat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqTstat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}



                  

## getD <- function(data, sample1, sample2, type = c("raw", "fit"),
##                  addPositions = FALSE, alpha = 0.95) {
##     d.conf <- function(p1, n1, p2, n2) {
##         ## Method 10 from Newcombe (Stat in Med, 1998)
##         getRoots <- function(p, n) {
##             mat <- cbind(p^2, -(2*p + z^2/n), 1+z^2/n)
##             roots <- t(Re(apply(mat, 1, polyroot)))
##             roots
##         }
##         roots1 <- roots2 <- matrix(NA_real_, nrow = length(p1), ncol = 2)
##         idx <- !is.na(p1)
##         roots1[idx,] <- getRoots(p1[idx], n1[idx])
##         idx <- !is.na(p2)
##         roots2[idx,] <- getRoots(p2[idx], n2[idx])
##         d <- p1 - p2
##         suppressWarnings({
##             lower <- d - z * sqrt(roots1[,1]*(1-roots1[,1])/n1 + roots2[,2]*(1-roots2[,2])/n2)
##             upper <- d + z * sqrt(roots1[,2]*(1-roots1[,2])/n1 + roots2[,1]*(1-roots2[,1])/n2)
##         })
##         return(data.frame(d = d, lower = lower, upper = upper))
##     }
##     type <- match.arg(type)
##     z <- abs(qnorm((1-alpha)/2, mean = 0, sd = 1))
##     switch(type,
##            raw = {
##                conf <- d.conf(p1 = data$M[, sample1] / data$Cov[, sample1],
##                               n1 = data$Cov[, sample1],
##                               p2 = data$M[, sample2] / data$Cov[, sample2],
##                               n2 = data$Cov[, sample2])
##                out <- conf
##            },
##            fit = {
##                p1 <- data$trans(data$coef[, sample1])
##                p2 <- data$trans(data$coef[, sample2])
##                d <- p1 - p2
##                ## Based on the delta method
##                se.d <- sqrt((data$se.coef[, sample1] * p1 * (1-p1))^2 +
##                             (data$se.coef[, sample2] * p2 * (1-p2))^2)
               
##                lower <- d - z * se.d
##                upper <- d + z * se.d
##                out <- data.frame(d = d, lower = lower, upper = upper)
##            })
##     if(addPositions) {
##         out$pos <- start(data$gr)
##         out$chr <- seqnames(data$gr)
##     }
##     out
## }
