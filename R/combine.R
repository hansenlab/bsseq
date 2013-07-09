setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    ## All of this assumes that we are place x and y "next" to each other,
    ##  ie. we are not combining the same set of samples sequenced at different times
    if (class(x) != class(y))
        stop(paste("objects must be the same class, but are ",
                   class(x), ", ", class(y), sep=""))
    if(hasBeenSmoothed(x) && hasBeenSmoothed(y) && !all.equal(x@trans, y@trans))
        stop("'x' and 'y' need to be smoothed on the same scale")
    pData <- combine(as(pData(x), "data.frame"), as(pData(y), "data.frame"))
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
            colnames(coef) <- rownames(pData)
            if(hasBeenSmoothed(x))
                coef[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "coef")[mm.x[,2],]
            if(hasBeenSmoothed(y))
                coef[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "coef")[mm.y[,2],]
            if(is.null(getBSseq(x, "se.coef")) && is.null(getBSseq(x, "se.coef")))
                se.coef <- NULL
            else {
                se.coef <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
                colnames(se.coef) <- sampleNames(pData)
                if(!is.null(getBSseq(x, "se.coef")))
                    se.coef[mm.x[,1], 1:ncol(x)] <- getBSseq(x, "se.coef")[mm.x[,2],]
                if(!is.null(getBSseq(y, "se.coef")))
                    se.coef[mm.y[,1], ncol(x) + 1:ncol(y)] <- getBSseq(y, "se.coef")[mm.y[,2],]
            }
        }
    }
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
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
    pData <- as(Reduce(combine, lapply(x, function(xx) as.data.frame(pData(xx)))), "DataFrame")
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
}

