collapseBSseq <- function(BSseq, columns) {
    ## columns is a vector of new names, names(columns) is sampleNames or empty
    stopifnot(is.character(columns))
    if(is.null(names(columns)) && length(columns) != ncol(BSseq))
        stop("if `columns' does not have names, it needs to be of the same length as `BSseq` has columns (samples)")
    if(!is.null(names(columns)) && !all(names(columns) %in% sampleNames(BSseq)))
        stop("if `columns` has names, they need to be sampleNames(BSseq)")
    if(is.null(names(columns)))
        columns.idx <- 1:ncol(BSseq)
    else
        columns.idx <- match(names(columns), sampleNames(BSseq))
    sp <- split(columns.idx, columns)
    M <- do.call(cbind, lapply(sp, function(ss) {
        rowSums(getBSseq(BSseq, "M"), cols = ss)
    }))
    Cov <- do.call(cbind, lapply(sp, function(ss) {
        rowSums(getBSseq(BSseq, "Cov"), cols = ss)
    }))
    BSseq(gr = getBSseq(BSseq, "gr"), M = M, Cov = Cov, sampleNames = names(sp))
}

chrSelectBSseq <- function(BSseq, seqnames = NULL, order = FALSE) {
    seqlevels(BSseq, force = TRUE) <- seqnames
    if(order)
        BSseq <- orderBSseq(BSseq, seqOrder = seqnames)
    BSseq
}


orderBSseq <- function(BSseq, seqOrder = NULL) {
    if(!is.null(seqOrder))
        seqlevels(BSseq, force = TRUE) <- seqOrder
    BSseq[order(granges(BSseq))]
}


getMeth <- function(BSseq, regions = NULL, type = c("smooth", "raw"),
                    what = c("perBase", "perRegion"), confint = FALSE, alpha = 0.95) {
    p.conf <- function(p, n, alpha) {
        z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
        upper <- (p + z^2/(2*n) + z*sqrt(  (p*(1-p) + z^2/(4*n)) / n)) /
            (1+z^2/n)
        lower <- (p + z^2/(2*n) - z*sqrt(  (p*(1-p) + z^2/(4*n)) / n)) /
            (1+z^2/n)
        return(list(meth = p, lower = lower, upper = upper))
    }
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    if(type == "smooth" & !hasBeenSmoothed(BSseq))
        stop("'type=smooth' requires the object to have been smoothed.")
    what <- match.arg(what)
    z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
    if(is.null(regions) && type == "smooth") {
        meth <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef"))
        if(confint) {
            upper <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef") +
                                      z * getBSseq(BSseq, type = "se.coef"))
            lower <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef") -
                                      z * getBSseq(BSseq, type = "se.coef"))
            return(list(meth = meth, lower = lower, upper = upper))
        } else {
            return(meth)
        }
    }
    if(is.null(regions) && type == "raw") {
        meth <- getBSseq(BSseq, type = "M") / getBSseq(BSseq, type = "Cov")
        if(confint) {
            return(p.conf(meth, n = getBSseq(BSseq, type = "Cov"), alpha = alpha))
        } else {
            return(meth)
        }
    }
    ## At this point, regions have been specified
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(is(regions, "GenomicRanges"))
    if(confint) stop("'confint = TRUE' is not supported by 'getMeth' when regions is given")
    grBSseq <- granges(BSseq)
    mm <- as.matrix(findOverlaps(regions, grBSseq))
    mmsplit <- split(mm[,2], mm[,1])
    if(what == "perBase") {
        if(type == "smooth") {
            out <- lapply(mmsplit, function(xx) {
                getBSseq(BSseq, "trans")(getBSseq(BSseq, "coef")[xx,,drop = FALSE])
            })
        }
        if(type == "raw") {
            out <- lapply(mmsplit, function(xx) {
                getBSseq(BSseq, "M")[xx,,drop = FALSE] / getBSseq(BSseq, "Cov")[xx,,drop = FALSE]
            })
        }
        outList <- vector("list", length(regions))
        outList[as.integer(names(mmsplit))] <- out
        return(outList)
    }
    if(what == "perRegion") {
        if(type == "smooth") {
            out <- lapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "trans")(getBSseq(BSseq, "coef")[xx,,drop = FALSE]), na.rm = TRUE)
            })
        }
        if(type == "raw") {
            out <- lapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "M")[xx,,drop = FALSE] / getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        out <- do.call(rbind, out)
        outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
        colnames(outMatrix) <- sampleNames(BSseq)
        outMatrix[as.integer(rownames(out)),] <- out
        return(outMatrix)
    }
}


getCoverage <- function(BSseq, regions = NULL, type = c("Cov", "M"),
                    what = c("perBase", "perRegionAverage", "perRegionTotal")) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    what <- match.arg(what)
    if(is.null(regions)) {
        if(what == "perBase")
            return(getBSseq(BSseq, type = type))
        if(what == "perRegionTotal")
            return(colSums(getBSseq(BSseq, type = type)))
        if(what == "perRegionAverage")
            return(colMeans(getBSseq(BSseq, type = type)))
    }
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(is(regions, "GenomicRanges"))
    grBSseq <- granges(BSseq)
    mm <- as.matrix(findOverlaps(regions, grBSseq))
    mmsplit <- split(mm[,2], mm[,1])
    if(what == "perBase") {
        if(type == "Cov") {
            out <- lapply(mmsplit, function(xx) {
                getBSseq(BSseq, "Cov")[xx,,drop = FALSE]
            })
        }
        if(type == "M") {
            out <- lapply(mmsplit, function(xx) {
                getBSseq(BSseq, "M")[xx,,drop = FALSE]
            })
        }
        outList <- vector("list", length(regions))
        outList[as.integer(names(mmsplit))] <- out
        return(outList)
    }
    if(what == "perRegionAverage") {
        if(type == "Cov") {
            out <- lapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        if(type == "M") {
            out <- lapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "M")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
    }
    if(what == "perRegionTotal") {
        if(type == "Cov") {
            out <- lapply(mmsplit, function(xx) {
                colSums(getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        if(type == "M") {
            out <- lapply(mmsplit, function(xx) {
                colSums(getBSseq(BSseq, "M")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
    }
    out <- do.call(rbind, out)
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)),] <- out
    return(outMatrix)
}


