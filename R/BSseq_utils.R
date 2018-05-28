chrSelectBSseq <- function(BSseq, seqnames = NULL, order = FALSE) {
    seqlevels(BSseq, pruning.mode = "coarse") <- seqnames
    if (order) BSseq <- orderBSseq(BSseq, seqOrder = seqnames)
    BSseq
}

orderBSseq <- function(BSseq, seqOrder = NULL) {
    if (!is.null(seqOrder)) {
        seqlevels(BSseq, pruning.mode = "coarse") <- seqOrder
    }
    BSseq[order(granges(BSseq))]
}

# TODO: getMeth() realises the result in memory iff regions is not NULL;
#       discuss with Kasper
# TODO: Whether or not colnames are added to returned value depends on whether
#       regions is non-NULL; discuss with Kasper
# TODO: Add parallel support
# TODO: Document withDimnames
getMeth <- function(BSseq, regions = NULL, type = c("smooth", "raw"),
                    what = c("perBase", "perRegion"), confint = FALSE,
                    alpha = 0.95, withDimnames = TRUE) {
    p.conf <- function(p, n, alpha) {
        z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
        upper <- (p + z ^ 2 / (2 * n) +
                      z * sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
            (1 + z ^ 2 / n)
        lower <- (p + z ^ 2 / (2 * n) -
                      z * sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
            (1 + z ^ 2 / n)
        return(list(meth = p, lower = lower, upper = upper))
    }

    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    if (type == "smooth" & !hasBeenSmoothed(BSseq)) {
        stop("'type=smooth' requires the object to have been smoothed.")
    }
    what <- match.arg(what)
    if (what == "perRegion" & is.null(regions)) {
        stop("'what=perRegion' but no 'regions' supplied")
    }
    z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
    if (is.null(regions) && type == "smooth") {
        coef <- getBSseq(BSseq, "coef", withDimnames)
        meth <- getBSseq(BSseq, "trans", withDimnames)(coef)
        if (confint) {
            upper <- meth + z * getBSseq(BSseq, "se.coef", withDimnames)
            lower <- meth - z * getBSseq(BSseq, "se.coef", withDimnames)
            return(list(meth = meth, lower = lower, upper = upper))
        } else {
            return(meth)
        }
    }
    if (is.null(regions) && type == "raw") {
        meth <- getBSseq(BSseq, "M", withDimnames) /
            getBSseq(BSseq, "Cov", withDimnames)
        if (confint) {
            return(p.conf(meth, getBSseq(BSseq, "Cov", withDimnames), alpha))
        } else {
            return(meth)
        }
    }

    ## At this point, regions have been specified
    if (class(regions) == "data.frame") {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))
    if (confint) {
        stop("'confint = TRUE' is not supported by 'getMeth' when regions is given")
    }
    grBSseq <- granges(BSseq)
    ov <- findOverlaps(grBSseq, regions)
    # NOTE: This realises a large object in memory (`meth`) - could do it in
    #       chunks if what = perRegion
    if (type == "smooth") {
        meth <- as.matrix(
            getBSseq(BSseq, "trans", withDimnames)(
                getBSseq(BSseq, "coef", withDimnames))[
                    queryHits(ov), , drop = FALSE])
    } else if (type == "raw") {
        meth <- as.matrix(
            (getBSseq(BSseq, "M", withDimnames) /
                 getBSseq(BSseq, "Cov", withDimnames))[
                     queryHits(ov), , drop = FALSE])
    }
    out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
    if (what == "perBase") {
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the results end up in the wrong order wrt to regions?
        outList <- vector("list", length(regions))
        outList[as.integer(names(out))] <- out
        return(outList)
    } else if (what == "perRegion") {
        out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the rows end up in the wrong order?
        outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
        if (withDimnames) colnames(outMatrix) <- sampleNames(BSseq)
        outMatrix[as.integer(rownames(out)), ] <- out
        outMatrix
    }
}

# TODO: getCoverage() realises the result in memory iff regions is not NULL;
#       discuss with Kasper
# TODO: Whether or not colnames are added to returned value depends on whether
#       regions is non-NULL; discuss with Kasper
# TODO: Document withDimnames
getCoverage <- function(BSseq, regions = NULL, type = c("Cov", "M"),
                        what = c("perBase", "perRegionAverage",
                                 "perRegionTotal"),
                        withDimnames = TRUE) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    what <- match.arg(what)
    if (is.null(regions)) {
        if (what == "perBase") {
            return(getBSseq(BSseq, type, withDimnames))
        }
        if (what == "perRegionTotal") {
            return(colSums2(getBSseq(BSseq, type, withDimnames)))
        }
        if (what == "perRegionAverage") {
            return(colMeans2(getBSseq(BSseq, type, withDimnames)))
        }
    }
    if (class(regions) == "data.frame") {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))
    grBSseq <- granges(BSseq)
    ov <- findOverlaps(grBSseq, regions)
    coverage <- getBSseq(BSseq, type, withDimnames)[
        queryHits(ov), , drop = FALSE]
    out <- lapply(split(coverage, subjectHits(ov)), matrix,
                  ncol = ncol(coverage))

    if (what == "perBase") {
        # TODO: Don't really understand the logic of the remaining code; how
        #       could the results end up in the wrong order wrt to regions?
        outList <- vector("list", length(regions))
        outList[as.integer(names(out))] <- out
        return(outList)
    } else if (what == "perRegionAverage") {
        out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
    } else if (what == "perRegionTotal") {
        out <- do.call(rbind, lapply(out, colSums2, na.rm = TRUE))
    }
    # TODO: Don't really understand the logic of the remaining code; how
    #       could the rows end up in the wrong order?
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    if (withDimnames) colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)), ] <- out
    outMatrix
}
