plotAnnoTrack <- function(gr, annoTrack, cex) {
    ## check may need to be modified
    if(!all(sapply(annoTrack, function(xx) is(xx, "GRanges"))))
        stop("all elements in 'annoTrack' needs to be 'GRanges'")
    plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
         ylim = c(0.5, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), xlab = "", ylab = "")
    lapply(seq(along = annoTrack), function(ii) {
        jj <- length(annoTrack) + 1- ii
        ir <- subsetByOverlaps(annoTrack[[ii]], gr)
        if(length(ir) > 0)
            rect(start(ir)-0.5, jj - 0.15, end(ir), jj + 0.15, col = "grey60", border = NA)
        mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1,
              cex = cex, adj = 0.3)
        })
}

# Based on plotAnnoTrack() and
# https://jhu-genomics.slack.com/archives/hansen_gtex/p1461760583000142
# TODO: Use standard Bioconductor object for storing gene information rather
#       than this ad hoc data.frame (e.g., based on TxDb.* packages)
plotGeneTrack <- function(gr, geneTrack, cex) {
    geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
    ol <- findOverlaps(geneTrack_gr, gr)
    genes <- geneTrack[queryHits(ol), ]
    plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
         ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)),
         xlab = "", ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
    if (nrow(genes) > 0) {
        for (g in 1:nrow(genes)) {
            geneind2 = which(geneTrack$gene_name == genes$gene_name[g])
            geneind2 = geneind2[which(geneTrack$isoforms[geneind2] == 1)]
            direction = unique(geneTrack$strand[geneind2])
            ES = geneTrack$start[geneind2]
            EE = geneTrack$end[geneind2]
            Exons = cbind(ES, EE)
            if (direction == "+") {
                lines(x = c(min(ES), max(EE)),
                      y = c(0.65, 0.65))
                apply(Exons, 1, function(x) {
                    polygon(c(x[1], x[2], x[2], x[1]),
                            c(0.45, 0.45, 0.85, 0.85), col = "darkgrey")
                    })
                text((max(start(gr), min(ES)) +
                          min(end(gr), max(EE))) / 2, 1.2,
                     genes$gene_name[g], cex = cex)
            } else {
                lines(x = c(min(ES), max(EE)),
                      y = c(-0.65, -0.65))
                apply(Exons, 1, function(x)
                    polygon(c(x[1], x[2], x[2], x[1]),
                            c(-0.45, -0.45, -0.85, -0.85),
                            col = "darkgrey"))
                text((max(start(gr), min(ES)
                ) + min(end(gr), max(EE)
                )) / 2, -1.2, genes$gene_name[g], cex = cex)
            }

        }
    }
}

plotManyRegions <- function(BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL,
                            annoTrack = NULL, cex.anno = 1,
                            geneTrack = NULL, cex.gene = 1.5,
                            col = NULL, lty = NULL, lwd = NULL,
                            BSseqStat = NULL, stat = "tstat.corrected", stat.col = "black",
                            stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                            mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
                            addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE, verbose = TRUE) {
    cat("[plotManyRegions] preprocessing ...")
    if(!is.null(regions)) {
        if(is(regions, "data.frame"))
            gr <- data.frame2GRanges(regions, keepColumns = FALSE)
        else
            gr <- regions
        if(!is(gr, "GRanges"))
            stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
    } else {
        gr <- granges(BSseq)
    }
    gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
    BSseq <- subsetByOverlaps(BSseq, gr)
    if(!is.null(BSseqStat))
        BSseqStat <- subsetByOverlaps(BSseqStat, gr)

    if(length(start(BSseq)) == 0)
        stop("No overlap between BSseq data and regions")
    if(!is.null(main) && length(main) != length(gr))
        main <- rep(main, length = length(gr))
    cat("done\n")
    for(ii in seq(along = gr)) {
        if(verbose) cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n", ii, nrow(regions)))
        plotRegion(BSseq = BSseq, region = regions[ii,], extend = extend,
                   col = col, lty = lty, lwd = lwd, main = main[ii], BSseqStat = BSseqStat,
                   stat = stat, stat.col = stat.col, stat.lwd = stat.lwd,
                   stat.lty = stat.lty, stat.ylim = stat.ylim,
                   addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth,
                   annoTrack = annoTrack, cex.anno = cex.anno,
                   geneTrack = geneTrack, cex.gene = cex.gene,
                   addTicks = addTicks, addPoints = addPoints,
                   pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    }
}

.bsPlotLines <- function(x, y, col, lty, lwd, plotRange) {
    if(sum(!is.na(y)) <= 1)
        return(NULL)
    xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 500)
    yy <- approxfun(x, y)(xx)
    lines(xx, yy, col = col, lty = lty, lwd = lwd)
}

.bsPlotPoints <- function(x, y, z, col, pointsMinCov) {
    points(x[z>pointsMinCov], y[z>pointsMinCov], col = col, pch = 16, cex = 0.5)
}

.bsHighlightRegions <- function(regions, gr, ylim, regionCol, highlightMain) {
    if(is.data.frame(regions))
        regions <- data.frame2GRanges(regions)
    if(highlightMain)
        regions <- c(regions, gr)
    if(is.null(regions)) return(NULL)
    ## regions <- pintersect(region, rep(gr, length(regions)))
    ## regions <- regions[width(regions) == 0]
    regions <- subsetByOverlaps(regions, gr)
    regions <- pintersect(regions, rep(gr, length(regions)))
    if(length(regions) == 0)
        return(NULL)
    rect(xleft = start(regions), xright = end(regions), ybottom = ylim[1],
         ytop = ylim[2], col = regionCol, border = NA)
}

.bsGetCol <- function(object, col, lty, lwd) {
    ## Assumes that object has pData and sampleNames methods
    if(is.null(col)) {
        if("col" %in% names(pData(object)))
            col <- pData(object)[["col"]]
        else
            col <- rep("black", nrow(pData(object)))
    }
    if(length(col) != ncol(object))
        col <- rep(col, length.out = ncol(object))
    if(is.null(names(col)))
        names(col) <- sampleNames(object)

    if(is.null(lty)) {
        if("lty" %in% names(pData(object)))
            lty <- pData(object)[["lty"]]
        else
            lty <- rep(1, ncol(object))
    }
    if(length(lty) != ncol(object))
        lty <- rep(lty, length.out = ncol(object))
    if(is.null(names(lty)))
        names(lty) <- sampleNames(object)

    if(is.null(lwd)) {
        if("lwd" %in% names(pData(object)))
            lwd <- pData(object)[["lwd"]]
        else
            lwd <- rep(1, nrow(pData(object)))
    }
    if(length(lwd) != ncol(object))
        lwd <- rep(lwd, length.out = ncol(object))
    if(is.null(names(lwd)))
        names(lwd) <- sampleNames(object)

    return(list(col = col, lty = lty, lwd = lwd))
}

.bsPlotTitle <- function(gr, extend, main, mainWithWidth) {
    if(is.data.frame(gr))
        gr <- data.frame2GRanges(gr)
    if(length(gr) > 1) {
        warning("plotTitle: gr has more than one element")
        gr <- gr[1]
    }
    plotChr <- as.character(seqnames(gr))
    plotRange <- c(start(gr), end(gr))
    regionCoord <- sprintf("%s: %s - %s", plotChr,
                           format(plotRange[1], big.mark = ",", scientific = FALSE),
                           format(plotRange[2], big.mark = ",", scientific = FALSE))
    if(mainWithWidth) {
        regionWidth <- sprintf("width = %s, extended = %s",
                               format(width(gr), big.mark = ",", scientific = FALSE),
                               format(extend, big.mark = ",", scientific = FALSE))
        regionCoord <- sprintf("%s (%s)", regionCoord, regionWidth)
    }
    if(main != "") {
        main <- sprintf("%s\n%s", main, regionCoord)
    } else {
        main <- regionCoord
    }
    main
}

.bsGetGr <- function(object, region, extend) {
    if(is.null(region)) {
        gr <- GRanges(seqnames = seqnames(object)[1],
                      ranges = IRanges(start = min(start(object)),
                      end = max(start(object))))
    } else {
        if(is(region, "data.frame"))
            gr <- data.frame2GRanges(region, keepColumns = FALSE)
        else
            gr <- region
        if(!is(gr, "GRanges") || length(gr) != 1)
            stop("'region' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
        gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
    }
    gr
}


.plotSmoothData <- function(BSseq, region, extend, addRegions, col, lty, lwd, regionCol,
                            addTicks, addPoints, pointsMinCov, highlightMain) {
    gr <- .bsGetGr(BSseq, region, extend)
    BSseq <- subsetByOverlaps(BSseq, gr)

    ## Extract basic information
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    positions <- start(BSseq)
    smoothPs <- getMeth(BSseq, type = "smooth")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)

    # Realise in memory data that are to be plotted
    if (addPoints) {
        rawPs <- as.array(rawPs)
        coverage <- as.array(coverage)
    }
    smoothPs <- as.array(smoothPs)

    ## get col, lwd, lty
    colEtc <- .bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)

    ## The actual plotting
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = c(0,1), xlim = c(start(gr), end(gr)), xlab = "", ylab = "Methylation")
    axis(side = 2, at = c(0.2, 0.5, 0.8))
    if(addTicks)
        rug(positions)

    .bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0,1),
                        regionCol = regionCol, highlightMain = highlightMain)

    if(addPoints) {
        sapply(1:ncol(BSseq), function(sampIdx) {
            abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
        })
    } # This adds vertical grey lines so we can see where points are plotted

    sapply(1:ncol(BSseq), function(sampIdx) {
        .bsPlotLines(positions, smoothPs[, sampIdx], col = colEtc$col[sampIdx],
                     lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx],
                     plotRange = c(start(gr), end(gr)))
    })

    if(addPoints) {
        sapply(1:ncol(BSseq), function(sampIdx) {
            .bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx],
                          col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov)
        })
    }
}


plotRegion <- function(BSseq, region = NULL, extend = 0, main = "", addRegions = NULL,
                       annoTrack = NULL, cex.anno = 1,
                       geneTrack = NULL, cex.gene = 1.5,
                       col = NULL, lty = NULL, lwd = NULL,
                       BSseqStat = NULL, stat = "tstat.corrected", stat.col = "black",
                       stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                       mainWithWidth = TRUE,
                       regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE,
                       pointsMinCov = 5, highlightMain = FALSE) {

    opar <- par(mar = c(0,4.1,0,0), oma = c(5,0,4,2), mfrow = c(1,1))
    on.exit(par(opar))
    if(is.null(BSseqStat) && is.null(geneTrack)) {
        layout(matrix(1:2, ncol = 1), heights = c(2,1))
    } else if (is.null(geneTrack)) {
        layout(matrix(1:3, ncol = 1), heights = c(2,2,1))
    } else {
        layout(matrix(1:4, ncol = 1), heights = c(2,2,1,0.3))
    }

    .plotSmoothData(BSseq = BSseq, region = region, extend = extend, addRegions = addRegions,
                    col = col, lty = lty, lwd = lwd, regionCol = regionCol,
                    addTicks = addTicks, addPoints = addPoints,
                    pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    gr <- .bsGetGr(BSseq, region, extend)

    if(!is.null(BSseqStat)) {
        BSseqStat <- subsetByOverlaps(BSseqStat, gr)
        if(is(BSseqStat, "BSseqTstat")) {
            stat.values <- as.array(getStats(BSseqStat)[, "tstat.corrected"])
            stat.values <- as.array(stat.values)
            stat.type <- "tstat"
        }
        if(is(BSseqStat, "BSseqStat")) {
            stat.type <- getStats(BSseqStat, what = "stat.type")
            if(stat.type == "tstat") {
                stat.values <- getStats(BSseqStat, what = "stat")
                stat.values <- as.array(stat.values)
            }
            if(stat.type == "fstat") {
                stat.values <- sqrt(getStats(BSseqStat, what = "stat"))
                stat.values <- as.array(stat.values)
            }
        }
        plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n",
             ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "", ylab = stat.type)
        axis(side = 2, at = c(-5,0,5))
        abline(h = 0, col = "grey60")
        .bsPlotLines(start(BSseqStat), stat.values, lty = stat.lty, col = stat.col, lwd = stat.lwd,
                     plotRange = c(start(gr), end(gr)))
    }

    if(!is.null(annoTrack))
        plotAnnoTrack(gr, annoTrack, cex.anno)

    if (!is.null(geneTrack))
        plotGeneTrack(gr, geneTrack, cex.gene)

    if(!is.null(main)) {
        main <- .bsPlotTitle(gr = region, extend = extend, main = main,
                             mainWithWidth = mainWithWidth)
        mtext(side = 3, text = main, outer = TRUE, cex = 1)
    }
    return(invisible(NULL))
}

##     plotP <- function(sample, label = "", col) {
##         plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
##              ylim = c(0,1), xlim = plotRange, xlab = "", ylab = "")
##         plotRects(c(0,1))
##         rawp <- getP(BSseq, sample = sample, type = "raw", addPositions = TRUE, addConfint = TRUE)
##         cols <- rep(alpha("black", 0.5), nrow(rawp))
##         segments(x0 = rawp$pos, y0 = rawp$lower, y1 = rawp$upper, col = cols)
##         points(positions, rawp$p, col = cols)
##         if(nrow(BSseq$coef) > 0) {
##             fitp <- getP(BSseq, sample = sample, type = "fit", addPosition = TRUE, addConfint = FALSE)
##             lines(fitp$pos, fitp$p, col = col,, lty = 1)
##             ## lines(fitp$pos, fitp$lower, col = col, lty = 2)
##             ## lines(fitp$pos, fitp$upper, col = col, lty = 2)
##             text(plotRange[1], 0.1, labels = label)
##         }
##     }

##     plotD <- function(sample1, sample2, col) {
##         plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
##              ylim = c(-1,1), xlim = plotRange, xlab = "", ylab = "")
##         plotRects(c(-1,1))
##         abline(h = 0, col = alpha("black", 0.6), lty = 2)
##         rawd <- getD(BSseq, sample1 = sample1, sample2 = sample2, type = "raw",
##                      addPositions = TRUE)
##         cols <- rep(alpha("black", 0.5), nrow(rawd))
##         segments(x0 = rawd$pos, rawd$lower, y1 = rawd$upper, col = cols)
##         points(rawd$pos, rawd$d, col = cols)
##         if(nrow(BSseq$coef) > 0) {
##             fitd <- getD(BSseq, sample1 = sample1, sample2 = sample2, type = "fit", addPositions = TRUE)
##             lines(fitd$pos, fitd$d, col = col, lty = 1)
##             ## lines(fitd$pos, fitd$lower, col = "blue", lty = 2)
##             ## lines(fitd$pos, fitd$upper, col = "blue", lty = 2)
##         }
##     }
