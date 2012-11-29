plotAnnoTrack <- function(gr, annoTrack) {
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
        mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1)
        })
}

plotManyRegions <- function(BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL, annoTrack = NULL, 
                            col = NULL, lty = NULL, lwd = NULL, BSseqTstat = NULL, mainWithWidth = TRUE,
                            regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE,
                            pointsMinCov = 5, highlightMain = FALSE, verbose = TRUE) {
    cat("preprocessing ...")
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
    if(!is.null(BSseqTstat))
        BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
    
    if(length(start(BSseq)) == 0)
        stop("No overlap between BSseq data and regions")
    if(!is.null(main) && length(main) != length(gr))
        main <- rep(main, length = length(gr))
    cat("done\n")
    for(ii in seq(along = gr)) {
        if(verbose) cat(sprintf("plotting region %d (out of %d)\n", ii, nrow(regions)))
        plotRegion(BSseq = BSseq, region = regions[ii,], extend = extend,
                   col = col, lty = lty, lwd = lwd, main = main[ii], BSseqTstat = BSseqTstat,
                   addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth,
                   annoTrack = annoTrack, addTicks = addTicks, addPoints = addPoints,
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

.bsPlotPoints <- function(x, y, z, col) {
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
    if(length(regions) == 0)
        return(NULL)
    rect(xleft = start(regions), xright = end(regions), ybottom = ylim[1],
         ytop = ylim[2], col = regionCol, border = NA)
}

.bsGetCol <- function(object, col, lty, lwd) {
    ## Assumes that object has pData and sampleNames methods
    if(is.null(col) && "col" %in% names(pData(object)))
        col <- pData(object)[["col"]]
    else
        col <- rep("black", nrow(pData(object)))
    if(is.null(names(col)))
        names(col) <- sampleNames(object)
    
    if(is.null(lty) && "lty" %in% names(pData(object)))
        lty <- pData(object)[["lty"]]
    else
        lty <- rep(1, nrow(pData(object)))
    if(is.null(names(lty)))
        names(lty) <- sampleNames(object)
    
    if(is.null(lwd) && "lwd" %in% names(pData(object)))
        lwd <- pData(object)[["lwd"]]
    else
        lwd <- rep(1, nrow(pData(object)))
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
                               format(width(gr) - 2*extend, big.mark = ",", scientific = FALSE),
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
    gr <- .bsGetGr(Bsseq, region, extend)
    BSseq <- subsetByOverlaps(BSseq, gr)
    
    ## Extract basic information
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    positions <- start(BSseq)
    smoothPs <- getMeth(BSseq, type = "smooth")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)
        
    ## get col, lwd, lty
    colEtc <- bsseq:::.bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)
    
    ## The actual plotting
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = c(0,1), xlim = c(start(gr), end(gr)), xlab = "", ylab = "Methylation")
    axis(side = 2, at = c(0.2, 0.5, 0.8))
    if(addTicks)
        rug(positions)

    .bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0,1),
                        regionCol = regionCol, highlightMain = highlightMain)
    
    if(addPoints) {
        sapply(sampleNames(BSseq), function(samp) {
            abline(v = positions[rawPs[, samp] > 0.1], col = "grey80", lty = 1)
        })
    } # This adds vertical grey lines so we can see where points are plotted

    sapply(sampleNames(BSseq), function(samp) {
        .bsPlotLines(positions, smoothPs[, samp], col = colEtc$col[samp],
                     lty = colEtc$lty[samp], lwd = colEtc$lwd[samp],
                     plotRange = c(start(gr), end(gr)))
    })

    if(addPoints) {
        sapply(sampleNames(BSseq), function(samp) {
            .bsPlotPoints(positions, rawPs[, samp], coverage[, samp],
                          col = colEtc$col[samp])
        })
    }
}


plotRegion <- function(BSseq, region = NULL, extend = 0, main = "", addRegions = NULL, annoTrack = NULL,
                          col = NULL, lty = NULL, lwd = NULL, BSseqTstat = NULL, mainWithWidth = TRUE,
                          regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE,
                          pointsMinCov = 5, highlightMain = FALSE) {
    
    opar <- par(mar = c(0,4.1,0,0), oma = c(5,0,4,2), mfrow = c(1,1))
    on.exit(par(opar))
    if(is.null(BSseqTstat))
        layout(matrix(1:2, ncol = 1), heights = c(2,1))
    else
        layout(matrix(1:3, ncol = 1), heights = c(2,2,1))

    bsseq:::.plotSmoothData(BSseq = BSseq, region = region, extend = extend, addRegions = addRegions,
                            col = col, lty = lty, lwd = lwd, regionCol = regionCol,
                            addTicks = addTicks, addPoints = addPoints,
                            pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    gr <- .bsGetGr(BSseq, region, extend)
    
    if(!is.null(BSseqTstat)) {
        if(!is.null(BSseqTstat))
            BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
        plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
             ylim = c(-8,8), xlim = plotRange, xlab = "", ylab = "t-stat")
        axis(side = 2, at = c(-5,0,5))
        abline(h = 0, col = "grey60")
        bsseq:::.bsPlotLines(start(BSseqTstat), BSseqTstat@stats[, "tstat"],
                             lty = 1, plotRange = plotRange, col = "red", lwd = 1)
        bsseq:::.bsPlotLines(start(BSseqTstat), BSseqTstat@stats[, "tstat.corrected"],
                             lty = 2, plotRange = plotRange, col = "red", lwd = 1)
        bsseq:::.bsPlotLines(start(BSseqTstat), 100*BSseqTstat@stats[, "tstat.sd"],
                             lty = 2, plotRange = plotRange, col = "blue", lwd = 1)
    }
    
    if(!is.null(annoTrack))
        bsseq:::plotAnnoTrack(gr, annoTrack)

    if(!is.null(main)) {
        main <- bsseq:::.bsPlotTitle(gr = region, extend = extend, main = main,
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
