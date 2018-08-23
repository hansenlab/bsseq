getSex <- function(BSseq, cutoff = 0){
	stopifnot(is(BSseq, "BSseq"))
	xChr <- match.arg(unique(as.character(seqnames(BSseq))), c("X","chrX"), several.ok = TRUE)
	yChr <- match.arg(unique(as.character(seqnames(BSseq))), c("Y","chrY"), several.ok = TRUE)
	
	BSseq.x <- chrSelectBSseq(BSseq, seqnames = xChr)
	BSseq.y <- chrSelectBSseq(BSseq, seqnames = yChr)
	
	xCov <- getCoverage(BSseq.x, what = "perRegionAverage")
	yCov <- getCoverage(BSseq.y, what = "perRegionAverage")
	metric <- xCov - yCov
	sexclust <- kmeans(metric, 2)
	
	if(all(sexclust$centers > cutoff)){
		message("[getSex] No samples with coverage of Y chromosome. Predicting all female samples.")
		df <- DataFrame(xCov = xCov, yCov = yCov, predictedSex = rep("F", ncol(BSseq)))
		rownames(df) <- sampleNames(BSseq)
		stop(return(df))
	}
	if(all(sexclust$centers < cutoff)){
		message("[getSex] No samples without coverage of Y chromosme. Predicting all male samples.")
		df <- DataFrame(xCov = xCov, yCov = yCov, predictedSex = rep("M", ncol(BSseq)))
		rownames(df) <- sampleNames(BSseq)
		stop(return(df))
	}
	predictedSex <- ifelse(metric > cutoff, "F", "M")
	
	xClust <- kmeans(xCov, centers=c(min(xCov), max(xCov)))
	yClust <- kmeans(yCov, centers=c(max(yCov), min(yCov)))
	
	predSex <- rep(NA, ncol(BSseq))
	names(predSex) <- sampleNames(BSseq)
	predSex[(xClust$cluster + yClust$cluster) ==  2] <- "M"
	predSex[(xClust$cluster + yClust$cluster) ==  4] <- "F"
	
	if(!identical(predictedSex, predSex))
		warning("[getSex] Possible discrepancy with Sample ID ",paste(sampleNames(BSseq)[(predictedSex!=predSex | is.na(predSex))],collapse=", "),". Check via plotSex().")

	df <- DataFrame(xCov = xCov, yCov = yCov, predictedSex = predictedSex)
    rownames(df) <- sampleNames(BSseq)
    df
}

addSex <- function(BSseq, sex = NULL) {
    if(is.null(sex))
        sex <- getSex(BSseq)$predictedSex
    if(is(sex, "DataFrame") && "predictedSex" %in% names(sex))
        sex <- sex$predictedSex
	pd <- pData(BSseq)
	if ("predictedSex" %in% names(pd))
		warning("Replacing existing 'predictedSex' column")
	pd$predictedSex <- sex
	pData(BSseq) <- pd
	BSseq
}

plotSex <- function(df, id = NULL) {
    stopifnot(all(c("predictedSex", "xCov", "yCov") %in% names(df)))
    if(is.null(id))
        id <- rownames(df)
    if(length(id) != length(df$predictedSex))
        stop("id length must match number of samples.")
    plot(df$xCov, df$yCov, type = "n",
         xlab = "X chr, Average coverage",
         ylab = "Y chr, Average coverage")
	pointcol <- rep("black",nrow(df))
	pointcol[df$predictedSex == "M"] <- "deepskyblue"
	pointcol[df$predictedSex == "F"] <- "deeppink3"
    text(df$xCov, df$yCov, id,
         col=ifelse(df$predictedSex == "M", "deepskyblue", "deeppink3"))
    legend("bottomleft", c("M","F"), col = c("deepskyblue", "deeppink3"), pch = 16)
}