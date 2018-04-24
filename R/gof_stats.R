poissonGoodnessOfFit <- function(BSseq, nQuantiles = 10^5) {
    Cov <- getBSseq(BSseq, type = "Cov")
    lambda <- rowMeans2(Cov)
    out <- list(chisq = NULL, df = ncol(Cov) - 1)
    out$chisq <- rowSums2((Cov - lambda)^2 / lambda)
    if(nQuantiles > length(out$chisq))
        nQuantiles <- length(out$chisq)
    probs <- seq(from = 0, to = 1, length.out = nQuantiles)
    out$quantiles <- quantile(out$chisq, prob = probs, na.rm = TRUE)
    class(out) <- "chisqGoodnessOfFit"
    return(out)
}

print.chisqGoodnessOfFit <- function(x, ...) {
    cat("Chisq Goodness of Fit object\n")
    cat(" ", length(x$chisq), "tests\n")
    cat(" ", x$df, "degrees of freedom\n")
}

binomialGoodnessOfFit <- function(BSseq, method = c("MLE"), nQuantiles = 10^5) {
    method <- match.arg(method)
    M <- getCoverage(BSseq, type = "M", what = "perBase")
    Cov <- getCoverage(BSseq, type = "Cov", what = "perBase")
    if(method == "MLE") # This is the MLE under iid
        p <- rowSums2(M) / rowSums2(Cov)
    ## p <- rowMeans(M / Cov, na.rm = TRUE)
    out <- list(chisq = NULL, quantiles = NULL, df = ncol(Cov) - 1)
    ## This gets tricky.  For the binomial model, we need Cov > 0, but
    ## this implies that a given location may have less then nCol(BSseq)
    ## number of observations, and hence the degrees of freedom will
    ## vary.
    out$chisq <- rowSums2((M - Cov*p)^2 / sqrt(Cov * p * (1-p)))
    probs <- seq(from = 0, to = 1, length.out = nQuantiles)
    ## Next line will remove all locations where not all samples have coverage
    out$quantiles <- quantile(out$chisq, prob = probs, na.rm = TRUE)
    class(out) <- "chisqGoodnessOfFit"
    return(out)
}

plot.chisqGoodnessOfFit <- function(x, type = c("chisq", "pvalue"),
                                    plotCol = TRUE, qqline = TRUE, pch = 16, cex = 0.75, ...) {
    type <- match.arg(type)
    switch(type,
           "chisq" = {
               yy <- x$quantiles
               xx <- qchisq(ppoints(yy), df = x$df)
           },
           "pvalue" = {
               yy <- sort(1 - qchisq(x$quantiles, df = x$df))
               xx <- qunif(ppoints(yy), min = 0, max = 1)
           })
    if(plotCol) {
        colors <- rep("black", length(yy))
        colors[floor(length(colors)*.95):length(colors)] <- "red"
        colors[floor(length(colors)*.99):length(colors)] <- "violet"
        colors[floor(length(colors)*.999):length(colors)] <- "orange"
    } else {
        colors <- "black"
    }
    qqplot(xx, yy, xlab = "Theoretical quantiles",
           ylab = "Observed quantiles", col = colors, pch = pch, cex = cex, ...)
    if(qqline) {
        yyqt <- quantile(yy, c(0.25, 0.75))
        xxqt <- quantile(xx, c(0.25, 0.75))
        slope <- diff(yyqt) / diff(xxqt)
        int <- yyqt[1L] - slope * xxqt[1L]
        abline(int, slope, col = "blue")
    }
    abline(0, 1, col = "blue", lty = 2)
}

