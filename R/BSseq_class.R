setClass("BSseq", contains = "SummarizedExperiment", 
         representation(trans = "function",
                        parameters = "list"))

setValidity("BSseq", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Cov", "M")))
    if(class(rowData(object)) != "GRanges")
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
    ## benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
    if(is.null(colnames(object)))
        msg <- validMsg(msg, "colnames (aka sampleNames) need to be set")
    if(min(assay(object, "M")) < 0)
        msg <- validMsg(msg, "the 'M' assay has negative entries")
    if(min(assay(object, "Cov")) < 0)
        msg <- validMsg(msg, "the 'Cov' assay has negative entries")
    if(max(assay(object, "M") - assay(object, "Cov")) > 0.5)
        msg <- validMsg(msg, "the 'M' assay has at least one entry bigger than the 'Cov' assay")
    if(!is.null(rownames(assay(object, "M"))) ||
       !is.null(rownames(assay(object, "Cov"))) ||
       ("coef" %in% assayNames(object) && !is.null(rownames(assay(object, "coef")))) ||
       ("se.coef" %in% assayNames(object) && !is.null(rownames(assay(object, "se.coef")))))
        warning("unnecessary rownames in object")
    if (is.null(msg)) TRUE else msg
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

setMethod("pData", "BSseq", function(object) {
    object@colData
})

setReplaceMethod("pData",
                 signature = signature(
                      object = "BSseq",
                      value = "data.frame"),
                 function(object, value) {
                     colData(object) <- as(value, "DataFrame")
                     object
                 })

setReplaceMethod("pData",
                 signature = signature(
                      object = "BSseq",
                      value = "DataFrame"),
                 function(object, value) {
                     colData(object) <- value
                     object
                 })

setMethod("sampleNames", "BSseq", function(object) {
    colnames(object)
})

setReplaceMethod("sampleNames",
                 signature = signature(
                     object = "BSseq",
                     value = "ANY"),
                 function(object, value) {
                     colnames(object) <- value
                     object
                 })

setMethod("length", "BSseq", function(x) {
    length(granges(x))
})

hasBeenSmoothed <- function(BSseq) {
    "coef" %in% assayNames(BSseq)
}

getBSseq <- function(BSseq, type = c("Cov", "M", "gr", "coef", "se.coef", "trans", "parameters")) {
    type <- match.arg(type)
    if(type %in% c("M", "Cov"))
        return(assay(BSseq, type))
    if(type %in% c("coef", "se.coef") && type %in% assayNames(BSseq))
        return(assay(BSseq, type))
    if(type %in% c("coef", "se.coef"))
       return(NULL)
    if(type == "trans")
        return(BSseq@trans)
    if(type == "parameters")
        return(BSseq@parameters)
    if(type == "gr")
        return(BSseq@rowData)
    
}

BSseq <- function(M = NULL, Cov = NULL, coef = NULL, se.coef = NULL,
                  trans = NULL, parameters = NULL, pData = NULL,
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
    if(!is(pData, "DataFrame"))
        pData <- as(pData, "DataFrame")
    if(is.null(sampleNames) && !is.null(pData) && !is.null(rownames(pData)))
        sampleNames <- rownames(pData)
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
    if(!is.null(coef)) {
        if(!is.matrix(coef) ||
           nrow(coef) != nrow(M) ||
           ncol(coef) != ncol(M))
            stop("'coef' does not have the right dimensions")
        if(is.null(colnames(coef)) || any(sampleNames != colnames(coef)))
            colnames(coef) <- sampleNames
        if(!is.null(rownames(coef)))
            rownames(coef) <- NULL
    }
    if(!is.null(se.coef)) {
        if(!is.matrix(se.coef) ||
           nrow(se.coef) != nrow(M) ||
           ncol(se.coef) != ncol(M))
            stop("'se.coef' does not have the right dimensions")
        if(is.null(colnames(se.coef)) || any(sampleNames != colnames(se.coef)))
            colnames(se.coef) <- sampleNames
        if(!is.null(rownames(se.coef)))
            rownames(se.coef) <- NULL
    }
    assays <- SimpleList(M = M, Cov = Cov, coef = coef, se.coef = se.coef)
    assays <- assays[!sapply(assays, is.null)]
    if(is.null(pData) || all(dim(pData) == c(0,0)))
        BSseq <- SummarizedExperiment(assays = assays, rowData = gr)
    else
        BSseq <- SummarizedExperiment(assays = assays, rowData = gr, colData = pData)
    BSseq <- as(BSseq, "BSseq")
    if(is.function(trans))
        BSseq@trans <- trans
    if(is.list(parameters))
        BSseq@parameters <- parameters
    BSseq
}


setMethod("updateObject", "BSseq",
          function(object, ...) {
               if(!is(try(object@assays, silent = TRUE), "try-error"))
                   return(object)
               BSseq(gr = object@gr, M = object@M, Cov = object@Cov,
                     coef = object@coef, se.coef = object@se.coef,
                     trans = object@trans, parameters = object@parameters,
                     pData = object@phenoData@data)
           })



                  

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
