fisherTests <- function(BSseq, group1, group2, lookup = NULL,
                        returnLookup = TRUE, mc.cores = 1,
                        verbose = TRUE) {
    stopifnot(is(BSseq, "BSseq"))
    if(is.character(group1)) {
        stopifnot(all(group1 %in% sampleNames(BSseq)))
        group1 <- match(group1, sampleNames(BSseq))
    }
    if(is.numeric(group1)) {
        stopifnot(min(group1) >= 1 & max(group1) <= ncol(BSseq))
    } else stop("problems with argument 'group1'")
    if(is.character(group2)) {
        stopifnot(all(group2 %in% sampleNames(BSseq)))
        group2 <- match(group2, sampleNames(BSseq))
    }
    if(is.numeric(group2)) {
        stopifnot(min(group2) >= 1 & max(group2) <= ncol(BSseq))
    } else stop("problems with argument 'group2'")
    stopifnot(length(intersect(group1, group2)) == 0)
    if(is.null(lookup)) {
        oldkeys <- ""
    } else {
        stopifnot(is.matrix(lookup) && ncol(lookup) == 2 &&
                  setequal(colnames(lookup), c("p.value", "log2OR")))
        oldkeys <- rownames(lookup)
    }
    if(verbose) cat("[fisherTests] setup ... ")
    ptime1 <- proc.time()
    M1 <- rowSums2(getCoverage(BSseq, type = "M")[, group1, drop = FALSE])
    M2 <- rowSums2(getCoverage(BSseq, type = "M")[, group2, drop = FALSE])
    Cov1 <- rowSums2(getCoverage(BSseq, type = "Cov")[, group1, drop = FALSE])
    Cov2 <- rowSums2(getCoverage(BSseq, type = "Cov")[, group2, drop = FALSE])
    keys <- paste(Cov1, Cov2, M1, M2, sep = "_")
    newkeys <- setdiff(keys, oldkeys)
    expand <- matrix(as.integer(do.call(rbind, strsplit(newkeys, "_"))), ncol = 4)
    colnames(expand) <- c("Cov1", "Cov2", "M1", "M2")
    expand2 <- expand
    colnames(expand2) <- c("U1", "U2", "M1", "M2")
    expand2[,"U1"] <- expand[,"Cov1"] - expand[,"M1"]
    expand2[,"U2"] <- expand[,"Cov2"] - expand[,"M2"]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(round(stime,1), "done in secs\n")
    if(verbose) cat(sprintf("[fisherTests] performing %d tests (using %d cores) ... ",
                            nrow(expand2), mc.cores))
    ptime1 <- proc.time()
    tests <- mclapply(1:nrow(expand2), function(ii) {
        unlist(fisher.test(matrix(expand2[ii,], ncol = 2))[c("p.value", "estimate")])
    }, mc.cores = mc.cores)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(round(stime,1), "done in secs\n")
    if(verbose) cat("[fisherTests] finalizing ... ")
    ptime1 <- proc.time()
    tests <- do.call(rbind, tests)
    tests[,2] <- log2(tests[,2])
    colnames(tests) <- c("p.value", "log2OR")
    rownames(tests) <- newkeys
    lookup <- rbind(tests, lookup)
    results <- lookup[keys,]
    rownames(results) <- NULL
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(round(stime,1), "done in secs\n")
    if(returnLookup)
        return(list(results = results, lookup = lookup))
    else
        return(results)
}


