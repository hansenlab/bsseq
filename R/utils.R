data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                 "DataFrame")
        mcols(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

.ON_DISK_SEEDS <- c("HDF5ArraySeed")
.ON_DISK_BACKENDS <- c("HDF5Array")

.areBackendsInMemory <- function(realization_backends) {
    if (is.null(realization_backends)) {
        return(TRUE)
    }
    vapply(realization_backends, function(realization_backend) {
        is.null(realization_backend) ||
            !realization_backend %in% .ON_DISK_BACKENDS
    }, logical(1L))
}

.getBSseqBackends <- function(x) {
    assay_backends <- lapply(assays(x, withDimnames = FALSE), function(assay) {
        if (is.matrix(assay)) return(NULL)
        seed_classes <- .getSeedClasses(assay)
        if (all(vapply(seed_classes, function(x) x == "matrix", logical(1)))) {
            return(NULL)
        }
        if (is.list(seed_classes)) {
            seed_packages <- lapply(seed_classes, attr, "package")
        } else {
            seed_packages <- attr(seed_classes, "package")
        }
        seed_packages <- unique(seed_packages)
        srb <- supportedRealizationBackends()
        srb[srb[["package"]] == seed_packages, "BACKEND"]
    })
    unique(unlist(assay_backends))
}

# TODO: https://github.com/Bioconductor/BiocParallel/issues/76
.isSingleMachineBackend <- function(BPPARAM) {
    if (is(BPPARAM, "SerialParam") || is(BPPARAM, "MulticoreParam")) {
        return(TRUE)
    } else if (is(BPPARAM, "SnowParam")) {
        if (is.numeric(bpworkers(BPPARAM)) &&
            BPPARAM$.clusterargs$type == "SOCK") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else if (is(BPPARAM, "DoparParam")) {
        # TODO: Can't figure this one out, so returning FALSE for now
        return(FALSE)
    } else {
        return(FALSE)
    }
}
