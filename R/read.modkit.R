read.modkit <- function(files,
                        colData = NULL,
                        rmZeroCov = FALSE,
                        strandCollapse = TRUE){
    gr_list <- list()
    sampleNames <- sub("\\.bed.gz$", "", basename(files))
    if (!is.null(colData)){
        rownames(colData) <- sampleNames
        }

    for (i in seq_along(files)){
        data <- read.table(files[i], header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="")
        data$V6[data$V6 == "."] <- "*"

        if (length(unique(data$V4)) == 2){
            gr <- GRanges(seqnames = data[data$V4 == "m", ]$V1,
                    ranges = IRanges(start = data[data$V4 == "m", ]$V2+1,
                                     end = data[data$V4 == "m", ]$V3),
                    strand = data[data$V4 == "m", ]$V6)

            mcols(gr)$m <- data[data$V4 == "m", ]$V12
            mcols(gr)$h <- data[data$V4 != "m", ]$V12
            mcols(gr)$u <- data[data$V4 == "m", ]$V13
            mcols(gr)$filter <- data[data$V4 == "m", ]$V16
        }else{
            gr <- GRanges(seqnames = data$V1,
                          ranges = IRanges(start = data$V2+1, end = data$V3),
                          strand = data$V6)

            mcols(gr)$m <- data$V12
            mcols(gr)$u <- data$V13
            mcols(gr)$filter <- data$V16
        }

        names(gr) <- sampleNames[i]
        gr_list[[sampleNames[i]]] <- gr
        }

    overlap_gr <- Reduce(subsetByOverlaps, gr_list)

    m_u_list <- lapply(gr_list, function(gr){
        overlap_data <- gr[gr %over% overlap_gr]
        if (!is.null(gr$h)) {
            data.frame(m = overlap_data$m, u = overlap_data$u,
                       h = overlap_data$h, filter = overlap_data$filter)
            } else {
                data.frame(m = overlap_data$m,  u = overlap_data$u,
                           filter = overlap_data$filter)
                }
        })

    m <- do.call(cbind, lapply(m_u_list, `[[`, "m"))
    u <- do.call(cbind, lapply(m_u_list, `[[`, "u"))
    h <- do.call(cbind, lapply(m_u_list, `[[`, "h"))
    filter <- do.call(cbind, lapply(m_u_list, `[[`, "filter"))

    if (!is.null(h)){
          bsseq_obj <- BSseq(M = as.matrix(m + h),
                             Cov = as.matrix(m + u + h),
                             Filtered = as.matrix(filter),
                             coef = NULL, se.coef = NULL,
                             pos = start(overlap_gr), trans = NULL,
                             parameters = NULL, pData = colData, gr = NULL,
                             chr = as.vector(seqnames(overlap_gr)),
                             sampleNames = sampleNames, rmZeroCov = rmZeroCov)
          if (strandCollapse) {
            bsseq_obj <- strandCollapse(bsseq_obj)
          }
    }else{
      bsseq_obj <- BSseq(M = as.matrix(m), Cov = as.matrix(u + m),
                         Filtered = as.matrix(filter),
                         coef = NULL, se.coef = NULL,
                         pos = start(overlap_gr), trans = NULL,
                         parameters = NULL, pData = colData, gr = NULL,
                         chr = as.vector(seqnames(overlap_gr)),
                         sampleNames = sampleNames, rmZeroCov = rmZeroCov)
      if (strandCollapse) {
        bsseq_obj <- strandCollapse(bsseq_obj)
      }
    }

    return(bsseq_obj)
}
