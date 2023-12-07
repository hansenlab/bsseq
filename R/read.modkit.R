read.modkit <- function(files,
                        colData = NULL,
                        rmZeroCov = FALSE,
                        strandCollapse = TRUE,
                        Other_mod_type = c("drop", "additional M",
                                           "additional Cov", "store separately")){
    gr_list <- list()
    sampleNames <- sub("\\.bed.gz$", "", basename(files))
    if (!is.null(colData)){
        rownames(colData) <- sampleNames
        }

    for (i in seq_along(files)){
        data <- read.table(files[i], header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="")
        gr <- GRanges(seqnames = data$V1,
                  ranges = IRanges(start = data$V2+1, end = data$V3))
        mcols(gr)$m <- data$V12
        mcols(gr)$u <- data$V13
        mcols(gr)$filter <- data$V16
        if (any(data$V14 != 0)){
            mcols(gr)$other_mod <- data$V14
            }
        names(gr) <- sampleNames[i]
        gr_list[[sampleNames[i]]] <- gr
        }

    overlap_gr <- Reduce(subsetByOverlaps, gr_list)

    m_u_list <- lapply(gr_list, function(gr){
        overlap_data <- gr[gr %over% overlap_gr]
        if (!is.null(gr$other_mod)) {
            data.frame(m = overlap_data$m, u = overlap_data$u,
                       filter = overlap_data$filter,
                       other_mod = overlap_data$other_mod)
            } else {
                data.frame(m = overlap_data$m,  u = overlap_data$u,
                           filter = overlap_data$filter)
                }
        })

    m <- do.call(cbind, lapply(m_u_list, `[[`, "m"))
    u <- do.call(cbind, lapply(m_u_list, `[[`, "u"))
    filter <- do.call(cbind, lapply(m_u_list, `[[`, "filter"))
    other_mod <- do.call(cbind, lapply(m_u_list, `[[`, "other_mod"))

    Other_mod_type <- match.arg(Other_mod_type)

    if (!is.null(other_mod)){
      if (Other_mod_type == "drop"){
        bsseq_obj <- BSseq(M = as.matrix(m), Cov = as.matrix(u + m),
                           Filtered = as.matrix(filter),
                           coef = NULL, se.coef = NULL,
                           pos = start(overlap_gr), trans = NULL,
                           parameters = NULL, pData = colData, gr = NULL,
                           chr = as.vector(seqnames(overlap_gr)),
                           sampleNames = sampleNames, rmZeroCov = rmZeroCov)
        if (strandCollapse) {strandCollapse(bsseq_obj)}
      } else if (Other_mod_type == "additional M"){
          bsseq_obj <- BSseq(M = as.matrix(m + other_mod),
                             Cov = as.matrix(m + u +other_mod),
                             Filtered = as.matrix(filter),
                             coef = NULL, se.coef = NULL,
                             pos = start(overlap_gr), trans = NULL,
                             parameters = NULL, pData = colData, gr = NULL,
                             chr = as.vector(seqnames(overlap_gr)),
                             sampleNames = sampleNames, rmZeroCov = rmZeroCov)
          if (strandCollapse) {strandCollapse(bsseq_obj)}
      } else if (Other_mod_type == "additional Cov"){
          bsseq_obj <- BSseq(M = as.matrix(m), Cov = as.matrix(m + u + other_mod),
                             Filtered = as.matrix(filter),
                             coef = NULL, se.coef = NULL,
                             pos = start(overlap_gr), trans = NULL,
                             parameters = NULL, pData = colData, gr = NULL,
                             chr = as.vector(seqnames(overlap_gr)),
                             sampleNames = sampleNames, rmZeroCov = rmZeroCov)
          if (strandCollapse) {strandCollapse(bsseq_obj)}
      } else if (Other_mod_type == "store separately"){
          bsseq_obj <- BSseq_mod(M = as.matrix(m), U = as.matrix(u),
                                 Filtered = as.matrix(filter),
                                 Other_mod = as.matrix(other_mod),
                                 coef = NULL, se.coef = NULL,
                                 pos = start(overlap_gr), trans = NULL,
                                 parameters = NULL, pData = colData, gr = NULL,
                                 chr = as.vector(seqnames(overlap_gr)),
                                 sampleNames = sampleNames, rmZeroCov = rmZeroCov)
          if (strandCollapse) {strandCollapse_mod(bsseq_obj)}
      }
    }else{
      bsseq_obj <- BSseq(M = as.matrix(m),Cov = as.matrix(u+m),
                         Filtered = as.matrix(filter),
                         coef = NULL, se.coef = NULL,
                         pos = start(overlap_gr), trans = NULL,
                         parameters = NULL, pData = colData, gr = NULL,
                         chr = as.vector(seqnames(overlap_gr)),
                         sampleNames = sampleNames, rmZeroCov = rmZeroCov)
      if (strandCollapse) {strandCollapse(bsseq_obj)}
    }

    return(bsseq_obj)
}
