read.modkit = function(files,
                       colData = NULL,
                       rmZeroCov = FALSE,
                       strandCollapse = TRUE){
  gr_list = list()
  sampleNames = sub("\\.bed$", "", basename(files))
  if (!is.null(colData)){
    rownames(colData) <- sampleNames
  }

  for (i in seq_along(files)){
    data = read.table(files[i], header = FALSE, sep="\t",
                      stringsAsFactors=FALSE, quote="")

    gr = GRanges(seqnames = data$V1,
                 ranges = IRanges(start = data$V2+1, end = data$V3))
    mcols(gr)$m = data$V12
    mcols(gr)$cov = data$V12 + data$V13
    mcols(gr)$filter = data$V16
    if (any(data$V14 != 0)){
        mcols(gr)$other_mod = data$V14
        mcols(gr)$cov_other_mod = data$V14 + data$V13}
    names(gr) = sampleNames[i]
    gr_list[[sampleNames[i]]] = gr
  }

  overlap_gr = Reduce(subsetByOverlaps, gr_list)

  m_cov_list = lapply(gr_list, function(gr){
    overlap_data = gr[gr %over% overlap_gr]
    if (!is.null(gr$other_mod)) {
      data.frame(m = overlap_data$m, cov = overlap_data$cov,
                 filter = overlap_data$filter,
                 other_mod = overlap_data$other_mod,
                 cov_other_mod = overlap_data$cov_other_mod)
      } else {
      data.frame(m = overlap_data$m, cov = overlap_data$cov,
                 filter = overlap_data$filter)
      }
    })

  m = do.call(cbind, lapply(m_cov_list, `[[`, "m"))
  cov = do.call(cbind, lapply(m_cov_list, `[[`, "cov"))
  filter = do.call(cbind, lapply(m_cov_list, `[[`, "filter"))
  other_mod = do.call(cbind, lapply(m_cov_list, `[[`, "other_mod"))
  cov_other_mod = do.call(cbind, lapply(m_cov_list, `[[`, "cov_other_mod"))

  bsseq_objs = list()
  bsseq_obj_m = BSseq(M = as.matrix(m), Cov = as.matrix(cov),
                    Filtered = as.matrix(filter),
                    coef = NULL, se.coef = NULL,
                    pos = start(overlap_gr), trans = NULL,
                    parameters = NULL, pData = colData, gr = NULL,
                    chr = as.vector(seqnames(overlap_gr)),
                    sampleNames = sampleNames, rmZeroCov = rmZeroCov)

  bsseq_objs[[1]] = bsseq_obj_m

  if (!is.null(other_mod)){
      bsseq_obj_other_mod = BSseq(M = as.matrix(other_mod),
                                  Cov = as.matrix(cov_other_mod),
                                  Filtered = as.matrix(filter),
                                  coef = NULL, se.coef = NULL,
                                  pos = start(overlap_gr), trans = NULL,
                                  parameters = NULL, pData = colData, gr = NULL,
                                  chr = as.vector(seqnames(overlap_gr)),
                                  sampleNames = sampleNames, rmZeroCov = rmZeroCov)

      bsseq_objs[[2]] = bsseq_obj_other_mod
    }

  if (strandCollapse){
      for (i in 1:length(bsseq_objs)) {
      bsseq:::strandCollapse(bsseq_objs[[i]])}
  }

  return(bsseq_objs)
}
