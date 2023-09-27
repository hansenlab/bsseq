read.modbam2bed = function(files,
                           colData = NULL,
                           rmZeroCov = FALSE,
                           strandCollapse = TRUE){
  gr_list = list()
  sampleNames = sub("\\.bed$","",basename(files))
  if (!is.null(colData)){
    rownames(colData) <- sampleNames
  }

  for (i in seq_along(files)){
    data = read.table(files[i],header = FALSE, sep="\t",
                      stringsAsFactors=FALSE, quote="")
    gr = GRanges(seqnames = data$V1,
                 ranges = IRanges(start = data$V2+1,end = data$V3))
    mcols(gr)$m = data$V13
    mcols(gr)$cov = data$V12 + data$V13
    mcols(gr)$filter = data$V14
    names(gr) = sampleNames[i]
    gr_list[[sampleNames[i]]] = gr
  }

  overlap_gr = Reduce(subsetByOverlaps, gr_list)

  m_cov_list = lapply(gr_list, function(gr){
    overlap_data = gr[gr %over% overlap_gr]
    data.frame(m = overlap_data$m, cov = overlap_data$cov,
               filter = overlap_data$filter)})

  m = do.call(cbind,lapply(m_cov_list,`[[`, "m"))
  cov = do.call(cbind, lapply(m_cov_list, `[[`, "cov"))
  filter = do.call(cbind, lapply(m_cov_list, `[[`, "filter"))

  bsseq_obj = BSseq(M = as.matrix(m), Cov = as.matrix(cov),
                    Filtered = as.matrix(filter),
                    coef = NULL,se.coef = NULL,
                    pos = start(overlap_gr),trans = NULL,
                    parameters = NULL, pData = colData, gr = NULL,
                    chr = as.vector(seqnames(overlap_gr)),
                    sampleNames = sampleNames,rmZeroCov = rmZeroCov)

  if (strandCollapse){bsseq:::strandCollapse(bsseq_obj)}
  return(bsseq_obj)
}
