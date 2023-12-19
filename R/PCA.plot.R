PCA <- function(                                         
  BSseq.obj,
  genome,
  tilewidth,
  cut.last.tile.in.chrom,
  CpG,
  nudge_x,
  nudge_y,
  size
)

{
  
  # Tile the UCSC mm10 mouse genome
  
  mm10.mouse.genome.1kb.tiles <- GenomicRanges::tileGenome(seqinfo(genome), tilewidth = tilewidth, cut.last.tile.in.chrom = cut.last.tile.in.chrom)
  
  # Filtering number of CpGs with 0 coverage across all samples
  
  coverage <- getCoverage(BSseq.obj)
  keep <- which(rowSums(coverage) !=0)
  BSseq.obj.filtered <- BSseq.obj[keep,]

  # Find which 1kb tile has more than a number of CpGs (default is 3)
  
  one.kb.tiled.genome.more.number.of.CpGs <- mm10.mouse.genome.1kb.tiles[which(GenomicRanges::countOverlaps(mm10.mouse.genome.1kb.tiles, BSseq.obj.filtered) >= CpG), ]

  # Find which CpG in BSseq.object overlap with 1kb tiled coordinates with more than number of CpGs and then keep those CpGs in the BSseq.obj
  
  BSseq.obj.CpGs.within.1kb.tilled.genome <- BSseq.obj.filtered[as.data.frame(GenomicRanges::findOverlaps(BSseq.obj.filtered, one.kb.tiled.genome.more.number.of.CpGs))$queryHits, ]
  
  # Get raw methylation values
  
  BSseq.obj.methylation <- getMeth(BSseq.obj.CpGs.within.1kb.tilled.genome, type = "raw")
  
  # Filter tiles with NA methylation values
  
  BSseq.obj.methylation.filtered <- BSseq.obj.methylation[rowSums(is.na(BSseq.obj.methylation)) == 0, ]
  
  # Perform PCA
  pca_data <- prcomp(t(BSseq.obj.methylation.filtered))
  pca_data_perc <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
  df_pca_data <- data.frame(
    PC1 = pca_data$x[, 1],
    PC2 = pca_data$x[, 2])
  
  pdf("PCA.function.2.pdf")  
  
  # Plot PCA
  p <- ggplot(df_pca_data, aes(PC1, PC2, label = row.names(df_pca_data))) +
    geom_point() +
    labs(
      x = paste0("PC1 (", pca_data_perc[1], ")"),
      y = paste0("PC2 (", pca_data_perc[2], ")")
    ) +
    geom_text(nudge_x = nudge_x, nudge_y = nudge_y, size = size) +
    theme_classic()
  
  print(p)
  
  dev.off()

}
