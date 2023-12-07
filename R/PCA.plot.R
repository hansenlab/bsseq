PCA <- function(                                         
  BSseq.obj,
  CpG,
  win_size,
  nudge_x,
  nudge_y,
  size
) {
  
  # Filtering number of CpGs with 0 coverage in all samples
  coverage <- getCoverage(BSseq.obj)
  keep <- which(rowSums(coverage) !=0)
  BSseq.obj.filtered <- BSseq.obj[keep,]
  
  # Tile the genome
  BSseq.obj.tiled <- tile_by_windows(bs = BSseq.obj.filtered, win_size = win_size)
  
  # Filter out CpG tiles based on cutoff of number of CpGs
  
  BSseq.obj.tiled.filtered <- BSseq.obj.tiled[rowSums(as.matrix(countOverlaps(BSseq.obj.tiled, BSseq.obj.filtered))) > CpG, ]
  
  # Get raw methylation values
  
  BSseq.obj.methylation <- getMeth(BSseq.obj.tiled.filtered, type = "raw")

  # Filter tiles with NA methylation values
  BSseq.obj.methylation.filtered <- BSseq.obj.methylation[rowSums(is.na(BSseq.obj.methylation)) == 0,]
  
  # Perform PCA
  pca_data <- prcomp(t(BSseq.obj.methylation.filtered))
  pca_data_perc <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
  df_pca_data <- data.frame(
    PC1 = pca_data$x[, 1],
    PC2 = pca_data$x[, 2])
  
  pdf("PCA.pdf")  
  
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
