\name{plotRegion}
\alias{plotRegion}
\alias{plotManyRegions}
\title{
  Plotting BSmooth methylation estimates
}
\description{
  Functions for plotting BSmooth methylation estimates.  Typically used
  to display differentially methylated regions.
}
\usage{
plotRegion(BSseq, region = NULL, extend = 0, main = "",
  addRegions = NULL, annoTrack = NULL, cex.anno = 1,
  geneTrack = NULL, cex.gene = 1.5, col = NULL, lty = NULL,
  lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected",
  stat.col = "black", stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8, 8),
  mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
  addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE)

plotManyRegions(BSseq, regions = NULL, extend = 0, main = "",
  addRegions = NULL, annoTrack = NULL, cex.anno = 1,
  geneTrack = NULL, cex.gene = 1.5, col = NULL, lty = NULL,
  lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected",
  stat.col = "black", stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8, 8),
  mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
  addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE,
  verbose = TRUE)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{region}{A \code{data.frame} (with start, end and chr columns)
    with 1 row or \code{GRanges} of length 1.  If \code{region} is
    \code{NULL} the entire \code{BSseq} argument is plotted.}
  \item{regions}{A \code{data.frame} (with start, end and chr columns)
    or \code{GRanges}.}
  \item{extend}{Describes how much the plotting region should be
  extended in either direction.  The total width of the plot is equal to
  the width of the region plus twice \code{extend}.}
  \item{main}{The plot title.  The default is to construct a title with
  information about which genomic region is being plotted.}
  \item{addRegions}{A set of additional regions to be highlighted on the
  plots.  As the \code{regions} argument.}
  \item{annoTrack}{A named list of \code{GRanges} objects.  Each
  component is a track and the names of the list are the track names.
  Each track will be plotted as solid bars, and we routinely display
  information such as CpG islands, exons, etc.}
  \item{cex.anno}{\code{cex} argument when plotting \code{annoTrack}.}
  \item{geneTrack}{\strong{EXPERIMENTAL}: A \code{data.frame} with columns:
  \code{chr}, \code{start}, \code{end}, \code{gene_ID}, \code{exon_number},
  \code{strand}, \code{gene_name}, \code{isoforms}. This interface is under
  active development and subject to change.}
  \item{cex.gene}{\code{cex} argument when plotting \code{geneTrack}.}
  \item{col}{The color of the methylation estimates, see details.}
  \item{lty}{The line type of the methylation estimates, see details.}
  \item{lwd}{The line width of the methylation estimates, see details.}
  \item{BSseqStat}{An object of class \code{BSseqStat}.  If present,
    a new panel will be shown with the t-statistics.}
  \item{stat}{Which statistics will be plotted (only used is
    \code{BSseqStat} is not \code{NULL}.)}
  \item{stat.col}{color for the statistics plot.}
  \item{stat.lwd}{line width for the statistics plot.}
  \item{stat.lty}{line type for the statistics plot.}
  \item{stat.ylim}{y-limits for the statistics plot.}
  \item{mainWithWidth}{Should the default title include information
  about width of the plot region.}
  \item{regionCol}{The color used for highlighting the region.}
  \item{addTicks}{Should tick marks showing the location of methylation
  loci, be added?}
  \item{addPoints}{Should the individual unsmoothed methylation
  estimates be plotted.  This usually leads to a very confusing plot,
  but may be useful for diagnostic purposes.}
  \item{pointsMinCov}{The minimum coverage a methylation loci need in
  order for the raw methylation estimates to be plotted.  Useful for
  filtering out low coverage loci.  Only used if \code{addPoints =
  TRUE}.}
  \item{highlightMain}{Should the plot region be highlighted?}
  \item{verbose}{Should the function be verbose?}
}
\details{
  The correct choice of aspect ratio depends on the width of the
  plotting region.  We tend to use \code{width = 10, height = 5}.

  \code{plotManyRegions} is used to plot many regions (hundreds or
  thousands), and is substantially quicker than repeated calls to
  \code{plotRegion}.

  This function has grown to be rather complicated over time.  For
  custom plotting, it is sometimes useful to use the function
  definition as a skeleton and directly modify the code.
}
\value{This function is invoked for its side effect: producing a plot.}
\seealso{
  The package vignette has an extended example.
}
\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}
