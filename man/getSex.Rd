\name{getSex}
\alias{getSex}
\alias{addSex}
\alias{plotSex}
\title{
  Estimating sample sex based on sex chromosome coverage
}
\description{
  Estimates sample sex based on sex chromosome coverage
}
\usage{
getSex(BSseq, cutoff = 0)
addSex(BSseq, sex = NULL)
plotSex(df, id = NULL)
}
\arguments{
  \item{BSseq}{An object of class \code{BSseq}.}
  \item{cutoff}{Value used to distinguish males and females when evaluating
  difference between coverage of X and Y chromosomes.}
  \item{sex}{An optional character vector of sex (with values \code{M}
    and \code{F}.}
  \item{df}{A data frame structured like the output of \code{getSex},
	with columns \code{predictedSex} \code{xCov} and \code{yCov}}
  \item{id}{Text used as plotting symbols in the \code{plotSex}
    function.  Used for sample identification on the plot.}
}
\details{
	Estimation of sex is based on coverage of the X and Y chromosomes.
	If a sample exhibits high coverage of the X chromosome and low (no)
	coverage of the Y chromosome, we predict female sex. In the opposite
	case, we predict male sex. 
	
	This function should be run on a \code{BSseq} object prior to any CpG-based
	filtering.
}
\value{
  For \code{getSex}, a \code{DataFrame} with columns \code{predictedSex}
  (a character with values \code{M} and \code{F}), \code{xCov} and
  \code{yCov} which are the average coverage values of the X and Y chromosomes,
  respectively.

  For \code{addSex}, an object of the same type as \code{BSseq} but
  with a column named \code{predictedSex} added to the pheno data.
}
\author{
  Shan V. Andrews \email{sandre17@jhu.edu}
}
\examples{
\dontrun{
#Not applicable to test data (only chr 22)
sex<-getSex(BSseq)
plotSex(sex)
}
}
