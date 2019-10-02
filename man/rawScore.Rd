\name{rawScore}
\alias{rawScore}
\title{Calculates the rawScore values}
\description{
Calculates rawScore values based on differences between the two background corrected arrays in a given GCS-score analysis (e.g. CEL_1 vs. CEL_2), using the internally generated Statistical Difference Threshold (SDT) values.
}
\usage{
rawScore(diff1, diff2, SDT1, SDT2)
}
\arguments{
  \item{diff1}{ The gc-background-corrected values for the probe intensities on the 1st array }
  \item{diff2}{ The gc-background-corrected values for the probe intensities on the 2nd array }
  \item{SDT1}{ The internally calculated Statistical Difference Threshold (SDT=4*rawQ*SF) for the 1st array }
  \item{SDT2}{ The internally calculated Statistical Difference Threshold (SDT=4*rawQ*SF) for the 2nd array }
  }
\value{
\code{rawScore} returns a numeric vector containing the raw, ungrouped scores for every probe grouping included in the analysis (as determined by \code{method}) }
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
rawScore(diff1, diff2, SDT1, SDT2)
}
}