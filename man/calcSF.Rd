\name{calcSF}
\alias{calcSF}
\title{calculates Scaling Factor (SF) values}
\description{
This internally called function calculates the scaling factor (SF) values for Affymetrix microarrays, for use in computing GCS-score values
}
\usage{
calcSF(diff, probetab, trim, clean.chip)
}
\arguments{
  \item{diff}{ The GC-content background corrected probe groupings for every probesetID or transcriptionclusterID on the given array type.  This is generated internally by the \code{computeSscore} function }
  \item{probetab}{ The internal datafile that contains the probe groupings and annotations for each array type and method type }
  \item{trim}{ The internal setting for the trimmed mean of every probe grouping on the array, as used in the calculation of SF.  For 3' IVT arrays, the \code{trim} is set to 0.02 by default. For all newer WT-type arrays, the \code{trim} is set to 0.04 by default }
  \item{clean.chip}{The clean chiptype name, based on the platform design package name.}
  }
\value{\code{calcSF} returns a numeric SF value for a given CEL file
  }
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
calcSF(diff, probetab, trim, clean.chip)
}
}