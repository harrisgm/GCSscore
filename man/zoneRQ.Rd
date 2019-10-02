\name{zoneRQ}
\alias{zoneRQ}
\title{Calculates the zone-based RawQ values}
\description{ This internally called function calculates the zone-based RawQ values. RawQ is a measure of the noise within a given zone on a microarray chip.  This noise value is used in the error model contained in the GCS-score algorithm
}
\usage{
zoneRQ(DT, affyCel, trim)
}
\arguments{
  \item{DT}{ Internally generated \code{data.table} containing the .CEL data, generated from the \code{list} that is created by the \code{affxparser} package }
  \item{affyCel}{ The .CEL file data, in \code{list} structure, as read in using the \code{readCel} function included in the \code{affxparser} package  }
  \item{trim}{ The internal setting for the trimmed mean of every probe grouping on the array.  For 3' IVT arrays, the \code{trim} is set to 0.02 by default. For all newer WT-type arrays, the \code{trim} is set to 0.04 by default }
  }
\value{
\code{zoneRQ} returns a numeric vector containing zone-based rawQ values for a given array 
}
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
zoneRQ(DT, affyCel, trim)
}
}