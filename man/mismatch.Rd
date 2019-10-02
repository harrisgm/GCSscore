\name{mismatch}
\alias{mismatch}
\title{Calculates mismatch values for probes}
\description{ This internally called function calculates the background correction for each probe based on the median intensity of all background probes with the same GC-content of as the target probe in question 
}
\usage{
mismatch(probes, bgp, intensity)
}
\arguments{
  \item{probes}{ probe indicies for target probes. Each probe index contains the gc-content of the probe }
  \item{bgp}{ contains probe location (indicies), GC-content, and annotations of the background probes of a given chip type.  For WT-type arrays, the \code{bgp} consists of 16,943 antigenomic background probes.  For 3' IVT arrays, the mismatch (MM) probes are used to calculate the \code{bgp} list in both methods  }
  \item{intensity}{ The intensities value of the \code{bgp} probes as read in from the .CEL file }
}
\details{
This internally called function calculates the probe background correction based on the median intensity of all background probes with the same gc-content of as the target probe in question 
}
\value{
\code{mismatch} returns a numeric vector containing the gc-content based background correction for every probe included in the analysis }
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
mismatch(probes, bgp, intensity)
}
}