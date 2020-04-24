\name{computeSscore}
\alias{computeSscore}
\title{Computes GCS-score values}
\description{
This internally called function computes the GCS-score values between two Affymetrix-style microarrays. The \code{computeSscore} function contains the majority of the GCS-score algorithm.
}
\usage{
computeSscore(cel1, cel2, probeFile, bgp, method, infoKey, SF1 = NULL, SF2 = NULL, 
              verbose = FALSE, trim = NULL, clean.chip)
}
\arguments{
  \item{cel1}{ The 1st Affymetrix CEL file, as read in by the \code{affxparser} package }
  \item{cel2}{ The 2nd Affymetrix CEL file, as read in by the \code{affxparser} package }
  \item{probeFile}{ The internal datafile that contains the probe groupings and annotations for each array type and method type }
  \item{bgp}{ The index of the probe location, GC-content, and annotations of the background probes of a given chip type.  For WT-type arrays, the \code{bgp} consists of 16,943 antigenomic background probes.  For 3' IVT arrays, the MisMatch (MM) probes are used to calculate the \code{bgp} list in both methods  }
  \item{method}{ Determines the method used to group and tally the probes when calculating GCS-score values }
    \item{infoKey}{ The \code{key} of how to group the probes together for the GCS-score calculations. Determines the method used to group and tally the probes when calculating GCS-score values.  For example, exon-level analysis groups probes by probeset_ids while gene-level groups probes by transcript_cluster_ids }
  \item{SF1}{ If the user has predetermined scaling factors, input user Scaling Factor (SF) for the 1st CEL file.  Otherwise, the \code{computeSscore} function will caluclate \code{SF1} directly from the 1st CEL file }
  \item{SF2}{ If the user has predetermined scaling factors, input user Scaling Factor (SF) for the 2nd CEL file.  Otherwise, the \code{computeSscore} function will caluclate \code{SF2} directly from the 2nd CEL file }
  \item{verbose}{If set to \code{TRUE}, additional information will be printed to the console while the algorithm is running.}
  \item{trim}{ Internal parameter determined by chip type .\code{trim}=0.04 for WT-type arrays and 0.02 for 3' IVT type arrays }
  \item{clean.chip}{The clean chiptype name, based on the platform design package name.}
  }
\details{
This internally called function computes the raw difference scores between the probes on each microarray, then groups the probes into probesets or transcript cluset ids, and normalizes the results to produce GCS-score values. The function returns the values to the main \code{GCscore2}, where BioConductor-based annotations are added to either the exon-level or gene-level probe groupings 
}
\value{A \code{data.table} object with GCS-Score values for the probe groupings (determined by the \code{method} argument)  }
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
computeSscore(cel1, cel2, probeFile, bgp, infoKey, method, SF1 = NULL, SF2 = NULL,
              verbose = FALSE, trim = NULL)
}
}