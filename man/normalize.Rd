\name{normalize}
\alias{normalize}
\title{Normalization of GCS-score values}
\description{
Normalizes the GCS-score values using all scores within 3*SD of the mean. This normalization step occurs after the probes have been tallied and grouped into probe_ids, according to the method (probeset_id for exon-level or transcription_cluster_id for gene-level }
\usage{
normalize(Score)
}
\arguments{
  \item{Score}{ The unnormalized GCS-score values (grouped and tallied according to the \code{method} selection) that are generated in the \code{computeScore} function }
  }
\value{
\code{normalize} Returns a numeric vector containing normalized GCS-score values for every probe_id grouping included in the analysis 
}
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

#Example of input, as the function would be called internally:
normalize(Score)
}
}