\name{GCSscore}
\alias{GCSscore}
\title{Main GCS-score Function}
\description{The main function used to call and run the GCS-score algorithm.}
\usage{
GCSscore(celFile1 = NULL, celFile2 = NULL, celTable = NULL,celTab.names = FALSE,
typeFilter = 0, method = 1, rm.outmask = FALSE, SF1 = NULL, SF2 = NULL, fileout = FALSE, 
gzip = FALSE, verbose = FALSE)
}
\arguments{
  \item{celFile1}{ If a one comparison run is desired, enter the filename and path to the 1st Affymetrix CEL file }
  \item{celFile2}{ If a one comparison run is desired, enter the filename and path to the 2nd Affymetrix CEL file }
  \item{celTable}{ If a batch run is desired, enter the filename and path to the CSV file containing the batch information }
  \item{celTab.names}{ If set to \code{TRUE}, then the GCS-score batch output is assigned the user-designated name (specified in the first column of the \code{celTable} CSV file (see examples)) }
  \item{typeFilter}{ If set to \code{0}, all available probe types are included in the calculation and normalization of the GCS-score values. If set to \code{1}, only probes  well-annotated probe_ids (from BioConductor \code{.db} packages) are included in the calculation and normalization of the GCS-score output }
  \item{method}{ This determines the method used to group and tally the probes_ids when calculating GCS-scores.  For Whole Transcriptome arrays, for gene-level (transcript_cluster_id-based) analysis, set \code{method = 1}, and for exon-level (probeset_id-based) analysis, set \code{method = 2}.  For the older generation arrays (3' IVT-style), if a PM-MM based background correction is desired, set \code{method = 1} (PM-MM gives identical results to the original \code{sscore} package).  If a GC-content based background correction is desired on the 3' IVT arrays, set \code{method = 2} }
  \item{rm.outmask}{ If set to \code{TRUE}, then probes that are flagged as MASKED or OUTLIER in either CEL file 1 or CEL file 2 will be removed from the analysis.  If set to \code{FALSE}, these probes are not filtered out and will be used in the GCS-score calculation }
  \item{SF1}{ Input a pre-determined Scaling Factor (SF) for the 1st CEL file }
  \item{SF2}{ Input a pre-determined Scaling Factor (SF) for the 2nd CEL file }
  \item{fileout}{ Determines if the resulting GCS-score output is written to disk in a CSV format following the completion of the
  function. By default, this is set to FALSE so unnecessary GCS-score outputs are not saved to disk after each run }
  \item{gzip}{ If set to \code{TRUE}, the \code{GCSscore} output that is written to disk is compressed.  This could prove useful if a large number runs are input using the batch submission }
  \item{verbose}{If set to \code{TRUE}, more information will be printed to the console during while the algorithm is running }
}
\details{
The input accepts individual CEL files or reads in a CSV file for batch runs.  The user also inputs parameters to determine the \code{method} used by the GCS-score algorithm to group and tally the individual probes on a given array.
}
\value{An \code{ExpressionSet} object with GCS-score values for the probe groupings (determined by the \code{method} argument) and the relevant annotation informtaion }
\examples{
if (length(list.files(path = ".", pattern = "*.CEL")) != 0){

######################## Single run example ###########################

# get the path to example CEL files provided with package:
celpath1 <- system.file("extdata/","MN_2_3.CEL", package = "GCSscore")
celpath2 <- system.file("extdata/","MN_4_1.CEL", package = "GCSscore")

# run GCSscore() function directly on the two .CEL files above:
GCSs.single <- GCSscore(celFile1 = celpath1, celFile2 = celpath2)

# convert GCSscore single-run from ExpressionSet to data.table:
GCSs.single.dt <- 
  data.table::as.data.table(cbind(GCSs.single@featureData@data,
                                  GCSs.single@assayData[["exprs"]]))

# show all column names included in the output:
colnames(GCSs.single.dt)

# show simplified output of select columns and rows:
GCSs.single.dt[10000:10005,
              c("transcriptclusterid","symbol",
                "ref_id","Sscore")]

######################## batch run example ############################

# get the path to example batch (.csv) file provided with package:
celtab_path <- system.file("extdata",
                           "GCSs_batch_ex.csv", 
                           package = "GCSscore")

# read in the .CSV file using fread():
celtab <- data.table::fread(celtab_path)

# view structure of 'celTable' input:
celtab

# add the path to the sample CEL files to the batch input:
#   NOTE: this step is not necessary if the .CEL files 
#         are in the working directory:

path <- system.file("extdata", package = "GCSscore")
celtab$CelFile1 <- celtab[,paste(path,CelFile1,sep="/")]
celtab$CelFile2 <- celtab[,paste(path,CelFile2,sep="/")]

# run GCSscore function on the batch input:
GCSs.batch <- GCSscore::GCSscore(celTable = celtab, celTab.names = TRUE)

# convert GCS-score output from 'ExpressionSet' to 'data.table':
GCSs.batch.dt <-
  data.table::as.data.table(cbind(GCSs.batch@featureData@data,
                                  GCSs.batch@assayData[["exprs"]]))
                                  
# show all columns included in the output:
colnames(GCSs.batch.dt)

# show simplified output of GCSscore batch example:
GCSs.batch.dt[10000:10005,
              c("transcriptclusterid","symbol",
                "example01","example02","example03")]
}
}