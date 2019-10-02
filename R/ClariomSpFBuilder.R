# Function to generate 'probeFile' probe-level annotations from platform design (pd) 
# packages from BioConductor.  The generated 'probeFile' data is as a data.table in 
# an .rda file.

# This function is for use on ClariomS Affymetrix arrays.

# These 'probeFiles' are then converted to a custom 'probePackage' using AnnotationForge.
# The necessary functions for data parsing / annotation package creation can be found on Github:
# https://github.com/harrisgm/GCSscore-probeFile-functions


# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/

# chip.pd <- "pd.clariom.s.mouse"
# chip.pd <- "pd.clariom.s.human"
# chip.df <- "pd.clariom.s.rat"
# probeFile.clariomsrat <- ClariomSpFBuilder(chip.pd = "pd.clariom.s.rat")

# # load necessary libaries:
# library(data.table)
# # library(devtools)
#   # NOTE: devtools not needed if not installing a package from GitHub
# library(stringr)
# library(RSQLite)

# Function: 
ClariomSpFBuilder <- function(chip.pd = NULL, species.pd = NULL) {
  # Install necessary pd.* package if not already installed:
  if (!requireNamespace(chip.pd, quietly = TRUE)){
    BiocManager::install(chip.pd)
    message("installing necessary platform design (pd) package from Bioconductor")}
  
  # get chip name (without periods/dots) from chip.pd:
  chip <- str_remove_all(strsplit(chip.pd,"pd.")[[1]][2],"[.]")
  
  # Load .sqlite file from pd.mta.1.0 package:
  chip.sqlite <- paste(chip.pd,".sqlite",sep ="")
  chip.db <- system.file("extdata",chip.sqlite,package = chip.pd)
  ## connect to db
  con <- dbConnect(drv=RSQLite::SQLite(), dbname=chip.db)
  ## list all tables
  tables <- dbListTables(con)
  ## exclude sqlite_sequence (contains table information)
  tables <- tables[tables != "sqlite_sequence"]
  chip.pdinfo <- vector("list", length=length(tables))
  names(chip.pdinfo) <- tables
  ## create a data.frame for each table
  for (i in seq(along=tables)) {
    chip.pdinfo[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
  }
  # Extract probe-level data:
  chip.pmfeature <- as.data.table(chip.pdinfo[["pmfeature"]])
  chip.featureSet <- as.data.table(chip.pdinfo[["featureSet"]])
  setkey(chip.pmfeature,fsetid)
  setkey(chip.featureSet,fsetid)
  chip.pmfeature <- chip.featureSet[chip.pmfeature]
  setkey(chip.pmfeature,fid)
  
  # Get PM sequence information from chip.pd: ----------------------------
  
  data("pmSequence", package = chip.pd, envir = environment())
  chip.pmSeq <- as.data.table(pmSequence)
  # must remove duplicated 'fid' rows or merge will add duplicated rows!!
  # NOTE:  unique vs. !duplicated --> ALWAYS use (!duplicated)
  # chip.test.unique<- chip.pmSeq[unique(fid)]
  # chip.test.notdup <- chip.pmSeq[!duplicated(fid)]
  # These objects have the SAME nrow, BUT only the chip.test.notdup will merge without extra rows!
  # Also the chip.unqiue is LARGER in megabytes --> implies NAs must have been introduced with the unique() command
  chip.pmSeq <- chip.pmSeq[!duplicated(fid)]
  setkey(chip.pmSeq,fid)
  
  chip.pmfeature <- chip.pmSeq[chip.pmfeature]
  # add in GC.count info for each sequence in the pmfeature object:
  chip.pmfeature[,GC.count := str_count(sequence, "G|C")]
  # Load in pd.chip.1.0 “netaffxTranscript" data for locustype/category---
  load(system.file("extdata","netaffxTranscript.rda",package = chip.pd))
  # creates "netaffxTranscript" object

  
  affynet.TC <- as.data.table(netaffxTranscript@data)
  # colnames(affynet.TC)
  # table(affynet.TC$locustype)
  # Coding NonCoding 
  # 44910     22848 
  
  # if (chip %like% "clarioms"){
    affynet.TC <- affynet.TC[,c("transcriptclusterid","category","locustype")]
    setkey(affynet.TC,transcriptclusterid)
    setkey(chip.pmfeature,man_fsetid)
    probeFile <- affynet.TC[chip.pmfeature]
  # }
  
  # key the probeFile to sort by fid (or other column type):
  setkey(probeFile,fid)
  
  # Remove unnecessary 'atom' column, if present:
  if (!is.na(match("atom", names(probeFile)))){
    probeFile[,c("atom") := NULL]
  }
  
  # GMH: set the intermediate files to output to write to a temporary directory:
  outdir <- tempdir()
  
  probe.tab.name <-  paste("GCSs.",chip,".probeFile.probe_tab",sep="")
  probe.tab.loc <- paste(outdir,probe.tab.name,sep="")
  fwrite(file = probe.tab.loc,probeFile,sep = "\t")
  arraytype <- chip
  # datafile <- "mouse4302.probeFile.probe_tab"
  # For structure of Data:
  # sapply(probeFile,class)
  
  # Running stock function from AnnotationForge package version 1.26.0
  makeProbePackage(
    arraytype = chip,
    outdir = outdir,
    species = species.pd,
    maintainer= "Guy Harris <harrisgm@vcu.edu>",
    version = "0.0.1",
    datafile = probe.tab.loc,
    importfun = "getClariomSprobefileData",
    check = FALSE)
  
  pkg.loc <- paste(outdir,"/",chip,".probeFile",sep="")
  
  # install the package, in the tempdir(), to the R library using 'devtools':
  devtools::install(pkg = pkg.loc)
  message(paste("GCSscore 'probeFile' created from BioConductor platform design (pd) package: ",chip.pd,sep=""))
  message(paste("GCSscore 'probeFile' package installed for chip: ",chip,sep=""))
  return(message("DONE"))	
}
