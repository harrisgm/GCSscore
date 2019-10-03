# Function to generate 'probeFile' probe-level annotations from platform 
# design (pd) packages from BioConductor.  The generated probeFile is stored 
# in a data.table within an .rda file. These 'probeFiles' are then converted 
# to a custom 'probePackage' using the makeProbePacage() function from
# the AnnotationForge package.

# This function is for use on WT gene/exon Affymetrix arrays.



# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/
# chip.pd <- "pd.mogene.1.0.st.v1"
# chip.pd <- "pd.moex.1.0.st.v1"
# example:
# probeFile.moex10stv1 <- WTgeneexonpFBuilder(chip.pd = "pd.moex.1.0.st.v1")

# Function: 
WTgeneexonpFBuilder <- function(chip.pd = NULL,species.pd = NULL) {
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
  
  # Create master probeFile for chip.pd:
  probeFile <- chip.pmfeature
  
  # key the probeFile to sort by fid (or other column type):
  setkey(probeFile,fid)
  
  # Remove unnecessary columns, if present:
  if (!is.na(match("atom", names(probeFile)))){
    probeFile[,c("atom") := NULL]
  }
  # If "row_names" column exists (such as in the HTA_2-0), remove it:
  if (!is.na(match("row_names", names(probeFile)))){
    probeFile[,c("row_names") := NULL]
  }
  
  # Remove probe sequences from geneST/exonST arrays to minimize file size:
  probeFile <- probeFile[,sequence := NULL]
  
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
    importfun = "getWTgeneexonprobefileData",
    check = FALSE)
  
  pkg.loc <- paste(outdir,"/",chip,".probeFile",sep="")
  
  # install the package, in the tempdir(), to the R library using 'devtools':
  devtools::install(pkg = pkg.loc)
  message(paste("GCSscore 'probeFile' created from BioConductor platform design (pd) package: ",chip.pd,sep=""))
  message(paste("GCSscore 'probeFile' package installed for chip: ",chip,sep=""))
  return(message("DONE"))	
}
