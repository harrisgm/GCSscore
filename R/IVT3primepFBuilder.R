# Function to generate 'probeFile' probe-level annotations from platform 
# design (pd) packages from BioConductor.  The generated probeFile is stored 
# in a data.table within an .rda file. These 'probeFiles' are then converted 
# to a custom 'probePackage' using the makeProbePacage() function from
# the AnnotationForge package.

# This function is for use on 3' IVT Affymetrix arrays.

# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/

# Example:
# probeFile.mo4302 <- ClariomSpFBuilder(chip.pd = "pd.mouse430.2")

# Function: 
IVT3primepFBuilder <- function(chip.pd = NULL, clean.chip = NULL, species.pd = NULL, pF.type = NULL) {
  # Install necessary pd.* package if not already installed:
  if (!requireNamespace(chip.pd, quietly = TRUE)){
    BiocManager::install(chip.pd)
    message(paste("installing necessary BioConductor platform design (pd) package: ",chip.pd,sep=""))}
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
  # Close the connections now that the data from the SQL file has been read:
  dbDisconnect(conn=con)
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
  chip.pmSeq <- chip.pmSeq[!duplicated(fid)]
  setkey(chip.pmSeq,fid)
  # merge data by key:
  chip.pmfeature <- chip.pmSeq[chip.pmfeature]
  # add in GC.count info for each sequence in the pmfeature object:
  chip.pmfeature[,GC.count := str_count(sequence, "G|C")]
  # NOTE: mmSequence is empty:
  # data("mmSequence", package = chip.pd, envir = environment())
  # chip.mmSeq <- as.data.table(mmSequence)
  
  probeFile <- chip.pmfeature
  # key the probeFile to sort by fid (or other column type):
  setkey(probeFile,fid)
  
  # Remove unnecessary 'atom' column, if present:
  if (!is.na(match("atom", names(probeFile)))){
    probeFile[,c("atom") := NULL]
  }
  # 03.08.20: Go ahead and remove the "sequence" column:
  # There is no point in including it in the probeFile:
  # Remove probe sequences from ClariomS-style arrays to minimize file size:
  probeFile <- probeFile[,sequence := NULL]
  
  
  # GMH: set the intermediate files to output to write to a temporary directory:
  outdir <- tempdir()
  message("writing intermediary .probe_tab file to temporary directory")
  probe.tab <-  paste("GCSs.",chip,".probeFile.probe_tab",sep="")
  probe.tab.loc <- paste(outdir,"/",probe.tab,sep="")
  fwrite(file = probe.tab.loc,probeFile,sep = "\t")
  probe.pkg.name <- paste(clean.chip,".probeFile",sep="")
  
  
  makeProbePackageGCSs(
    arraytype = clean.chip,
    outdir = outdir,
    species = species.pd,
    chip.pd = chip.pd,
    maintainer= "Guy Harris <harrisgm@vcu.edu>",
    version = "0.0.5",
    pkgname = probe.pkg.name,
    datafile = probe.tab.loc,
    pF.type = pF.type,
    importfun = "get3primeIVTprobefileData",
    check = FALSE)
  
  pkg.loc <- paste(outdir,"/",probe.pkg.name,sep="")
  # install the package, in the tempdir(), to the R library using 'devtools':
  devtools::install(pkg = pkg.loc,upgrade = "never")
  message(paste("GCSscore 'probeFile' created from BioConductor platform design (pd) package: ",chip.pd,sep=""))
  message(paste("GCSscore 'probeFile' package installed for chip: ",chip,sep=""))
  return(message("DONE"))	
}
