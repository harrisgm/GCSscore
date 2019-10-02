# Function to generate 'probeFile' probe-level annotations from platform design (pd) 
# packages from BioConductor.  The generated 'probeFile' data is as a data.table in 
# an .rda file.

# This function is for use on ClariomD/XTA Affymetrix arrays.

# These 'probeFiles' are then converted to a custom 'probePackage' using AnnotationForge.
# The necessary functions for data parsing / annotation package creation can be found on Github:
# https://github.com/harrisgm/GCSscore-probeFile-functions


# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/

# probeFile.hta20 <- ClariomDXTApFBuilder(chip.pd = "pd.hta.2.0")
# probeFile.mta10 <- ClariomDXTApFBuilder(chip.pd = "pd.mta.1.0")
# probeFile.rta10 <- ClariomDXTApFBuilder(chip.pd = "pd.rta.1.0")
# probeFile.clariomdhuman <- ClariomDXTApFBuilder(chip.pd = "pd.clariom.d.human")

# # load necessary libaries:
# library(data.table)
# # library(devtools)
#   # NOTE: devtools not needed if not installing a package from GitHub
# library(stringr)
# library(RSQLite)

# Function:
ClariomDXTApFBuilder <- function(chip.pd = NULL, species.pd = NULL) {
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
  # Remove unnecessary 'atom' column, if present:
  if (!is.na(match("atom", names(chip.pmfeature)))){
    chip.pmfeature[,c("atom") := NULL]
  }
  
  # NOTE: chip.probefile.BioC doesn't contain the named "features", it ONLY has the numbered ones
  # NOTE: chip.core.mps has TCid names in it!
  chip.core.mps <- as.data.table(chip.pdinfo[["core_mps"]])
  # GMH 9.23.19: removing 'transcript_cluster_id' from core.mps:
  #REASON: the 'transcript_cluster_id' is incomplete (for JUC probesets) relative to the 'featureSet' object:
  if (!is.na(match("transcript_cluster_id", names(chip.core.mps)))){
    chip.core.mps[,c("transcript_cluster_id") := NULL]
  }
  # Additionally, delete 'row.names' (in hta20) column right now, rather than at the end.
  if (!is.na(match("row_names", names(chip.core.mps)))){
    chip.core.mps[,c("row_names") := NULL]
  }
  # Add PSR/JUC ids here as well:
  chip.featureSet <- as.data.table(chip.pdinfo[["featureSet"]])
  # Cut featureSet down to just 2 columns:
  # GMH 9.23.19: adding transcript_cluster_id back to featureSet (back to 3 columns):
      #REASON: was losing the TCid info for JUC probesets  
  chip.featureSet <- chip.featureSet[,c("fsetid","man_fsetid","transcript_cluster_id")]
  
  # Remove the additional labeling (present in some of the XTA featureSets):
  # Done on a conditional basis due to subtle variaitons.  User may need to customize these commands to fit future array releases:
  # Users may need to add additional array types to the list in the if statement:
  
  # for hta20:remove extra ID repeats, add '.1' back to the end of the transcript_cluster_ids:
  if (chip.pd == "pd.hta.2.0"){
    chip.featureSet[,transcript_cluster_id := sapply(strsplit(transcript_cluster_id, "///"), "[", 1 )]
    # add '.1' to the end of singly-parsed TCid:
    chip.featureSet[,transcript_cluster_id := sapply(paste(transcript_cluster_id,".1",sep = ""),"[",1)]
    # change 'NA.1' values created by the above command back to regular NA values:
    chip.featureSet[transcript_cluster_id=="NA.1","transcript_cluster_id"] <- NA 
  }
  
  
  setkey(chip.pmfeature,fsetid)
  setkey(chip.core.mps,fsetid)
  setkey(chip.featureSet,fsetid)
  
  # Add in info from core.mps and featureSet:
  chip.pmfeature <- chip.core.mps[chip.pmfeature]
  chip.pmfeature <- chip.featureSet[chip.pmfeature]
  # setkey of 'pmfeature' to 'fid' for matching probes with the corresponding 'sequence':
  setkey(chip.pmfeature,fid)

  # Get PM sequence information from chip.pd: ----------------------------
  
  data("pmSequence", package = chip.pd, envir = environment())
  chip.pmSeq <- as.data.table(pmSequence)
  # Necessary for hta2.0: must remove duplicated 'fid' rows or merge will add duplicated rows!!
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
  affynet.TC <- as.data.table(netaffxTranscript@data)
  load(system.file("extdata","netaffxProbeset.rda",package = chip.pd))
  affynet.PSR <- as.data.table(netaffxProbeset@data)

  # if (chip %like% "mta10"| chip %like% "hta20" | chip %like% "rta10" | chip %like% "clariomdhuman"){
    affynet.TC <- affynet.TC[,c("transcriptclusterid","category")]
    setkey(affynet.TC,transcriptclusterid)
    setkey(chip.pmfeature,transcript_cluster_id)
    chip.pmfeature <- affynet.TC[chip.pmfeature]
# }
    
    affynet.PSR <- affynet.PSR[,c("probesetid","locustype","probesettype")]
    # table(table(affynet.PSR$locustype))
    # NOTE: for hta2.0, many rows affynet.PSR has multiple locustype entries per row:
    # Remove all but first 'locustype' and 'transcriptclusterid' entry:
    affynet.PSR[,locustype := sapply(strsplit(affynet.PSR$locustype, "///"), "[", 1 )]
    # affynet.PSR[,transcriptclusterid := sapply(strsplit(affynet.PSR$transcriptclusterid, "///"), "[", 1 )]
    # Assign key of affynet.PSR to 'probesetid' for PSR/JUC 'locustype' and 'probesettype' merge:
    setkey(affynet.PSR,probesetid)
    # Change key of pmfeature to 'man_fsetid' (PSR/JUC) for 'locustype' and 'probesettype' merge:
    setkey(chip.pmfeature,man_fsetid)
    # to test, run:
        # View(chip.pmfeature[man_fsetid %!in% affynet.PSR$probesetid])
    # Create master probeFile for using chip.pd:
    probeFile <- affynet.PSR[chip.pmfeature]
  
  # key the probeFile to sort by fid (or other column type):
  setkey(probeFile,fid)
  
  # Remove probe sequences from XTA-style arrays to minimize file size:
  probeFile <- probeFile[,sequence := NULL]
  
  # GMH: set the intermediate files to output to write to a temporary directory:
  outdir <- tempdir()
  message("writing intermediary .probe_tab file to temporary directory")
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
    importfun = "getXTAprobefileData",
    check = FALSE)
  
  pkg.loc <- paste(outdir,"/",chip,".probeFile",sep="")
  
  # install the package, in the tempdir(), to the R library using 'devtools':
  devtools::install(pkg = pkg.loc)
  message(paste("GCSscore 'probeFile' created from BioConductor platform design (pd) package: ",chip.pd,sep=""))
  message(paste("GCSscore 'probeFile' package installed for chip: ",chip,sep=""))
  return(message("DONE"))	
}
