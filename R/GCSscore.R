# MAIN FUNCTION/WRAPPER FOR THE GCSscore ALGORITHM:
GCSscore <- function(celFile1 = NULL, celFile2 = NULL, celTable = NULL,
                     celTab.names = FALSE, typeFilter = 0, method = 1, 
                     rm.outmask = FALSE, SF1 = NULL, SF2 = NULL, 
                     fileout = FALSE, gzip = FALSE, verbose = FALSE) {
  
  # Version checks: Lastest versions of the GCSscore-built packages and BioC:
  probepkg.vers = "0.0.6"
  annot.vers = "0.0.5"
  bioC.latest = "3.10"
  # Define inverse operators:
  `%!in%` <- Negate(`%in%`)
  `%!like%` <- Negate(`%like%`)
  
  # # START OF TESTING VARIABLES
  # chip <- "RTA-1_0"
  # chip <- "MTA-1_0"
  # chip <- "Clariom_S_Mouse"
  # chip <- "Mouse430_2"
  # clean.chip <- tolower(gsub("-|_", "",chip))
  # clean.chip <- tolower(gsub("v1", "",clean.chip))
  # probepkg <- paste(clean.chip,".probeFile",sep="")
  # pdpkg <- paste("pd.",tolower(gsub("-|_", ".",chip)),sep = "")
  # chip.pd <- pdpkg
  # # END OF TESTING VARIABLES

# BioConductor Version checks ---------------------------------------------

  # Need to make sure Bioconductor is at least "3.10"
  # So, bioconductor has to be installed to download it from the repository,
  # But, the package will be on github as well.
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  } else if (BiocManager::version() < bioC.latest){
    # This command will update all Bioc packages to version 3.10
    message(paste("Please ugrade to the latest version of bioconductor: ",bioC.latest,sep=""))
    # NOTE: you should only update the packages that are necessary for the GCSscore to run:
    # I had 83 packages to update: most of them do not matter, and some require compilation!
    BiocManager::install(version = bioC.latest)
  } else {message("\n  **latest version of bioconductor is already installed**")}
  
  # NOW, the installation works like a charm!!
  # BiocManager::install("GCSscore", update = FALSE)
  

# CEL file support and logic checks ---------------------------------------

  # Basic testing that all entered parameters are set to usable values:
  # if (!(typeFilter == 0 | typeFilter == 1))stop("typeFilter must be (0) or (1). (1) is the recommended setting")
  # if (!(method == 1 | method == 2))stop("method must be (1) or (2)")
  # # check that either celFile1/2 or celTable is selected:
  # if (is.null(celTable)) stopifnot(!is.null(celFile1), !is.null(celFile2))
  # else if (is.null(celFile1) & is.null(celFile2)) stopifnot(!is.null(celTable))
  # else stop("input either: 2 .CEL files (single-run) or 'celTable' (batch-run)")
  
  # NEW
  # Basic testing that all entered parameters are set to usable values:
  if (!(typeFilter == 0 | typeFilter == 1))stop("typeFilter must be (0) or (1). (0) is the recommended setting")
  if (!(method == 1 | method == 2))stop("method must be (1) or (2)")
  # check that either celFile1/2 or celTable is selected:
  if (is.null(celTable)) {stopifnot(!is.null(celFile1), !is.null(celFile2))
  } else if (is.null(celFile1) & is.null(celFile2)) {stopifnot(!is.null(celTable))
  } else stop("input either: 2 .CEL files (single-run) or 'celTable' (batch-run)")
  # END NEW
  
  # DEFINE CHIP-TYPES THAT THE SSCORE2 CODE IS COMPATIABLE WITH:
  # All chipXTA and chipClarS have been checked (all chips as of 09.30.2019):
  # EDIT 02.17.20: added support for U133A_2 chip type
  # EDIT 03.03.20: removed support for "Rhesus chip type.
  # Focus is providing support for the 2 most widely used platforms for the 3' IVT. 

  chip3IVT <- c("Mouse430_2","HG-U133A_2")
  chipXTA <- c("MTA-1_0","HTA-2_0","RTA-1_0","Clariom_D_Human")
  chipClarS <- c("Clariom_S_Mouse","Clariom_S_Human","Clariom_S_Rat")
  # List of all supported chip types (as of 03.10.20):
  chipType <- c(chip3IVT, chipXTA, chipClarS)
  
  if (!is.null(celTable)) {
    if (!is.data.table(celTable)) celTable <- as.data.table(celTable)
    chipType1 <- chipType2 <- vector("character", length = nrow(celTable))
    for (i in 1:nrow(celTable)) {
      stopifnot (!is.null(celTable[[2]][i]), !is.null(celTable[[3]][i]))  
      chipType1[i] <- readCelHeader(celTable[[2]][i])$chiptype
      chipType2[i] <- readCelHeader(celTable[[3]][i])$chiptype
    }
    
    if (all.equal(chipType1, chipType2)){
      # add in 'celHead1' for batch runs of the 3'IVT chips:
      celHead1 <- readCelHeader(celTable[[2]][1])
      chip <- chipType1[1] 
      clean.chip <- tolower(gsub("-|_", "",chip))
      } else stop("Not all .CEL files in the 'batch.csv' are the same chip-type")
    # cleanear chip (to remove 'v1' from cetain gene/exon arrays):
    clean.chip <- tolower(gsub("v1", "",clean.chip))
    
  } else if (!is.null(celFile1) & !is.null(celFile2)) {
    celHead1 <- readCelHeader(celFile1)
    celHead2 <- readCelHeader(celFile2)
    if (identical(celHead1$chiptype, celHead2$chiptype)) {
      assign("chip", celHead1$chiptype)
      clean.chip <- tolower(gsub("-|_", "",chip))
      # cleanear chip (to remove 'v1' from cetain gene/exon arrays):
      clean.chip <- tolower(gsub("v1", "",clean.chip))
    }
  }

# CHECK IF CHIP-TYPE IS SUPPORTED -----------------------------------------

  if (chip %in% chipType){
    message(paste("GCSscore supports the selected chip-type: ",chip,sep=""))
  } else stop(paste("***The selected chip-type is not currently supported. \n       ***Please contact the maintaner for adding support for chip-type: ",chip,sep=""))
  # consider adding a list of the supported chip-types to the printed message
  
  pdpkg <- paste("pd.",tolower(gsub("-|_", ".",chip)),sep = "")
  chip.pd <- pdpkg

# 1. Start by checking if the pd.pkg for the chip-type is installed:
if (!requireNamespace(pdpkg, quietly = TRUE)){
  message(paste("Bioconductor platform design (.pd) package needs to be installed for chip-type: ",chip,sep=""))
  message(paste("installing .pd package: ",pdpkg,sep=""))
  # next command is throwing an error for  update (is.logical(update) is not TRUE), 
  # Bio requires update to be FALSE and not "never"
  BiocManager::install(pdpkg,ask = FALSE,update = FALSE)
} else {message(paste("Bioconductor platform design (.pd) package already installed for: ", chip,sep=""))}

# CHECK FOR INSTALLATION OF GCSscore-BUILT PACKAGE: xxx.probeFile --------
  # check if probeFile pacakge is installed:
  # 1. Pre-set flag for GCSscore built package (.probeFile) to be installed
  install.pF <- 1
  probepkg <- paste(clean.chip,".probeFile",sep="")
  message(paste("GCS-score analysis initiated for chip-type: ", chip,sep=""))
  
  # 2. Check if probeFile is installed and up-to-date:
if (!requireNamespace(probepkg, quietly = TRUE)){
  message(paste("GCS-score 'probeFile' package needs to be installed for chip-type: ", chip,sep=""))
  # set install flag to 1
  install.pF <- 1
} else if (loadNamespace(probepkg)[[".__NAMESPACE__."]][["spec"]][["version"]] < probepkg.vers) {
  message(paste("GCS-score 'probeFile' package for chip-type: (", chip,") needs to be updated to: ",probepkg.vers, sep=""))
  # set install flag to 1
  install.pF <- 1
} else {message(paste("The latest verison (", probepkg.vers,") of the GCS-score 'probeFile' package already installed for chip-type: ", chip,sep=""))
  # set install flag to 0
  install.pF <- 0
}

  # THE TYPEFILTER STUFF:
  if (typeFilter==0){
    # probeFile <- probeFile
    message("typeFilter set to (0) by default for best .CEL file scaling and normalization statistics \n all probe_id types (including control probe_ids and bgp probe_ids) will be for GCS-score calculation.")
    message("it is generally recommended to leave the typeFilter option set to (0)")
  } else{
    # probeFile <- probeFile[probesetid %in% annot$probe_id]
    message("typeFilter option (1) has been disabled in this release of the GCSscore package.")
    message("Pre-filtering for well annotated genes is deterimetnal to the DEG results. \n")
    message("typeFilter is automatically set to (0) for best .CEL file scaling and normalization statistics \n all probe_id types (including control probe_ids and bgp probe_ids) will be for GCS-score calculation.")
  }
  
# SETUP / PACAKGE CREATION FOR XTA-ARRAYS ---------------------------------
  if (chip %in% chipXTA){
    pF.type <- "XTA"
    message(paste("** Checking if Bioconductor annotation (.db) packages are installed chip-type: ", chip," **  \n",sep=""))
    # check if transcriptcluser.db annotation package is installed:
    packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
    annotName <- paste(clean.chip, "transcriptcluster", sep = "")
    if (!requireNamespace(packageName, quietly = TRUE)){
      message(paste("transcriptcluster annotation (.db) package not installed for chip-type: ", chip,sep=""))
      message(paste("installing (.db) package: ",packageName,sep=""))
      BiocManager::install(packageName,ask = FALSE, update = FALSE)
    }  else{message(paste("Annotation package (",packageName,") already installed for chip-type: ", chip,sep=""))}
    # check if probeset.db annotation package is installed:
    packageName <- paste(clean.chip, "probeset.db", sep = "")
    annotName <- paste(clean.chip, "probeset", sep = "")
    if (!requireNamespace(packageName, quietly = TRUE)){
      message(paste("probeset annotation (.db) package not installed for chip-type: ", chip,sep=""))
      message(paste("installing (.db) package: ",packageName,sep=""))
      BiocManager::install(packageName,ask = FALSE, update = FALSE)
    }  else{message(paste("Annotation package (",packageName,") already installed for chip-type: ", chip,sep=""))}
    
    # While performing checks, get species information for packageName:
    species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
    # replace ' ' with '_' to match makeProbePackage():
    species.pd <- gsub(" ", "_",species.pd)
    
    if (install.pF==1){
      # run probeFile builder functions for the given 'species.pd':
      message(paste("\n ** Generating probeFile package for chip-type: ",chip," **",sep = ""))
      message(" *** This may take a few minutes for XTA-style arrays! ***")
      message("*This package will only be generated once for each chip-type*\n*Or if probeFile version for the chip-type needs to be updated*")
      ClariomDXTApFBuilder(chip.pd = chip.pd, clean.chip = clean.chip, species.pd = species.pd, pF.type = pF.type)
    }
    
    annotTCid <- paste(clean.chip, ".TC.netaffx.annot", sep = "")
    if (!requireNamespace(annotTCid, quietly = TRUE)){
      message(paste("transcript-level netaffx-based annotation package needs to be installed for chip-type: ", chip,sep=""))
      message(paste("\n ** Generating transcript-level annotation package for chip-type: ",chip," **",sep = ""))
      packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
      annotName <- paste(clean.chip, "transcriptcluster", sep = "")
      netaffxAnnotBuilderXTA(chip.pd = chip.pd, clean.chip = clean.chip, species.pd = species.pd,packageName = packageName,annotName = annotName)
      
    } else if (loadNamespace(annotTCid)[[".__NAMESPACE__."]][["spec"]][["version"]] < annot.vers) {
      message(paste("transcript-level netaffx-based annotation package (", chip,") needs to be updated to: ",annot.vers, sep=""))
      message(paste("\n ** Generating transcript-level annotation package for chip-type: ",chip," **",sep = ""))
      packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
      annotName <- paste(clean.chip, "transcriptcluster", sep = "")
      netaffxAnnotBuilderXTA(chip.pd = chip.pd ,clean.chip = clean.chip,species.pd = species.pd,packageName = packageName,annotName = annotName)
      
    } else {message(paste("The latest verison (", annot.vers,") of the transcript-level netaffx-based annotation package already installed for chip-type: ", chip,sep=""))}
    
    # 2. check probesetid-level annotation package:
    annotPSR <- paste(clean.chip, ".PSR.netaffx.annot", sep = "")
    if (!requireNamespace(annotPSR, quietly = TRUE)){
      message(paste("transcript-level netaffx-based annotation package needs to be installed for chip-type: ", chip,sep=""))
      message(paste("\n ** Generating probeset-level annotation package for chip-type: ",chip," **",sep = ""))
      packageName <- paste(clean.chip, "probeset.db", sep = "")
      annotName <- paste(clean.chip, "probeset", sep = "")
      netaffxAnnotBuilderXTA(chip.pd = chip.pd ,clean.chip,species.pd,packageName,annotName)
      
    } else if (loadNamespace(annotPSR)[[".__NAMESPACE__."]][["spec"]][["version"]] < annot.vers) {
      message(paste("transcript-level netaffx-based annotation package (", chip,") needs to be updated to: ",annot.vers, sep=""))
      message(paste("\n ** Generating probeset-level annotation package for chip-type: ",chip," **",sep = ""))
      packageName <- paste(clean.chip, "probeset.db", sep = "")
      annotName <- paste(clean.chip, "probeset", sep = "")
      netaffxAnnotBuilderXTA(chip.pd = chip.pd,clean.chip,species.pd,packageName,annotName)
      
    } else {message(paste("The latest verison (", annot.vers,") of the probeset-level netaffx-based annotation package already installed for chip-type: ", chip,sep=""))}
    
    
    
    # SETUP METHODS-BASED FOR XTA-ARRAYS ----------------------------------------
    
    # Values shared by TCid (method = 1) and PSRid (method = 2) analysis:
    trim <- 0.04
    # LOAD PROBEFILE FOR CURRENT CHIP:
    probeFile <- eval(parse(text = paste(clean.chip,".probeFile::",clean.chip,".probeFile",sep="")))
    message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
    #Parse bgp probe list (same for all XTA probeFiles):  
    bgp <- probeFile[probesetid %like% "AFFX-BkGr-GC"]
    
    # REMOVED THIS >4 probes FILTER:------------------
    # # 03.08.20: Filter the probeFile for "probesetids" with >=4 probes.
    # info <- probeFile[,.(nProbes = .N), keyby = "probesetid"]
    # setkey(info,nProbes)
    # info <- info[nProbes >= 4]
    # probeFile <- probeFile[probesetid %in% info$probesetid]
    # # END OF 03.08.20 edits
    # REMOVED THIS >4 probes FILTER:------------------
    
    if (method == 1){
      message(" ** performing transcriptclusterid-level (gene-centric) level analysis ** \n")
      methodTag <- "gene_level"
      method <- "transcriptclusterid"
      # EDIT 03.08.20:
      # WE NEED TO GO AHEAD AND REMOVE THE JUCS FROM THE INFO CALCULATION FOR TC-ID METHOD:
      info <- probeFile[,.(nProbes = .N), keyby = method]
      # setkey(info,nProbes)
      infoKey <- key(info)
      # Now, create a TCid-only probeFile to get the 'info' for later without the JUCs:
      probeFile.TCid <- probeFile[probesetid %!like% "JUC"]
      # Get the 'info' for just the TCid-only probes (no JUC probes)
      info <- probeFile.TCid[,.(nProbes = .N), keyby = method]
      # setkey(info,nProbes)
      # END OF EDITS ON 03.08.20:
      
      # Load the correct annotation package for method:
      netaffx.annot <- eval(parse(text = paste(clean.chip,".TC.netaffx.annot::",clean.chip,".TC.netaffx.annot",sep="")))
      message(paste("loading annotations from package: ", clean.chip,".TC.netaffx.annot",sep = ""))
      
      # Ensure keys are correctly set:
      if (!identical(key(info), method)) setkeyv(info, method)
      
      if (!identical(key(netaffx.annot), key(info))) {
        setkeyv(netaffx.annot, key(info))
      }
      infoKey <- key(info)
      
      # Merge the info and the netaffx.annot together:
      info <- netaffx.annot[info]
      # Reset the key of info back to the 'method' OR back to its assigned key:
      setkeyv(info, infoKey)
      # determine last column of 'info', so Sscores can be appended to 'info'
      cIdx <- length(colnames(info))
    } else if (method == 2){
      message(" ** performing probeset-level (exon-centric) level analysis ** \n")
      methodTag <- "exon_level"
      method <- "probesetid"
      info <- probeFile[,.(nProbes = .N), keyby = method]
      # setkey(info,nProbes)
      infoKey <- key(info)
      
      # Load the correct annotation package for method:
      netaffx.annot <- eval(parse(text = paste(clean.chip,".PSR.netaffx.annot::",clean.chip,".PSR.netaffx.annot",sep="")))
      message(paste("loading annotations from package: ", clean.chip,".PSR.netaffx.annot",sep = ""))
      
      # Ensure keys are correctly set:
      if (!identical(key(info), method)) setkeyv(info, method)
      
      if (!identical(key(netaffx.annot), key(info))) {
        setkeyv(netaffx.annot, key(info))
      }
      infoKey <- key(info)
      
      # Merge the info and the netaffx.annot together:
      info <- netaffx.annot[info]
      # Reset the key of info back to the 'method' OR back to its assigned key:
      setkeyv(info, infoKey)
      # determine last column of 'info', so Sscores can be appended to 'info'
      cIdx <- length(colnames(info))
    } else stop("Select valid 'method' for XTA-type arrays: 1 (gene-Level) and 2(exon-level)")
    
    # BOTH METHODS: Recheck keys, and set key of bgp and probeFile to the method:
    
    if (!identical(key(info), method)) setkeyv(info, method)
    
    if (!identical(key(probeFile), key(info))) {
      setkeyv(probeFile, key(info))
      setkeyv(bgp, key(info))
    }
    infoKey <- key(info)
  }
  
# SETUP / PACAKGE CREATION FOR ClariomS-ARRAYS -------------------------
    # See script: "package_setup_ClariomS.R"
# GCSscore package setup for clariomS arrays (03.08.20):

if (chip %in% chipClarS){
  pF.type <- "clariomS"
  message(paste("** Checking if Bioconductor annotation (.db) packages are installed chip-type: ", chip," **  \n",sep=""))
  # check if transcriptcluser.db annotation package is installed:
  packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
  annotName <- paste(clean.chip, "transcriptcluster", sep = "")
  if (!requireNamespace(packageName, quietly = TRUE)){
    message(paste("transcriptcluster annotation (.db) package not installed for chip-type: ", chip,sep=""))
    message(paste("installing (.db) package: ",packageName,sep=""))
    BiocManager::install(packageName,ask = FALSE, update = FALSE)
  }  else{message(paste("Annotation package (",packageName,") already installed for chip-type: ", chip,sep=""))}
  
  # While performing checks, get species information for packageName:
  # species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
  # THERE IS NO NEED TO LOAD THE WHOLE ANNOTATION PACKAGE HERE
  species.pd <- eval(parse(text= paste(packageName,"::",annotName,"ORGANISM",sep = "")))
  # replace ' ' with '_' to match makeProbePackage():
  species.pd <- gsub(" ", "_",species.pd)
  
  # check install the updated probeFile package, if install.pF==1
  if (install.pF==1){
    # run probeFile builder functions for the given 'species.pd':
    message(paste("\n ** Generating probeFile package for chip-type: ",chip," **",sep = ""))
    # message(" *** This may take a few minutes for XTA-style arrays! ***")
    message("*This package will only be generated once for each chip-type*\n*Or if probeFile version for the chip-type needs to be updated*")
    ClariomSpFBuilder(chip.pd = pdpkg, clean.chip = clean.chip, species.pd = species.pd, pF.type = pF.type)
  }
  
  # Now check for each netaffx.annot package:
  # 1. check transcript-level annotation package:
  annotTCid <- paste(clean.chip, ".TC.netaffx.annot", sep = "")
  if (!requireNamespace(annotTCid, quietly = TRUE)){
    message(paste("transcript-level netaffx-based annotation package needs to be installed for chip-type: ", chip,sep=""))

    message(paste("\n ** Generating transcript-level annotation package for chip-type: ",chip," **",sep = ""))
    packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
    annotName <- paste(clean.chip, "transcriptcluster", sep = "")
    netaffxAnnotBuilderClarS(chip.pd = pdpkg ,clean.chip,species.pd,packageName,annotName)
    
  } else if (loadNamespace(annotTCid)[[".__NAMESPACE__."]][["spec"]][["version"]] < annot.vers) {
    message(paste("transcript-level netaffx-based annotation package (", chip,") needs to be updated to: ",annot.vers, sep=""))

    message(paste("\n ** Generating transcript-level annotation package for chip-type: ",chip," **",sep = ""))
    packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
    annotName <- paste(clean.chip, "transcriptcluster", sep = "")
    netaffxAnnotBuilderClarS(chip.pd = pdpkg ,clean.chip,species.pd,packageName,annotName)
    
  } else {message(paste("The latest verison (", annot.vers,") of the transcript-level netaffx-based annotation package already installed for chip-type: ", chip,sep=""))
  }
  # QUICK VIEW OF ANNOT FILE:
  # test <- clariomsmouse.TC.netaffx.annot::clariomsmouse.TC.netaffx.annot
  
  trim <- 0.04
  message(" *NOTE: 'method' is automatically set to (1) (transcript_cluster_id-level) for all ClariomS arrays")
  methodTag <- "gene_level"
  method <- "transcriptclusterid"
  message(" *Performing gene-level analysis (by 'transcript_cluster_id')\n")
  
  probeFile <- eval(parse(text = paste(clean.chip,".probeFile::",clean.chip,".probeFile",sep="")))
  message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
  #Parse bgp probe list (same for all XTA probeFiles):  
  bgp <- probeFile[transcriptclusterid %like% "AFFX-BkGr-GC"]
  info <- probeFile[,.(nProbes = .N), keyby = method]
  infoKey <- key(info)
  
  # Remove TCids with less than 4 probes:
  # setkey(info,nProbes)
  # info <- info[nProbes >= 4]
  # probeFile <- probeFile[transcriptclusterid %in% info$transcriptclusterid]
  # setkey(info,"transcriptclusterid")
  
  # Load the correct annotation package for method:
  netaffx.annot <- eval(parse(text = paste(clean.chip,".TC.netaffx.annot::",clean.chip,".TC.netaffx.annot::",sep="")))
  message(paste("loading annotations from package: ", clean.chip,".TC.netaffx.annot",sep = ""))
  
  # Ensure keys are correctly set:
  if (!identical(key(info), method)) setkeyv(info, method)
  
  if (!identical(key(netaffx.annot), key(info))) {
    setkeyv(netaffx.annot, key(info))
  }
  infoKey <- key(info)
  
  # Reset the key of info back to the 'method' OR back to its assigned key:
  setkeyv(info, infoKey)
  # Merge the info and the netaffx.annot together:
  info <- netaffx.annot[info]
  # determine last column of 'info', so Sscores can be appended to 'info'
  cIdx <- length(colnames(info))
  
  # Recheck keys, and set key of bgp and probeFile to the method:
  if (!identical(key(info), method)) setkeyv(info, method)
  
  if (!identical(key(probeFile), key(info))) {
    setkeyv(probeFile, key(info))
    setkeyv(bgp, key(info))
  }
  infoKey <- key(info)
}
#  INSERT SECTION FOR 3'IVT ARRAYS:

#  INSERT SECTION FOR 3' IVT ARRAYS---------------
# See script: "package_setup_3primeIVT.R"
# GCSscore package setup for 3' IVT arrays (03.08.20):  

if (chip %in% chip3IVT){
  trim <- 0.02
  pF.type <- "3primeIVT"
  # detailed 'method' for 3prime IVT arrays (PM-bkg_gc OR PM-MM)
  if (method == 1){
    method <- methodTag <- "gc"
    message(" *Performing gene-level analysis with GC-based bkg correction (GCS-score)\n")
  } else if (method == 2){
    method <- methodTag <- "pmmm"
    message(" *Performing gene-level analysis with PM-MM bkg correction (S-score legacy)\n")
  }
  
  # check to make sure proper annotation package is installed:
  packageName <- paste(clean.chip, ".db", sep = "")
  annotName <- paste(clean.chip, "", sep = "")
  if (!requireNamespace(packageName, quietly = TRUE)){
    message(paste("annotation (.db) package not installed for chip-type: ", chip,sep=""))
    message(paste("installing (.db) package: ",packageName,sep=""))
    BiocManager::install(packageName,ask = FALSE, update = FALSE)
  } else{message(paste("annotation (.db) package (",packageName,") already installed for chip-type: ", chip,sep=""))}
  
  # While performing checks, get species information for packageName:
  species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
  # replace ' ' with '_' to match makeProbePackage():
  species.pd <- gsub(" ", "_",species.pd)
  
  # check install the updated probeFile package, if install.pF==1
  if (install.pF==1){
    # run probeFile builder functions for the given 'species.pd':
    message(paste("\n ** Generating probeFile package for chip-type: ",chip," **",sep = ""))
    # message(" *** This may take a few minutes for XTA-style arrays! ***")
    message("*This package will only be generated once for each chip-type*\n*Or if probeFile version for the chip-type needs to be updated*")
    IVT3primepFBuilder(chip.pd = pdpkg, clean.chip = clean.chip, species.pd = species.pd, pF.type = pF.type)
  }
  
  probeFile <- eval(parse(text = paste(clean.chip,".probeFile::",clean.chip,".probeFile",sep="")))
  message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
  probeFile[,MM_fid := (fid+celHead1$cols)]
  
  # For method = 'gc', create 'bgp' file to simplify code:
  bgp <- probeFile[,c("probesetid","x","y","GC.count","MM_fid")]
  # rename "MM_fid" to "fid" to avoid coding edits:
  bgp[,fid := MM_fid]
  bgp[,MM_fid := NULL]
  # set keys to fid to align the PM probes and the BGP probes (for PM-MM):
  setkey(bgp,fid)
  setkey(probeFile,fid)
  
  # get info (gene symbol and gene name) from annotation package:
  annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
  annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
  # set keys to 'probe_id' and merge the two annotations together:
  setkey(annot.symbol,probe_id)
  setkey(annot.genename,probe_id)
  annot <- annot.symbol[annot.genename]
  message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
  
    # Filter down top probeFile to probesets that are annotated in the .db package:
    # probeFile <- probeFile[probesetid %in% annot$probe_id]
    # # Filter down bgp [equivalent to probeFile for MM probes] to probesets that are annotated in the .db package:
    # bgp <- bgp[probesetid %in% annot$probe_id]
  
  # For all 3prime IVT arrays, infokey is automatically set to: probesetID (man_fsetid)
  info <- probeFile[,.N, keyby = .(probesetid)]
  infoKey <- key(info)
  
  # Add annotations to 'info' before running GCS-score algorithm:
  # Have it such that The full list of TCids in the GCSscore is output, and the ones with genenames/symbols are annotated and the others are '---' or NA
  setkey(annot,probe_id)
  info <- annot[info]
  # Rename first column from 'probe_id' back to the infoKey' object:
  names(info)[which(names(info)=="probe_id")] <- infoKey
  # Reset the key of info back to the 'method' OR back to its assigned key:
  # setkeyv(info, method)
  setkeyv(info, infoKey)
  # determine last column of 'info', so Sscores can be appended to 'info'
  cIdx <- length(colnames(info))
}
# Go ahead and create a seperate annotation 'info' for later expressionSet.
  info.annot <- info
  setkeyv(info.annot,infoKey)
  
  # Disable rm.outmask for all non 3'-IVT arrays:
  if (chip %!in% chip3IVT){
    rm.outmask <- FALSE
    message("removal of outliers and masked probes is disabling in non 3'-IVT chip types")
  }

# mark the start of the GCS-score calculations:
message("")
message("*********************************")
message("**** begin GCS-score analysis ***")
message("*********************************")
message("")
# setup single run or batch run:
if (is.character(method)) {
  if (!is.null(celTable))  {
    for (i in 1:nrow(celTable)) {
      message("reading .CEL files")
      cel1 <- readCel(celTable[[2]][i], readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE,readOutliers = TRUE, readMasked = TRUE)
      cel2 <- readCel(celTable[[3]][i], readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE,readOutliers = TRUE, readMasked = TRUE)
      if (rm.outmask== TRUE){
        # Create list of all 'fid' probe indices which have been flagged in either celFile1, celFile2, or both:
        probes.rm <- rbind(as.data.table(cel1[["masked"]]),
                           as.data.table(cel1[["outliers"]]),
                           as.data.table(cel2[["outliers"]]),
                           as.data.table(cel2[["masked"]]))
        names(probes.rm) <- "probes.rm"
        # Remove the duplicated values to create clean 'probes.rm' object:
        probes.rm <- probes.rm[!duplicated(probes.rm)]
        # remove outlier and masked probes from the probeFile (call it probeFile.temp, so it can be overwritten inside the loop!)
        # OKAY, for the mouse 430_2, leave the 'MM_fid' in the probeFile and assign all of the "pmmm" method math to be within the probeFile!
        # I can keep the filtered 'bgp' list for all of the other array types (and the GC_BKG method for the 4302)
        # if ("MM_fid" %in% colnames(probeFile)){}
        if (method == "pmmm"){
          # remove all probe fids that are either MM or PM probes from 3' IVT arrays:
          probeFile.temp <- probeFile[!(fid %in% probes.rm$probes.rm | MM_fid %in% probes.rm$probes.rm)]
          # set dummy variable for 'bgp', since it is unused for the PM-MM method:
          bgp.temp <- bgp[!(fid %in% probes.rm$probes.rm)]
        }
        # the 'else' condition:   
        if (method != "pmmm"){ 
          
          probeFile.temp <- probeFile[!(fid %in% probes.rm$probes.rm)]
          # also remove outlier probes from the 'bgp' file:
          bgp.temp <- bgp[!(fid %in% probes.rm$probes.rm)]
        }
      }
      if (rm.outmask == FALSE){
        # Don't filter out any of the outlier/masked probes before analysis:
        probeFile.temp <- probeFile
        bgp.temp <- bgp
      }
      # Score <- computeSscore(cel1, cel2, probeFile, bgp, method, infoKey, SF1 = SF1, SF2 = SF2, verbose = verbose, trim = trim)
      Score <- computeSscore(cel1, cel2, probeFile = probeFile.temp, bgp = bgp.temp, method, infoKey, SF1 = SF1, SF2 = SF2, verbose = verbose, trim = trim, clean.chip = clean.chip)
      info <- info[Score, on = infoKey]
      if (celTab.names){
        message("")
        message("*custom run name (column name) for BATCH job:")
        message("")
        colnames(info)[cIdx + i] <- as.character(celTable[i,1])
        message(paste("    ",as.character(celTable[i,1]),"\n",sep=""))
      }
      else{
        message("")
        message("*standard run name (column name) for BATCH job:")
        # Add filenames to colname (CEL_1 vs CEL_2)
        celName1 <- strsplit(cel1$header$filename, "/")[[1]]; celName1 <- celName1[length(celName1)]
        celName2 <- strsplit(cel2$header$filename, "/")[[1]]; celName2 <- celName2[length(celName2)]
        colnames(info)[cIdx + i] <- paste(celName1, "vs", celName2)
        message("")
        message(paste("    ",celName1,"vs", celName2, "\n"))
        message("")
      }
      message(" **** run completed ****")
      message(" ***********************")
      message("")
    }
  } else if (!is.null(celFile1) & !is.null(celFile2)) {
    message("reading .CEL files")
    cel1 <- readCel(celFile1, readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE, readOutliers = TRUE, readMasked = TRUE)
    cel2 <- readCel(celFile2, readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE, readOutliers = TRUE, readMasked = TRUE)
    if (rm.outmask==TRUE){
      # Create list of all 'fid' probe indices which have been flagged in either celFile1, celFile2, or both:
      probes.rm <- rbind(as.data.table(cel1[["masked"]]),
                         as.data.table(cel1[["outliers"]]),
                         as.data.table(cel2[["outliers"]]),
                         as.data.table(cel2[["masked"]]))
      names(probes.rm) <- "probes.rm"
      # Remove the duplicated values to create clean 'probes.rm' object:
      probes.rm <- probes.rm[!duplicated(probes.rm)]
      # remove outlier and masked probes from the probeFile (call it probeFile.temp, so it can be overwritten inside the loop!)
      # OKAY, for the mouse 430_2, leave the 'MM_fid' in the probeFile and assign all of the "pmmm" method math to be within the probeFile!
      # I can keep the filtered 'bgp' list for all of the other array types (and the GC_BKG method for the 4302)
      # if ("MM_fid" %in% colnames(probeFile)){}
      if (method == "pmmm"){
        # remove all probe fids that are either MM or PM probes from 3' IVT arrays:
        probeFile.temp <- probeFile[!(fid %in% probes.rm$probes.rm | MM_fid %in% probes.rm$probes.rm)]
        # set dummy variable for 'bgp', since it is unused for the PM-MM method:
        bgp.temp <- bgp[!(fid %in% probes.rm$probes.rm)]
      }
      # the 'else' condition:   
      if (method != "pmmm"){ 
        
        probeFile.temp <- probeFile[!(fid %in% probes.rm$probes.rm)]
        # also remove outlier probes from the 'bgp' file:
        bgp.temp <- bgp[!(fid %in% probes.rm$probes.rm)]
      }
    }
    if (rm.outmask == FALSE){
      # Don't filter out any of the outlier/masked probes before analysis:
      probeFile.temp <- probeFile
      bgp.temp <- bgp
    }
    Score <- computeSscore(cel1, cel2, probeFile = probeFile.temp, bgp = bgp.temp, method, infoKey, SF1 = SF1, SF2 = SF2, verbose = verbose, trim = trim, clean.chip = clean.chip)
    info <- info[Score, on = infoKey]
    colnames(info)[cIdx + 1] <- "Sscore"
    # Add message pad for one line space:
    message("")
    # Verbose output for single-runs:
    #   if (verbose){
    #   View (info)
    ## Check if Sscore depend on nProbes, also see the skew between files:
    #   plot(info$nProbes,info$Sscore)
    #  hist(info$Sscore,xlim = c(-6,6),breaks = seq(-50,50,.25))
    # }
  }
  
  # For array types with i.nProbes (ClariomS/XTA)"
  # in info.annot: replace nProbes with tallied i.nProbes (for AFFX ids):
  if (!is.na(match("i.nProbes", names(info.annot)))){
    info.annot[,nProbes := i.nProbes]
    info.annot[,i.nProbes := NULL]
  }
  # REPLACE ALL BLANK VALUES WITH "NA", IN BOTH INFO AND INFO.ANNOT:
  info[info==""] <- NA
  info.annot[info.annot==""] <- NA
  # ADDITIONAL: MTA 1.0 (TCID-LEVEL) HAS MANY SYMBOLS == "---" (FOR ~33,000 TCIDS)
  info.annot[info.annot=="---"] <- NA
  # ADD IN CODE TO OUTPUT GCS-score to Biobase data structure: ExpressionSet
  # Create matrix with all the GCS-score values
  # NOTE: first 4 columns contain 'featureSet' data, the remaining columns contain the differential expression values
  # GCSs.exprs <- as.matrix(info[,5:ncol(info)], header=TRUE,as.is=TRUE)
  GCSs.exprs <- as.matrix(info[,(cIdx+1):ncol(info)], header=TRUE,as.is=TRUE)
  # GCSs.fsetData <- new("AnnotatedDataFrame", data=info[,1:4])
  GCSs.fsetData <- new("AnnotatedDataFrame", data=as.data.frame(info.annot))
  GCSs.out <- ExpressionSet(assayData = GCSs.exprs,
                            featureData  = GCSs.fsetData,
                            annotation = packageName)
  # Make the 'featureNames' equal to the 'probe_id' (column #1 of info) for the relevant structures within the ExpressionSet:
  featureNames(GCSs.out) <- GCSs.out@featureData@data[[1]]
  # head(featureNames(GCSs.out))
  
  # Write GCS-score output to .CSV file using 'data.table' function: fwrite
  if (fileout) {
    #ID tag: a custom formatted timestamp applied to all written files
    ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
    fileName <- paste(chip, "_", methodTag, "_", "Sscore2_", ID, ".csv", sep = "")
    message("*writing the GCS-score output to .CSV file in the working directory")
    fwrite(info, fileName)
  }
  if (gzip) {
    message("*compressing the GSC-score output into .gz file")
    system(paste("gzip", fileName))
  }
}
message("** GCS-score analysis complete **")
message("*********************************")

# options(warn = 0)
# return(info)
# return data in a Bioconductor-friendly data structure.
return(GCSs.out)
}
