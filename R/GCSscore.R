GCSscore <- function(celFile1 = NULL, celFile2 = NULL, celTable = NULL,
                     celTab.names = FALSE, typeFilter = 1, method = 1, 
                     rm.outmask = FALSE, SF1 = NULL, SF2 = NULL, 
                     fileout = FALSE, gzip = FALSE, verbose = FALSE) {

  # Basic testing that all entered parameters are set to usable values:
  if (!(typeFilter == 0 | typeFilter == 1))stop("typeFilter must be (0) or (1). (1) is the recommended setting")
  if (!(method == 1 | method == 2))stop("method must be (1) or (2)")
  # check that either celFile1/2 or celTable is selected:
  if (is.null(celTable)) stopifnot(!is.null(celFile1), !is.null(celFile2))
  else if (is.null(celFile1) & is.null(celFile2)) stopifnot(!is.null(celTable))
  else stop("input either: 2 .CEL files (single-run) or 'celTable' (batch-run)")
  # DEFINE CHIP-TYPES THAT THE SSCORE2 CODE IS COMPATIABLE WITH:
  # All chipXTA and chipClarS have been checked (all chips as of 09.30.2019):
  chip3IVT <- c("Rhesus", "Mouse430_2")
  chipGene <- c("MoGene-2_1-st", "MoGene-1_0-st-v1", "DroGene-1_0-st")
  chipExon <- c("MoEx-1_0-st-v1")
  chipXTA <- c("MTA-1_0","HTA-2_0","RTA-1_0","Clariom_D_Human")
  chipClarS <- c("Clariom_S_Mouse","Clariom_S_Human","Clariom_S_Rat")
  # chipTrans <- c("MTA-1_0", "Clariom_S_Mouse")
  
  chipType <- c(chip3IVT, chipGene, chipXTA, chipClarS)
  
  if (!is.null(celTable)) {
    if (!is.data.table(celTable)) celTable <- as.data.table(celTable)
    chipType1 <- chipType2 <- vector("character", length = nrow(celTable))
    for (i in 1:nrow(celTable)) {
      stopifnot (!is.null(celTable[[2]][i]), !is.null(celTable[[3]][i]))  
      chipType1[i] <- readCelHeader(celTable[[2]][i])$chiptype
      chipType2[i] <- readCelHeader(celTable[[3]][i])$chiptype
    }
    
    if (all.equal(chipType1, chipType2)){
      chip <- chipType1[1] 
      clean.chip <- tolower(gsub("-|_", "",chip))}
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
  # check if probeFile pacakge is installed:
  probepkg <- paste(clean.chip,".probeFile",sep="")
  pdpkg <- paste("pd.",tolower(gsub("-|_", ".",chip)),sep = "")
  message(paste("GCS-score analysis initiated for chip-type: ", clean.chip,sep=""))
  
  # if probeFile package is not installed:
  if (!requireNamespace(probepkg, quietly = TRUE)){
    # devtools::install_github(paste("harrisgm/",probepkg,sep=""),upgrade = "never")
    message(paste("GCS-score 'probeFile' package needs to be installed for chip-type: ", clean.chip,sep=""))
    # check if platform design package for chip type has been installed:
    if (!requireNamespace(pdpkg, quietly = TRUE)){
      message("the platform design (.pd) for this chip-type will also need to be installed")
      message(paste("installing .pd package: ",pdpkg,sep=""))
      BiocManager::install(pdpkg,ask = FALSE)
    } else {message(paste("platform design (.pd) package already installed for: ", clean.chip,sep=""))}
  } else {message(paste("GCS-score 'probeFile' package already installed for chip-type: ", clean.chip,sep=""))}
  
  
  
  if (chip %in% chipXTA){
    trim <- 0.04
    # First, check if .db annotation package is installed:
    # Select 'method' for probe tally (gene-level(1) or exon-level(2)):
    if (method == 1){
      methodTag <- "gene_level"
      method <- "transcriptclusterid"
      message(" *performing gene-level analysis (by 'transcript_cluster_id')\n")
      
      # check to make sure proper annotation package is installed:
      packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
      annotName <- paste(clean.chip, "transcriptcluster", sep = "")
      if (!requireNamespace(packageName, quietly = TRUE)){
        BiocManager::install(packageName,ask = FALSE)
        message("the annotation (.db) for this chip-type will also need to be installed")
        message(paste("installing .db package: ",packageName,sep=""))}
      
      # Check if 'probeFile' package is built and installed,
      # Build the probe package if not already installed:
      if (!requireNamespace(probepkg, quietly = TRUE)){
        species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
        # replace ' ' with '_' to match makeProbePackage():
        species.pd <- gsub(" ", "_",species.pd)
        # run probeFile builder functions for the given 'species.pd':
        message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
        message(" *** This may take several minutes for XTA-style arrays! ***")
        message("*This package will only be generated once for each chip-type* \n")
        ClariomDXTApFBuilder(chip.pd = pdpkg, species.pd = species.pd)
      } else {message(paste("probeFile from package: '", clean.chip,".probeFile' is already installed",sep = ""))}
      # Read in 'probeFile' probe-level data file from annotationForge-generated package:
      probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
      message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
      #Set bgp probe list:  
      bgp <- probeFile[probesetid %like% "AFFX-BkGr-GC"]
      
      # get info (gene symbol and gene name) from annotation package:
      annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
      annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
      # set keys to 'probe_id' and merge the two annotations together:
      setkey(annot.symbol,probe_id)
      setkey(annot.genename,probe_id)
      annot <- annot.symbol[annot.genename]
      message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
      
      # Filter probeFile down to annotated TCids from .db file:
      # Account for one conditional if chip type = "Clariom_D_Human":
      # In general, for gene-level XTA arrays: key by 'meta_fsetid' gives best match to .db files:
      # for 'Clariom_D_Human' --> must key by 'transcriptclusterid' or it will not work.
      if (chip == "Clariom_D_Human"){
        method <- "transcriptclusterid"
        probeFile <- probeFile[transcriptclusterid %in% annot$probe_id]
      } else {
        method <- "meta_fsetid"
        probeFile <- probeFile[meta_fsetid %in% annot$probe_id]
      }
    }
    # (exon-level) probesetid:
    else if (method == 2){
      methodTag <- "exon_level"
      method <- "probesetid"
      message(" *Performing exon-level analysis (by 'probeset_id')\n")

      # check to make sure annotation package is installed:
      packageName <- paste(clean.chip, "probeset.db", sep = "")
      annotName <- paste(clean.chip, "probeset", sep = "")
      if (!requireNamespace(packageName, quietly = TRUE)){
        BiocManager::install(packageName,ask = FALSE)}
      
      # Check if 'probeFile' package is built and installed,
      # Build the probe package if not already installed:
      if (!requireNamespace(probepkg, quietly = TRUE)){
        species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
        # replace ' ' with '_' to match makeProbePackage():
        species.pd <- gsub(" ", "_",species.pd)
        # run probeFile builder functions for the given 'species.pd':
        message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
        message(" *** This may take a few minutes to complete ***")
        message("*This package will only be generated once for each chip-type* \n")
        ClariomDXTApFBuilder(chip.pd = pdpkg, species.pd = species.pd)
      }
      # Read in 'probeFile' probe-level data file from annotationForge-generated package:
      probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
      message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
      #Set bgp probe list:  
      bgp <- probeFile[probesetid %like% "AFFX-BkGr-GC"]
      
      # get info (gene symbol and gene name) from annotation package:
      annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
      annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
      # set keys to 'probe_id' and merge the two annotations together:
      setkey(annot.symbol,probe_id)
      setkey(annot.genename,probe_id)
      annot <- annot.symbol[annot.genename]
      message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
      # Filter down to PSRids that are annotated in the .db package:
      probeFile <- probeFile[probesetid %in% annot$probe_id]
    }
    # else stop("Select valid 'method' for XTA-type arrays: 1 (gene-Level) and 2(exon-level)")
    
    # Key the probeFile to the method:
    info <- probeFile[,.(nProbes = .N), keyby = method]
    infoKey <- key(info)
    # Ensure keys are correctly set:
    if (!identical(key(info), method)) setkeyv(info, method)
    
    if (!identical(key(probeFile), key(info))) {
      setkeyv(probeFile, key(info))
      setkeyv(bgp, key(info))
    }
    infoKey <- key(info)
  }
  
  if (chip %in% chipClarS){
    trim <- 0.04
    methodTag <- "gene_level"
    method <- "transcriptclusterid"
    message(" *Performing gene-level analysis (by 'transcript_cluster_id')\n")
    if (method == 2){
    message(" *NOTE: 'method' is automatically set to (1) (transcript_cluster_id-level) for all ClariomS arrays")
    }
    # check to make sure annotation package is installed:
    packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
    annotName <- paste(clean.chip, "transcriptcluster", sep = "")
    if (!requireNamespace(packageName, quietly = TRUE)){
      BiocManager::install(packageName,ask = FALSE)}
    
    # Check if 'probeFile' package is built and installed,
    # Build the probe package if not already installed:
    if (!requireNamespace(probepkg, quietly = TRUE)){
      species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
      # replace ' ' with '_' to match makeProbePackage():
      species.pd <- gsub(" ", "_",species.pd)
      # run probeFile builder functions for the given 'species.pd':
      message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
      message(" *** This may take a few minutes to complete ***")
      message("*This package will only be generated once for each chip-type* \n")
      ClariomSpFBuilder(chip.pd = pdpkg, species.pd = species.pd)
    }
    
    # Read in 'probeFile' probe-level data file from annotationForge-generated package:
    probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
    message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
    #Set bgp probe list:  
    bgp <- probeFile[transcriptclusterid %like% "AFFX-BkGr-GC"]
    
    # get info (gene symbol and gene name) from annotation package:
    annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
    annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
    # set keys to 'probe_id' and merge the two annotations together:
    setkey(annot.symbol,probe_id)
    setkey(annot.genename,probe_id)
    annot <- annot.symbol[annot.genename]
    message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
    
    # Select the filter option:
    if (typeFilter==0){
      probeFile <- probeFile
      message("typeFilter set to (0) all probe_id types (including control probe_ids and bgp probe_ids) on chip will be for GCS-score calculation")
      message("it is generally recommended to leave the typeFilter option set to (1)")
    }
    else{ 
      # Filter down to TCids that are annotated in the .db package:
      probeFile <- probeFile[transcriptclusterid %in% annot$probe_id]
    }
    # Key the probeFile to the method:
    info <- probeFile[,.(nProbes = .N), keyby = method]
    infoKey <- key(info)
    # Ensure keys are correctly set:
    if (!identical(key(info), method)) setkeyv(info, method)
    
    if (!identical(key(probeFile), key(info))) {
      setkeyv(probeFile, key(info))
      setkeyv(bgp, key(info))
    }
    infoKey <- key(info)
  }
  
  
  
  if (chip %in% chipGene | chip %in% chipExon){
    trim <- 0.04
    probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
    message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
    #Set bgp probe list (use type_dict to select bgp probes -->  [7: control->bgp->antigenomic]):
    bgp <- probeFile[type == 7]
    # Select probe tally method:
    if (method ==1){
      # Select probe tally method:
      methodTag <- "gene_level"
      method <- "transcriptclusterid"
      message(" *Performing gene-level analysis (by 'transcript_cluster_id')\n")
      
      # check to make sure proper annotation package is installed:
      packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
      annotName <- paste(clean.chip, "transcriptcluster", sep = "")
      if (!requireNamespace(packageName, quietly = TRUE)){
        BiocManager::install(packageName,ask = FALSE)}
      
      # Check if 'probeFile' package is built and installed,
      # Build the probe package if not already installed:
      if (!requireNamespace(probepkg, quietly = TRUE)){
        species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
        # replace ' ' with '_' to match makeProbePackage():
        species.pd <- gsub(" ", "_",species.pd)
        # run probeFile builder functions for the given 'species.pd':
        message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
        message(" *** This may take a few minutes to complete ***")
        message("*This package will only be generated once for each chip-type* \n")
        WTgeneexonpFBuilder(chip.pd = pdpkg, species.pd = species.pd)
      }
      
      # Read in 'probeFile' probe-level data file from annotationForge-generated package:
      probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
      message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
      #Set bgp probe list (use type_dict to select bgp probes -->  [7: control->bgp->antigenomic]):
      bgp <- probeFile[type == 7]
      
      # get info (gene symbol and gene name) from annotation package:
      annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
      annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
      # set keys to 'probe_id' and merge the two annotations together:
      setkey(annot.symbol,probe_id)
      setkey(annot.genename,probe_id)
      annot <- annot.symbol[annot.genename]
      # For Gene and Exon Arrays, the 'probe_id' needs to be converted to 'integer' to match the probeFile:
      annot[,probe_id := as.integer(probe_id)]
      message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
      
      # Filter probeFile down to annotated TCids from .db file:
      if (typeFilter==0){
        probeFile <- probeFile
        message("*typeFilter set to (0) all probe_id types (including control probe_ids and bgp probe_ids) on chip will be for GCS-score calculation")
        message("**it is generally recommended to leave the typeFilter option set to (1)\n")
      }
      else{
        # Filter down to TCids that are annotated in the .db package:
        probeFile <- probeFile[transcriptclusterid %in% annot$probe_id]
      }
    }
    # (exon-level) probesetid:
    else if (method == 2){
      methodTag <- "exon_level"
      method <- "fsetid"
      message(" *Performing exon-level analysis (by 'probeset_id')\n")
      # check to make sure annotation package is installed:
      packageName <- paste(clean.chip, "probeset.db", sep = "")
      annotName <- paste(clean.chip, "probeset", sep = "")
      if (!requireNamespace(packageName, quietly = TRUE)){
        BiocManager::install(packageName,ask = FALSE)}
      
      # Check if 'probeFile' package is built and installed,
      # Build the probe package if not already installed:
      if (!requireNamespace(probepkg, quietly = TRUE)){
        species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
        # replace ' ' with '_' to match makeProbePackage():
        species.pd <- gsub(" ", "_",species.pd)
        # run probeFile builder functions for the given 'species.pd':
        message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
        message(" *** This may take a few minutes to complete ***")
        message("*This package will only be generated once for each chip-type* \n")
        WTgeneexonpFBuilder(chip.pd = pdpkg, species.pd = species.pd)
      }
      
      # Read in 'probeFile' probe-level data file from annotationForge-generated package:
      probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
      message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
      #Set bgp probe list (use type_dict to select bgp probes -->  [7: control->bgp->antigenomic]):
      bgp <- probeFile[type == 7]
      
      # get info (gene symbol and gene name) from annotation package:
      annot.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
      annot.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
      # set keys to 'probe_id' and merge the two annotations together:
      setkey(annot.symbol,probe_id)
      setkey(annot.genename,probe_id)
      annot <- annot.symbol[annot.genename]
      # For Gene and Exon Arrays, the 'probe_id' needs to be converted to 'integer' to match the probeFile:
      annot[,probe_id := as.integer(probe_id)]
      message(paste("loading SYMBOL and GENENAME annotations from BioConductor package: ", packageName,sep = ""))
      # Filter probeFile down to annotated probesetids from .db file:
      if (typeFilter==0){
        probeFile <- probeFile
        message(" *typeFilter set to (0) all probe_id types (including control probe_ids and bgp probe_ids) on chip will be for GCS-score calculation")
        message(" *it is generally recommended to leave the typeFilter option set to (1)")
      }
      # Filter down to probesetids  that are annotated in the .db package:
      else{probeFile <- probeFile[fsetid %in% annot$probe_id]}
      
    }
    # Key the probeFile to the method:
    info <- probeFile[,.(nProbes = .N), keyby = method]
    infoKey <- key(info)
    # Ensure keys are correctly set:
    if (!identical(key(info), method)) setkeyv(info, method)
    
    if (!identical(key(probeFile), key(info))) {
      setkeyv(probeFile, key(info))
      setkeyv(bgp, key(info))
    }
    infoKey <- key(info)
  }
  
  
  # 3prime IVT sytle arrays:
  if (chip %in% chip3IVT){
    trim <- 0.02
    # check to make sure proper annotation package is installed:
    packageName <- paste(clean.chip, ".db", sep = "")
    annotName <- paste(clean.chip, "", sep = "")
    if (!requireNamespace(packageName, quietly = TRUE)){
      BiocManager::install(packageName,ask = FALSE)}
    
    # Check if 'probeFile' package is built and installed,
    # Build the probe package if not already installed:
    if (!requireNamespace(probepkg, quietly = TRUE)){
      species.pd <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"ORGANISM)",sep = "")))
      # replace ' ' with '_' to match makeProbePackage():
      species.pd <- gsub(" ", "_",species.pd)
      # run probeFile builder functions for the given 'species.pd':
      message(paste("\n ** Generating probeFile package for chip-type: ",clean.chip," **",sep = ""))
      message(" *** This may take a few minutes to complete ***")
      message("*This package will only be generated once for each chip-type* \n")
      IVT3primepFBuilder(chip.pd = pdpkg, species.pd = species.pd)
    }
    
    # Read in 'probeFile' probe-level data file from annotationForge-generated package:
    probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
    message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
    
    probeFile <- eval(parse(text = paste("as.data.table(",clean.chip,".probeFile::",clean.chip,".probeFile)",sep="")))
    message(paste("loading probeFile from package: ", clean.chip,".probeFile",sep = ""))
    # get nrow/ncol (celHead1$rows/celHead1$cols):
    # Create the 'MM_fid' probeFile, which references the MM probe found 1 row beneath the PM probe:
    probeFile[,MM_fid := (fid+celHead1$cols)]
    
    # Select 'method' for 3prime IVT arrays (PM-bkg_gc OR PM-MM)
    if (method == 1){
      method <- methodTag <- "gc"
      message(" *Performing gene-level analysis with GC-based bkg correction (GCS-score)\n")
    } else if (method == 2){
      method <- methodTag <- "pmmm"
      message(" *Performing gene-level analysis with PM-MM bkg correction (S-score legacy)\n")
    }
    
    # For method = 'gc', create 'bgp' file to simplify code:
    bgp <- probeFile[,c("man_fsetid","x","y","GC.count","MM_fid")]
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
    
    # Filter probeFile down to annotated TCids from .db file:
    if (typeFilter==0){
      probeFile <- probeFile
      message("typeFilter set to (0) all probe_id types (including control probe_ids and bgp probe_ids) on chip will be for GCS-score calculation")
      message("it is generally recommended to leave the typeFilter option set to (1)")
    }
    if (typeFilter==1){
      # Filter down top probeFile to probesets that are annotated in the .db package:
      probeFile <- probeFile[man_fsetid %in% annot$probe_id]
      # Filter down bgp [equivalent to probeFile for MM probes] to probesets that are annotated in the .db package:
      bgp <- bgp[man_fsetid %in% annot$probe_id]
    }
    # For all 3prime IVT arrays, infokey is automatically set to: probesetID (man_fsetid):
    info <- probeFile[,.N, keyby = .(man_fsetid)]
    infoKey <- key(info)
  }
  
  # Add annotations to 'info' before running GCS-score algorithm:
  # Have it such that The full list of TCids in the GCSscore is output, and the ones with genenames/symbols are annotated and the others are '---' or NA
  setkey(annot,probe_id)
  info <- annot[info]
  # Rename first column from 'probe_id' back to the 'method' using save 'infoKey' object:
  # names(info)[which(names(info)=="probe_id")] <- method
  names(info)[which(names(info)=="probe_id")] <- infoKey
  # Reset the key of info back to the 'method' OR back to its assigned key:
  # setkeyv(info, method)
  setkeyv(info, infoKey)
  # determine last column of 'info', so Sscores can be appended to 'info'
  cIdx <- length(colnames(info))
  
  # mark the start of the GCS-score calculations:
  message("")
  message("*********************************")
  message("**** begin GCS-score analysis ***")
  message("")
  # setup single run or batch run:
  if (is.character(method)) {
    if (!is.null(celTable))  {
      for (i in 1:nrow(celTable)) {
        message("reading .CEL files")
        cel1 <- readCel(celTable[[2]][i], readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE,readOutliers = FALSE, readMasked = FALSE)
        cel2 <- readCel(celTable[[3]][i], readIntensities = TRUE, readStdvs = TRUE, readPixels = TRUE, readXY = TRUE,readOutliers = FALSE, readMasked = FALSE)
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
        Score <- computeSscore(cel1, cel2, probeFile = probeFile.temp, bgp = bgp.temp, method, infoKey, SF1 = SF1, SF2 = SF2, verbose = verbose, trim = trim)
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
      Score <- computeSscore(cel1, cel2, probeFile = probeFile.temp, bgp = bgp.temp, method, infoKey, SF1 = SF1, SF2 = SF2, verbose = verbose, trim = trim)
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
    # ADD IN CODE TO OUTPUT GCS-score to Biobase data structure: ExpressionSet
    # Create matrix with all the GCS-score values
    # NOTE: first 4 columns contain 'featureSet' data, the remaining columns contain the differential expression values
    GCSs.exprs <- as.matrix(info[,5:ncol(info)], header=TRUE,as.is=TRUE)
    GCSs.fsetData <- new("AnnotatedDataFrame", data=info[,1:4])
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
