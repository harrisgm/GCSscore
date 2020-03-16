# Function to generate 'netaffx' annotations from platform 
# design (pd) packages from BioConductor.  The generated probeFile is stored 
# in a data.table within an .rda file. These annotations are then converted 
# to a custom 'probePackage' using the makeProbePacage() function from
# the AnnotationForge package.

# This function is for use on ClariomS Affymetrix arrays.

# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/

# # START OF TESTING VARIABLES
# chip <- "RTA-1_0"
# chip <- "MTA-1_0"
# chip <- "Mouse430_2"
# chip <- "Clariom_S_Mouse"
# clean.chip <- tolower(gsub("-|_", "",chip))
# clean.chip <- tolower(gsub("v1", "",clean.chip))
# probepkg <- paste(clean.chip,".probeFile",sep="")
# pdpkg <- paste("pd.",tolower(gsub("-|_", ".",chip)),sep = "")
# chip.pd <- pdpkg
# packageName <- paste(clean.chip, "transcriptcluster.db", sep = "")
# annotName <- paste(clean.chip, "transcriptcluster", sep = "")
# # packageName <- paste(clean.chip, "probeset.db", sep = "")
# # annotName <- paste(clean.chip, "probeset", sep = "")
# species.pd <- eval(parse(text= paste(packageName,"::",annotName,"ORGANISM",sep = "")))
# species.pd <- gsub(" ", "_",species.pd)
# # END OF TESTING VARIABLES, START OF FUNCTION TEST:
# 
# # netaffxAnnotBuilderClarS(chip.pd, clean.chip,species.pd,packageName,annotName)
# netaffxAnnotBuilderXTA(chip.pd, clean.chip,species.pd,packageName,annotName)
# END OF FUNCTION TEST:

netaffxAnnotBuilderClarS <- function(chip.pd = NULL, clean.chip = NULL, species.pd = NULL,
                                     packageName = NULL, annotName = NULL) {
  # Install necessary pd.* package if not already installed:
  # 03.02.20: is install.pd should not be needed in the function as it is
  # already in the GCSscore.R function
  # if (!requireNamespace(chip.pd, quietly = TRUE)){
  #   BiocManager::install(chip.pd)
  #   message("installing necessary platform design (pd) package from Bioconductor")}
  # 
  
  # # get chip name (without periods/dots) from chip.pd:
  # clean.chip <- stringr::str_remove_all(strsplit(chip.pd,"pd.")[[1]][2],"[.]")
  # # 03.02.20: the above line seems redundant, since clean.chip is the same thing (I think)
  
  # try the following conditional:
  if (packageName %like% "transcriptcluster"){
    message("Building full transcriptclusterid-based annotation file from Bioconductor sources")
    load(system.file("extdata","netaffxTranscript.rda",package = chip.pd))
    netaffx.TCid.all <- data.table(netaffxTranscript@data)
    
    netaffx.TCid.annot <- netaffx.TCid.all[,c("transcriptclusterid")]
    netaffx.TCid.annot[,symbol := sapply(strsplit(netaffx.TCid.all$geneassignment, " // "), "[", 2)]
    netaffx.TCid.annot[,ref_id := sapply(strsplit(netaffx.TCid.all$geneassignment, " // "), "[", 1)]
    netaffx.TCid.annot[,name := sapply(strsplit(netaffx.TCid.all$geneassignment, " // "), "[", 3)]
    netaffx.TCid.annot[,chr := sapply(strsplit(netaffx.TCid.all$seqname, " // "), "[", 1)]
    netaffx.TCid.annot[,start := netaffx.TCid.all$start]
    netaffx.TCid.annot[,stop := netaffx.TCid.all$stop]
    netaffx.TCid.annot[,nProbes := netaffx.TCid.all$totalprobes]
    netaffx.TCid.annot[,locustype := sapply(strsplit(netaffx.TCid.all$locustype, "///"), "[", 1 )]
    netaffx.TCid.annot[,category := netaffx.TCid.all$category]
    
    # Now intergrate the well-annotated symbols/names for the respective .db package:
    db.symbol <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"SYMBOL)",sep = "")))
    db.genename <- eval(parse(text= paste("as.data.table(",packageName,"::",annotName,"GENENAME)",sep = "")))
    # set keys to 'probe_id' and merge the two annotations together:
    setkey(db.symbol,probe_id)
    setkey(db.genename,probe_id)
    db.TCid.annot <- db.symbol[db.genename]
    # rename columns to db.sybmol and db.name:
    names(db.TCid.annot) <- c("transcriptclusterid","db.symbol","db.name")
    # reset the keys:
    setkey(db.TCid.annot,transcriptclusterid)
    setkey(netaffx.TCid.annot,transcriptclusterid)
    netaffx.db.annot <- db.TCid.annot[netaffx.TCid.annot]
    # outdir <- "~/Desktop/"
    # merge them into a new data.table:
    
    outdir <- tempdir()
    # datafile <- paste(chip,".netaffx.TCid.annot.tab",sep="")
    netaffx.annot.tab <- paste(clean.chip,".netaffx.annot.probe_tab",sep="")
    netaffx.annot.loc <- paste(outdir,"/",netaffx.annot.tab,sep="")
    fwrite(netaffx.db.annot,file=netaffx.annot.loc,sep = "\t")
    annot.pkg.name = paste(clean.chip,".TC.netaffx.annot",sep="")
    
    # # Running stock function from AnnotationForge package version 1.26.0
    # AnnotationForge::makeProbePackage(
    makeProbePackageGCSs(
      # arraytype = annot.pkg.name,
      arraytype = clean.chip,
      outdir = outdir,
      species = species.pd,
      chip.pd = chip.pd,
      maintainer= "Guy Harris <harrisgm@vcu.edu>",
      version = "0.0.5",
      pkgname = annot.pkg.name,
      pkg.db = packageName,
      datafile = netaffx.annot.loc,
      importfun = "getTCnetaffxdatXTA",
      check = FALSE)
    
    # If time, try making custom makeProbePackage (that has updated descriptions in a local directory):
    # MORE IMPORTANT: look at attaching 'ProbePkg-template' from AnnotationForge into the GCSscore code:
    # If keep multiple templates: you can have correct package descriptions when they are created.
    # SEE where 'createRes' section of makePrboePackage():
    # originDir = system.file("ProbePkg-template", package=thispkg)
    # makeGCSannotPackage(
    #   arraytype = annot.pkg.name,
    #   outdir = outdir,
    #   species = species.pd,
    #   maintainer= "Guy Harris <harrisgm@vcu.edu>",
    #   version = "0.0.1",
    #   datafile = netaffx.annot.loc,
    #   importfun = "getTCnetaffxdat",
    #   check = FALSE)
    
    message("This transcriptcluster/TCid annotation package contains data parsed from the following source: ")
    message(paste("  (1) The nettaffx annotations from BioConductor platform design (pd) package: ",chip.pd,sep=""))
    message(paste("  (2) The concise list of well annotated genes from BioConductor annotation data (.db) package: ", packageName,sep=""))
    message(paste("transcriptcluster/TCid annotation package installing for chip: ",clean.chip,sep=""))
    # install the package, in the tempdir(), to the R library using 'devtools':
    pkg.loc <- paste(outdir,"/",annot.pkg.name,sep="")
    # Prevent the install from prompting the user to update packages when installing
    devtools::install(pkg = pkg.loc,upgrade = "never")
  }
  else stop("Did not detect a 'transcriptclusterid' column. Please recheck that the chip type is 'ClariomS' before continuing")
  
  # INSTALL PACKAGE (SHOULD BE SAME FOR BOTH METHODS)
  # # Brief checks of newly created annotation package (using clariomsmouse:
  # x <- clariomsmouse.annot.TCid.netaffx::clariomsmouse.annot.TCid.netaffx
  # class(x$transcriptclusterid)
  # class(x[,symbol])
  # class(x[,stop])
  return(message("DONE"))	
}
