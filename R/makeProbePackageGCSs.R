# Modified makeProbePackage() function for use with GCSscore-built 
# annotation and probeFile pacakges 

# The function and the accompanying ProbePkg-template files 
# were modified from source code contained within the Bioconductor
# package: AnnotationForge (1.28.0) (see original copyright below)
##----------------------------------------------------------------------
## Copyright R. Gentleman and W. Huber, 2003, all rights reserved
##----------------------------------------------------------------------
makeProbePackageGCSs <- function(arraytype,
                                 importfun = "getProbeDataAffy",
                                 maintainer,
                                 version,
                                 species,
                                 pkgname,
                                 chip.pd,
                                 pF.type,
                                 outdir,
                                 quiet = FALSE, check = TRUE, build = TRUE, unlink = TRUE, ...)
{
  ## Bureucracy: check arguments
  if (missing(maintainer) || !is.character(maintainer))
    stop(paste("'maintainer' is missing or invalid. Please specify the maintainer of the",
               "package that you want to create in the form: Your name <you@domain>", sep="\n"))
  if (missing(version) || !is.character(version))
    stop(paste("'version' is missing or invalid. Please select a version number larger",
               "than those used for any previous versions of this package.", sep="\n"))
  if(!exists(importfun) || !is.function(get(importfun)))
    stop("'importfun' must be a function.")
  if (missing(species) || !is.character(species))
    stop(paste("'species' is missing or invalid. Please specify the species that the",
               "package will pertain to using the form: Genus_species (e.g., Homo_sapiens).", sep = "\n"))
  
  ## Call the import function
  ## importRes is a list with three elements:
  ## $pkgname : package name
  ## $dataEnv : environment containing data objects
  ## $symVal  : named list with symbol-value substitutions
  if (!quiet) cat("Importing the data.\n")
  # importRes <- do.call(importfun, c(arraytype = arraytype, pkgname = pkgname, chip.pd = chip.pd, list(...)))
  importRes <- do.call(importfun, c(arraytype = arraytype, pkgname = pkgname, chip.pd = chip.pd, list(...)))
  
  
  pkgname <- importRes$pkgname
  thispkg <- "AnnotationForge"
  desc    <- packageDescription(thispkg)
  
  stopifnot(desc$Package ==thispkg)
  thispkgVers <- desc$Version
  
  symbolValues <- c(importRes$symVal,
                    list(
                      VERSION            = version,
                      CREATOR            = paste("package", thispkg, "version", thispkgVers),
                      ANNOTATIONDBIVERSION = thispkgVers,                                      
                      MAINTAINER         = maintainer,
                      SPECIES            = species,
                      CHIP.PD            = chip.pd))
  
  if (pF.type == "clariomS"){
    createRes <- createPackage(pkgname,
                               destinationDir = outdir,
                               # originDir = system.file("ProbePkg-template", package=thispkg),
                               # originDir = "~/Dropbox (Personal)/GCSscore_SUBMISSION_10162019/BioC_source/GCSscore/inst/ProbePkg-template-ClariomS",
                               # originDir = system.file("ProbePkg-template-ClariomS", package = "GCSscore"),
                               originDir = "~/Desktop/GitHub_syncing_rev2_2020/GitHub_new_clone/GCSscore/inst/ProbePkg-template-ClariomS/",
                               symbolValues = symbolValues,
                               unlink = unlink, quiet = quiet)
  } else if (pF.type == "XTA") {
    createRes <- createPackage(pkgname,
                               destinationDir = outdir,
                               # originDir = system.file("ProbePkg-template", package=thispkg),
                               # originDir = "~/Dropbox (Personal)/GCSscore_SUBMISSION_10162019/BioC_source/GCSscore/inst/ProbePkg-template-XTA",
                               originDir = system.file("ProbePkg-template-XTA", package = "GCSscore"),
                               symbolValues = symbolValues,
                               unlink = unlink, quiet = quiet)
  } else if (pF.type == "3primeIVT") {
    createRes <- createPackage(pkgname,
                               destinationDir = outdir,
                               # originDir = system.file("ProbePkg-template", package=thispkg),
                               originDir = "~/Dropbox (Personal)/GCSscore_SUBMISSION_10162019/BioC_source/GCSscore/inst/ProbePkg-template-3IVT",
                               # originDir = system.file("ProbePkg-template-3IVT", package = "GCSscore"),
                               symbolValues = symbolValues,
                               unlink = unlink, quiet = quiet)
  } else stop("probeFile cannot be built.  This chip-type/generation is not supported")
  
  ## Write the data objects
  if (!quiet) cat("Writing the data.\n")
  save(list  = ls(importRes$dataEnv),
       file  = file.path(createRes$pkgdir, "data", paste(pkgname, ".rda", sep="")),
       envir = importRes$dataEnv,
       compress = TRUE)
  
  R_exe <- file.path(R.home(), "bin", "R")
  ## R CMD check
  cdir <- getwd()
  setwd(outdir)
  on.exit(setwd(cdir))
  if (check) {
    if (!quiet)
      cat("Checking the package.\n")
    ## Capture output to avoid spewing on screen, then read from log
    checkOut <- system(paste(R_exe, "CMD check", pkgname), intern=TRUE)
    logFile <- file.path(paste(pkgname, "Rcheck", sep="."), "00check.log")
    if (!file.exists(logFile)) {
      stop(paste("Expected but did not find the log-file", logFile, "after R CMD check"))
    } else {
      thelines <- readLines(logFile)
      warns <- grep("WARNING", thelines, value=TRUE)
      errs  <- grep("ERROR", thelines, value=TRUE)
      if (length(warns)>0)
        cat("*** WARNINGS ***\n", warns)
      if (length(errs)>0)
        stop(errs)
    }
    if (unlink)
      unlink(paste(pkgname, ".Rcheck", sep=""), recursive = TRUE)
  }
  
  ## R CMD build
  if (build) {
    if (!quiet)
      cat("Building the package.\n")
    buildOut <- system(paste(R_exe, "CMD build"), intern=TRUE)
  }
  setwd(cdir)
  return(pkgname)
}

