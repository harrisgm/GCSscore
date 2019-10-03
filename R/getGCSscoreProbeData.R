# Functions to read the generated .probe_tab probe-level annotation files:

## A function that reads tab-delimited probe sequence
## (and other stuff) files from Affymetrix

# Code from "AnnotationForge" package: 
# https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html

# modified from original AnnotationForge function name: getProbeDataAffy()

# Once the functions have been sourced, see run examples:
      # arraytype <- "mta10"
      # datafile <- "mta10probeFile.probe_tab"
      # # For structure of Data:
      # sapply(probeFile,class)
      # 
      # makeProbePackage(
      #   arraytype = "mta10",
      #   species = "Mus_musculus",
      #   maintainer= "Guy Harris <harrisgm@vcu.edu>",
      #   version = "0.0.1",
      #   datafile = "mta10probeFile.probe_tab",
      #   importfun = "getXTAprobefileData",
      #   check = FALSE)



##----------------------------------------------------------------------
## Copyright R. Gentleman and W. Huber, 2003, all rights reserved
##----------------------------------------------------------------------
# 

get3primeIVTprobefileData <- function(arraytype, datafile,
                                     pkgname = NULL, comparewithcdf = FALSE)
{
  # loadNamespace("affy")
  
  if(missing(datafile)) {
    datafile <- paste(arraytype, "_probe_tab", sep="")
  } else {
    if (is(datafile, "character")) {
      datafile <- file(datafile, "r")
      on.exit(close(datafile))
    }
  }
  
  # arraytype = affy::cleancdfname(arraytype, addcdf=FALSE)
  # cdfname   = affy::cleancdfname(arraytype)
  # if (is.null(pkgname))
  pkgname = paste(arraytype, ".probeFile", sep="")
  int <- integer(0)
  chr <- ""
  # what = list(chr, int, int, int, chr, chr)
  # For ClariomS probeFiles:
  what = list(int, chr, int, int, chr, int, int, int)
  ## If datafile is a connection, scan is like readLines(), so we don't want to skip a line.
  if(inherits(datafile, "connection")){
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what)
  }else{
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what, skip=1)
  }
  # Comment out the "mess" part of the code since this is totally custom:
  
  # if((any(unlist(head) != c("Probe Set Name", "Probe X", "Probe Y",
  #                           "Probe Interrogation Position", "Probe Sequence", "Target Strandedness")))
  #    &&
  #    (any(unlist(head) != c("Probe Set ID", "probe x", "probe y",
  #                           "probe interrogation position", "probe sequence", "target strandedness")))
  # ) {
  #   mess = paste("The data file", datafile, "does not have the expected column names",
  #                "in its header line. Please make sure it is the right data file. If you are",
  #                "positive, you may need to write a customized data import function",
  #                "to replace 'getProbeDataAffy'. You may use 'getProbeDataAffy' as a template.",
  #                "Please see the help files for the functions 'getProbeDataAffy' and",
  #                "'makeProbePackage', and the makeProbePackage vignette in package AnnotationForge.\n")
  #   stop(mess)
  # }
  
  for (i in which(sapply(what, class)=="numeric")) {
    z = which(is.na(dat[[i]]))
    if(length(z)>0)
      stop(paste("Corrupted data file: found non-number in line ", z[1],
                 " of column ", head[i], ": ", dat[z[1], i]), sep="") 
  }
  
  ## data frame with the probe data:
  # pt = data.frame(sequence = I(dat[[5]]),           ## character
  #                 x        = dat[[2]],  ## integer
  #                 y        = dat[[3]],  ## integer
  #                 Probe.Set.Name               = I(dat[[1]]),          ## character 
  #                 Probe.Interrogation.Position = dat[[4]], ## integer
  #                 Target.Strandedness          = dat[[6]])             ## factor
  
  ## data frame with the custom probeFile data:
  pt = data.frame(fid = dat[[1]],
                  sequence = dat[[2]],
                  fsetid = dat[[3]],
                  strand = dat[[4]],
                  man_fsetid = dat[[5]],
                  x = dat[[6]],
                  y = dat[[7]],
                  GC.count = dat[[8]])
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  if(is.character(datafile))
    datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
                       ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}

getWTgeneexonprobefileData <- function(arraytype, datafile,
                                pkgname = NULL, comparewithcdf = FALSE)
{
  # loadNamespace("affy")
  
  if(missing(datafile)) {
    datafile <- paste(arraytype, "_probe_tab", sep="")
  } else {
    if (is(datafile, "character")) {
      datafile <- file(datafile, "r")
      on.exit(close(datafile))
    }
  }
  
  # arraytype = affy::cleancdfname(arraytype, addcdf=FALSE)
  # cdfname   = affy::cleancdfname(arraytype)
  # if (is.null(pkgname))
  pkgname = paste(arraytype, ".probeFile", sep="")
  int <- integer(0)
  chr <- ""
  # what = list(chr, int, int, int, chr, chr)
  # For WT gene/exon probeFiles:
  # if NO remove 'seqence' column:
  what = list(int, int, int, int, int, int, int, int, int, int, int, int, int, int)
  # if DO NOT remove 'sequence' column, then:
  # what = list(int, chr, int, int, int, int, int, int, int, int, int, int, int, int, int)
  ## If datafile is a connection, scan is like readLines(), so we don't want to skip a line.
  if(inherits(datafile, "connection")){
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what)
  }else{
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what, skip=1)
  }
  # Comment out the "mess" part of the code since this is totally custom:
  
  # if((any(unlist(head) != c("Probe Set Name", "Probe X", "Probe Y",
  #                           "Probe Interrogation Position", "Probe Sequence", "Target Strandedness")))
  #    &&
  #    (any(unlist(head) != c("Probe Set ID", "probe x", "probe y",
  #                           "probe interrogation position", "probe sequence", "target strandedness")))
  # ) {
  #   mess = paste("The data file", datafile, "does not have the expected column names",
  #                "in its header line. Please make sure it is the right data file. If you are",
  #                "positive, you may need to write a customized data import function",
  #                "to replace 'getProbeDataAffy'. You may use 'getProbeDataAffy' as a template.",
  #                "Please see the help files for the functions 'getProbeDataAffy' and",
  #                "'makeProbePackage', and the makeProbePackage vignette in package AnnotationForge.\n")
  #   stop(mess)
  # }
  
  for (i in which(sapply(what, class)=="numeric")) {
    z = which(is.na(dat[[i]]))
    if(length(z)>0)
      stop(paste("Corrupted data file: found non-number in line ", z[1],
                 " of column ", head[i], ": ", dat[z[1], i]), sep="") 
  }
  
  ## data frame with the probe data:
  # pt = data.frame(sequence = I(dat[[5]]),           ## character
  #                 x        = dat[[2]],  ## integer
  #                 y        = dat[[3]],  ## integer
  #                 Probe.Set.Name               = I(dat[[1]]),          ## character 
  #                 Probe.Interrogation.Position = dat[[4]], ## integer
  #                 Target.Strandedness          = dat[[6]])             ## factor
  
  ## data frame with the custom probeFile data (if REMOVE the 'sequence' column:
  pt = data.frame(fid = dat[[1]],
                  fsetid = dat[[2]],
                  strand = dat[[3]],
                  start = dat[[4]],
                  stop = dat[[5]],
                  transcriptclusterid = dat[[6]],
                  exon_id = dat[[7]],
                  crosshyb_type = dat[[8]],
                  level = dat[[9]],
                  chrom = dat[[10]],
                  type = dat[[11]],
                  x = dat[[12]],
                  y = dat[[13]],
                  GC.count = dat[[14]])
  
  ## data frame with the custom probeFile data (including the 'sequence' column):
  # pt = data.frame(fid = dat[[1]],
  #                 sequence = dat[[2]],
  #                 fsetid = dat[[3]],
  #                 strand = dat[[4]],
  #                 start = dat[[5]],
  #                 stop = dat[[6]],
  #                 transcriptclusterid = dat[[7]],
  #                 exon_id = dat[[8]],
  #                 crosshyb_type = dat[[9]],
  #                 level = dat[[10]],
  #                 chrom = dat[[11]],
  #                 type = dat[[12]],
  #                 x = dat[[13]],
  #                 y = dat[[14]],
  #                 GC.count = dat[[15]])
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  if(is.character(datafile))
    datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
                       ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}

getClariomSprobefileData <- function(arraytype, datafile,
                             pkgname = NULL, comparewithcdf = FALSE)
{
  # loadNamespace("affy")
  
  if(missing(datafile)) {
    datafile <- paste(arraytype, "_probe_tab", sep="")
  } else {
    if (is(datafile, "character")) {
      datafile <- file(datafile, "r")
      on.exit(close(datafile))
    }
  }
  
  # arraytype = affy::cleancdfname(arraytype, addcdf=FALSE)
  # cdfname   = affy::cleancdfname(arraytype)
  # if (is.null(pkgname))
  pkgname = paste(arraytype, ".probeFile", sep="")
  int <- integer(0)
  chr <- ""
  # what = list(chr, int, int, int, chr, chr)
  # For ClariomS probeFiles:
  what = list(chr, chr, chr, int, chr, int, int, int, int, int)
  ## If datafile is a connection, scan is like readLines(), so we don't want to skip a line.
  if(inherits(datafile, "connection")){
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what)
  }else{
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what, skip=1)
  }
# Comment out the "mess" part of the code since this is totally custom:
  
  # if((any(unlist(head) != c("Probe Set Name", "Probe X", "Probe Y",
  #                           "Probe Interrogation Position", "Probe Sequence", "Target Strandedness")))
  #    &&
  #    (any(unlist(head) != c("Probe Set ID", "probe x", "probe y",
  #                           "probe interrogation position", "probe sequence", "target strandedness")))
  # ) {
  #   mess = paste("The data file", datafile, "does not have the expected column names",
  #                "in its header line. Please make sure it is the right data file. If you are",
  #                "positive, you may need to write a customized data import function",
  #                "to replace 'getProbeDataAffy'. You may use 'getProbeDataAffy' as a template.",
  #                "Please see the help files for the functions 'getProbeDataAffy' and",
  #                "'makeProbePackage', and the makeProbePackage vignette in package AnnotationForge.\n")
  #   stop(mess)
  # }
  
  for (i in which(sapply(what, class)=="numeric")) {
    z = which(is.na(dat[[i]]))
    if(length(z)>0)
      stop(paste("Corrupted data file: found non-number in line ", z[1],
                 " of column ", head[i], ": ", dat[z[1], i]), sep="") 
  }
  
  ## data frame with the probe data:
  # pt = data.frame(sequence = I(dat[[5]]),           ## character
  #                 x        = dat[[2]],  ## integer
  #                 y        = dat[[3]],  ## integer
  #                 Probe.Set.Name               = I(dat[[1]]),          ## character 
  #                 Probe.Interrogation.Position = dat[[4]], ## integer
  #                 Target.Strandedness          = dat[[6]])             ## factor
  
  ## data frame with the custom probeFile data:
  pt = data.frame(transcriptclusterid = dat[[1]],
                  category = dat[[2]],
                  locustype = dat[[3]],
                  fid = dat[[4]],
                  seqname = dat[[5]],
                  fsetid = dat[[6]],
                  type = dat[[7]],
                  x = dat[[8]],
                  y = dat[[9]],
                  GC.count = dat[[10]])
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  if(is.character(datafile))
    datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
                       ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}

getXTAprobefileData <- function(arraytype, datafile,
                                     pkgname = NULL, comparewithcdf = FALSE)
{
  # loadNamespace("affy")
  
  if(missing(datafile)) {
    datafile <- paste(arraytype, "_probe_tab", sep="")
  } else {
    if (is(datafile, "character")) {
      datafile <- file(datafile, "r")
      on.exit(close(datafile))
    }
  }
  
  # arraytype = affy::cleancdfname(arraytype, addcdf=FALSE)
  # cdfname   = affy::cleancdfname(arraytype)
  # if (is.null(pkgname))
  pkgname = paste(arraytype, ".probeFile", sep="")
  int <- integer(0)
  chr <- ""
  # what = list(chr, int, int, int, chr, chr)
  # For ClariomS probeFiles:
  what = list(chr, chr, chr, chr, chr, int, int, chr, int, int, int)
  ## If datafile is a connection, scan is like readLines(), so we don't want to skip a line.
  if(inherits(datafile, "connection")){
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what)
  }else{
    head <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, nlines=1, what="character")
    dat  <- scan(datafile, sep="\t", quiet=TRUE, multi.line = FALSE, what=what, skip=1)
  }
  # Comment out the "mess" part of the code since this is totally custom:
  
  # if((any(unlist(head) != c("Probe Set Name", "Probe X", "Probe Y",
  #                           "Probe Interrogation Position", "Probe Sequence", "Target Strandedness")))
  #    &&
  #    (any(unlist(head) != c("Probe Set ID", "probe x", "probe y",
  #                           "probe interrogation position", "probe sequence", "target strandedness")))
  # ) {
  #   mess = paste("The data file", datafile, "does not have the expected column names",
  #                "in its header line. Please make sure it is the right data file. If you are",
  #                "positive, you may need to write a customized data import function",
  #                "to replace 'getProbeDataAffy'. You may use 'getProbeDataAffy' as a template.",
  #                "Please see the help files for the functions 'getProbeDataAffy' and",
  #                "'makeProbePackage', and the makeProbePackage vignette in package AnnotationForge.\n")
  #   stop(mess)
  # }
  
  for (i in which(sapply(what, class)=="numeric")) {
    z = which(is.na(dat[[i]]))
    if(length(z)>0)
      stop(paste("Corrupted data file: found non-number in line ", z[1],
                 " of column ", head[i], ": ", dat[z[1], i]), sep="") 
  }
  
  ## data frame with the probe data:
  # pt = data.frame(sequence = I(dat[[5]]),           ## character
  #                 x        = dat[[2]],  ## integer
  #                 y        = dat[[3]],  ## integer
  #                 Probe.Set.Name               = I(dat[[1]]),          ## character 
  #                 Probe.Interrogation.Position = dat[[4]], ## integer
  #                 Target.Strandedness          = dat[[6]])             ## factor
  
  ## data frame with the custom probeFile data:
  pt = data.frame(probesetid = dat[[1]],
                  locustype = dat[[2]],
                  probesettype = dat[[3]],
                  transcriptclusterid = dat[[4]],
                  category = dat[[5]],
                  fid = dat[[6]],
                  fsetid = dat[[7]],
                  meta_fsetid = dat[[8]],
                  x = dat[[9]],
                  y = dat[[10]],
                  GC.count = dat[[11]])
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  if(is.character(datafile))
    datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
                       ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}


# NOTE: if maximum compression of the probe package is desired,
# source this modification of the makeProbePackage() function.
# The function has been altered to set 'compress' equal to 'xz'

# source: makeProbePackage() [ changed 'compress' to compress = "xz"]:
# Code from "AnnotationForge" package: 
# https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html
##----------------------------------------------------------------------
## Copyright R. Gentleman and W. Huber, 2003, all rights reserved
##----------------------------------------------------------------------
# 
# makeProbePackage <- function(arraytype,
#                              importfun = "getProbeDataAffy",
#                              maintainer,
#                              version,
#                              species,
#                              pkgname = NULL,
#                              outdir  = ".",
#                              quiet = FALSE, check = TRUE, build = TRUE, unlink = TRUE, ...)
# {
#   ## Bureucracy: check arguments
#   if (missing(maintainer) || !is.character(maintainer))
#     stop(paste("'maintainer' is missing or invalid. Please specify the maintainer of the",
#                "package that you want to create in the form: Your name <you@domain>", sep="\n"))
#   if (missing(version) || !is.character(version))
#     stop(paste("'version' is missing or invalid. Please select a version number larger",
#                "than those used for any previous versions of this package.", sep="\n"))
#   if(!exists(importfun) || !is.function(get(importfun)))
#     stop("'importfun' must be a function.")
#   if (missing(species) || !is.character(species))
#     stop(paste("'species' is missing or invalid. Please specify the species that the",
#                "package will pertain to using the form: Genus_species (e.g., Homo_sapiens).", sep = "\n"))
# 
#   ## Call the import function
#   ## importRes is a list with three elements:
#   ## $pkgname : package name
#   ## $dataEnv : environment containing data objects
#   ## $symVal  : named list with symbol-value substitutions
#   if (!quiet) cat("Importing the data.\n")
#   importRes <- do.call(importfun, c(arraytype = arraytype, pkgname = pkgname, list(...)))
# 
#   pkgname <- importRes$pkgname
#   thispkg <- "AnnotationForge"
#   desc    <- packageDescription(thispkg)
# 
#   stopifnot(desc$Package ==thispkg)
#   thispkgVers <- desc$Version
# 
#   symbolValues <- c(importRes$symVal,
#                     list(
#                       VERSION            = version,
#                       CREATOR            = paste("package", thispkg, "version", thispkgVers),
#                       ANNOTATIONDBIVERSION = thispkgVers,
#                       MAINTAINER         = maintainer,
#                       SPECIES            = species))
# 
#   ## Create package
#   createRes <- createPackage(pkgname,
#                              destinationDir = outdir,
#                              originDir = system.file("ProbePkg-template", package=thispkg),
#                              symbolValues = symbolValues,
#                              unlink = unlink, quiet = quiet)
# 
#   ## Write the data objects
#   if (!quiet) cat("Writing the data.\n")
#   save(list  = ls(importRes$dataEnv),
#        file  = file.path(createRes$pkgdir, "data", paste(pkgname, ".rda", sep="")),
#        envir = importRes$dataEnv,
#        # compress = TRUE
#        compress = "xz")
# 
#   R_exe <- file.path(R.home(), "bin", "R")
#   ## R CMD check
#   cdir <- getwd()
#   setwd(outdir)
#   on.exit(setwd(cdir))
#   if (check) {
#     if (!quiet)
#       cat("Checking the package.\n")
#     ## Capture output to avoid spewing on screen, then read from log
#     checkOut <- system(paste(R_exe, "CMD check", pkgname), intern=TRUE)
#     logFile <- file.path(paste(pkgname, "Rcheck", sep="."), "00check.log")
#     if (!file.exists(logFile)) {
#       stop(paste("Expected but did not find the log-file", logFile, "after R CMD check"))
#     } else {
#       thelines <- readLines(logFile)
#       warns <- grep("WARNING", thelines, value=TRUE)
#       errs  <- grep("ERROR", thelines, value=TRUE)
#       if (length(warns)>0)
#         cat("*** WARNINGS ***\n", warns)
#       if (length(errs)>0)
#         stop(errs)
#     }
#     if (unlink)
#       unlink(paste(pkgname, ".Rcheck", sep=""), recursive = TRUE)
#   }
# 
#   ## R CMD build
#   if (build) {
#     if (!quiet)
#       cat("Building the package.\n")
#     buildOut <- system(paste(R_exe, "CMD build"), intern=TRUE)
#   }
#   setwd(cdir)
#   return(pkgname)
# }
