
getTCnetaffxdatXTA <- function(arraytype, datafile,
                               pkgname, pkg.db, chip.pd, comparewithcdf = FALSE){
  
  dat <- data.table::fread(datafile)
  pt = data.table::data.table(transcriptclusterid = dat$transcriptclusterid,
                              symbol = dat$symbol,
                              name = dat$name,
                              ref_id = dat$ref_id,
                              db.symbol = dat$db.symbol,
                              db.name = dat$db.name,
                              chr = dat$chr,
                              start = dat$start,
                              stop = dat$stop,
                              nProbes = dat$nProbes,
                              locustype = dat$locustype,
                              category = dat$category,
                              stringsAsFactors = TRUE)
  
  class(pt) = c("netaffx.db.annot", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = paste("This package contains a detailed transcriptclusterid-level annotation file for use with the GCSscore R package.  It was created by parsing annotation data from the 'netaffxTranscript.rda' data contained within the Bioconductor chip-type platform design (pd) package: ",chip.pd, "and from the Bioconductor transcriptclusterid-level annotation data (db) package: ",pkg.db,".  This package was created using a customized version of the makeProbePackage function and a custom 'ProbePkg-template', both of which are sourced from the AnnotationForge package (version 1.28.0).  These modificed files are contained in the GCSscore R package, which is available on Github and Bioconductor",sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)),
                CHIP.PD    = chip.pd,
                PKG.DB     = pkg.db)
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}
