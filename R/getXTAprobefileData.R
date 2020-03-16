
getXTAprobefileData <- function(arraytype, datafile,
                                pkgname, chip.pd, comparewithcdf = FALSE)
{

  dat <- data.table::fread(datafile)

  # *** To keep as 'factor' (smaller file size, and faster), must add: 'stringsAsFactors'
  pt = data.table::data.table(probesetid = dat$probesetid,
                               locustype = dat$locustype,
                               probesettype = dat$probesettype,
                               transcriptclusterid = dat$transcriptclusterid,
                               category = dat$category,
                               fid = dat$fid,
                               fsetid = dat$fsetid,
                               meta_fsetid = dat$meta_fsetid,
                               x = dat$x,
                               y = dat$y,
                               GC.count = dat$GC.count,
                              stringsAsFactors = TRUE)

  class(pt) = c("probeFile", class(pt))
  
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = paste("This probeFile for use with the GCSscore R package.  It was created by parsing data from the Bioconductor chip-type platform design (pd) package: ", chip.pd,sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)),
                CHIP.PD    = chip.pd)
  
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}
