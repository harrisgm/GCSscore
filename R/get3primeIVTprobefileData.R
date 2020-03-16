
get3primeIVTprobefileData <- function(arraytype, datafile,
                                      pkgname, chip.pd, comparewithcdf = FALSE)
{

  dat <- data.table::fread(datafile)
  ## data table with the custom probeFile data:
  pt = data.table::data.table(probesetid = dat$man_fsetid,
                              fsetid = dat$fsetid,
                              fid = dat$fid,
                              x = dat$x,
                              y = dat$y,
                              GC.count = dat$GC.count,
                              stringsAsFactors = TRUE)
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = paste("This probeFile for use with the GCSscore R package.  It was created by parsing data from the Bioconductor chip-type platform design (pd) package: ", chip.pd,sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)),
                CHIP.PD    = chip.pd)
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}