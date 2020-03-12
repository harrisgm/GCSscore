# 03.08.20: create data.table version of: getClariomSprobefileData

getClariomSprobefileData <- function(arraytype, datafile,
                                     pkgname = NULL, comparewithcdf = FALSE)
{
  pkgname <- arraytype
  
  dat <- data.table::fread(datafile)
  
  pt = data.table::data.table(transcriptclusterid = dat$transcriptclusterid,
                              category = dat$category,
                              locustype = dat$locustype,
                              fid = dat$fid,
                              fsetid = dat$fsetid,
                              type = dat$type,
                              x = dat$x,
                              y = dat$y,
                              GC.count = dat$GC.count,
                              stringsAsFactors = TRUE)
  
  class(pt) = c("probeFile", class(pt))
  
  ## assign
  dataEnv = new.env(parent=emptyenv())
  assign(pkgname, pt, envir=dataEnv)
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  # if(is.character(datafile))
  #   datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
  #                      ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}
