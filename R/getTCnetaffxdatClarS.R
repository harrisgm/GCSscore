# Build netaffx annotation package using makeProbePackage:
getTCnetaffxdatClarS <- function(arraytype, datafile,
                                 pkgname=NULL, comparewithcdf = FALSE){
  
  pkgname <- arraytype
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
  
  datasource = "This probeFile for use with the GCSscore R package.  It was created by parsing the data contained within the chiptype platform design (pd.packge) on Bioconductor."
  # if(is.character(datafile))
  # datasource = paste(datasource, " The file name was ", gsub("_", "\\\\_", datafile),
  # ".", sep="")
  
  symVal = list(ARRAYTYPE  = arraytype,
                DATASOURCE = datasource,
                NROW       = as.character(nrow(pt)),
                NCOL       = as.character(ncol(pt)))
  
  # if(comparewithcdf) .lgExtraParanoia(pt, cdfname)
  
  # return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
  return(list(pkgname = pkgname, symVal = symVal, dataEnv = dataEnv))
}