
getXTAprobefileData <- function(arraytype, datafile,
                                pkgname = NULL, comparewithcdf = FALSE)
{
  pkgname <- arraytype
  dat <- data.table::fread(datafile)
  # compare 'pt' as a data.frame vs a data.table:
  # pt = data.frame(probesetid = dat$probesetid,
  #                 locustype = dat$locustype,
  #                 probesettype = dat$probesettype,
  #                 transcriptclusterid = dat$transcriptclusterid,
  #                 category = dat$category,
  #                 fid = dat$fid,
  #                 fsetid = dat$fsetid,
  #                 meta_fsetid = dat$meta_fsetid,
  #                 x = dat$x,
  #                 y = dat$y,
  #                 GC.count = dat$GC.count)
  
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
