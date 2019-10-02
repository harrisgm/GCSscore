mismatch <- function(probes, bgp, intensity) {
  # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
  GC.count <- fid <- '.' <- NULL
  
	med.gcList <- bgp[,.(MedianGC = median(intensity[fid])), keyby = .(GC.count)]
	mmIndex <- probes[,as.numeric(GC.count) - (min(as.numeric(GC.count)) - 1)]
	MM <- med.gcList[[2]][mmIndex]
	return(MM)
}
