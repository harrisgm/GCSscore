mismatch <- function(probes, bgp, intensity) {
  # Function for calculating GC-based bkg correction:
	med.gcList <- bgp[,.(MedianGC = median(intensity[fid])), keyby = .(GC.count)]
	mmIndex <- probes[,as.numeric(GC.count) - (min(as.numeric(GC.count)) - 1)]
	MM <- med.gcList[[2]][mmIndex]
	return(MM)
}
