zoneRQ <- function(DT, affyCel, trim) {
  # function to calculate the RawQ of 16 sectored zones
	rows <- affyCel$header$rows; cols <- affyCel$header$cols
	DT[,xZone := cut(X, seq(0, cols, by = cols/4), dig.lab = 10, include.lowest = TRUE)]
	DT[,yZone := cut(Y, seq(0, rows, by = rows/4), dig.lab = 10, include.lowest = TRUE)]
	sectorRQ <- setkey(DT, Intensities)[,.(RawQ = sum(STDVS[1:(trim * .N)] /
					sqrt(nPixels[1:(trim * .N)]))/(trim * .N)), by = .(xZone, yZone)]
	rq <- mean(sectorRQ[[3]])
	return(rq)
}
