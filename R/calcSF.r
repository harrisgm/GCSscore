calcSF <- function(diff, probetab, trim, clean.chip) {
  # function for calculating scaling factors
	diff[diff <= 0] <- 2^-20
	probetab[,probeDiff := diff]
	XTA.chips <- c("mta10","rta10","hta20","clariomdhuman")
	# REAL LINE (02.24.20):
	# EDIT 03.08.20:  
	# ADD CONDITIONAL FOR "XTA" ARRAYS TO USE THE PSR/JUC PROBESETS FOR THE SF CALCULATION
	# INSTEAD OF USING THE "TCID"s TO CALCULATE THE SF
	# THIS WAY THE CHIPS ARE TRULY SCALED BASED ON (NEARLY) ALL OF THE PROBES:
	# if clean.chip is one of the 4 XTA-style arrays:
	if (clean.chip %in% XTA.chips){
	  # Set the key to be "probesetid" regardless of the method.
	  # This will use all PSR/JUC and internal controls that have 4 or more probes:
	  # If we add in the gene/exon arrays, then the exon1.0ST arrays will go here:
	  tukey <- probetab[,.(SF = tbrm(probeDiff, C = 5)), keyby = "probesetid"]
	} else{
	  # If it is not an XTA-style array, then use the regular method:
	  # This includes the ClariomS arrays, the 3' IVT arrays, (and potentially the gene1.0ST arrays):
	tukey <- probetab[,.(SF = tbrm(probeDiff, C = 5)), keyby = key(probetab)]
	}
	
	SF <- tukey[,500/mean(SF, trim)]
	probetab[,probeDiff := NULL]
	return(SF)
}
