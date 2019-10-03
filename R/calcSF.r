calcSF <- function(diff, probetab, trim) {
  # function for calculating scaling factors
	diff[diff <= 0] <- 2^-20
	probetab[,probeDiff := diff]
	tukey <- probetab[,.(SF = tbrm(probeDiff, C = 5)), keyby = key(probetab)]
	SF <- tukey[,500/mean(SF, trim)]
	probetab[,probeDiff := NULL]
	return(SF)
}
