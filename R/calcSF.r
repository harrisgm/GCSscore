calcSF <- function(diff, probetab, trim) {
  # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
  probeDiff <- '.' <- NULL
  
	diff[diff <= 0] <- 2^-20
	probetab[,probeDiff := diff]
	tukey <- probetab[,.(SF = tbrm(probeDiff, C = 5)), keyby = key(probetab)]
	SF <- tukey[,500/mean(SF, trim)]
	probetab[,probeDiff := NULL]
	return(SF)
}
