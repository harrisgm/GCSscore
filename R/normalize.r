normalize <- function (Score) {
  # function to normalize the GCSscore values:
	x <- Score
	Sx <- sum(x)
	Sxx <- sum(x*x)
	Sstdev <- sqrt((Sxx-Sx*Sx/length(Score))/length(Score))
	# ON: 02.24.20 (Should meanSx be calculated before or after the SD threshold???)
	meanSx <- Sx/length(Score)
	x <- Score-meanSx;
	# x <- x[abs(x) < 3*Sstdev]
	# ON: 02.24.20 (For MTA normalization: try using higher SD threshold)
	# x <- x[abs(x) < 3*Sstdev]
	x <- x[abs(x) < 3*Sstdev]
	Sx <- sum(x)
	Sxx <- sum(x*x)
	num <- length(x)
	Sstdev <- ((Sxx-Sx*Sx/num)/num)
	# Cutoff calculation removed from GCSscore algorithm:
	# if (Sstdev < 0.01) {
	# 	Sstdev <- 1.0
	# } else {
	# 	Sstdev <- sqrt(Sstdev)
	# }
	#
	Sstdev <- sqrt(Sstdev)
	meanSx <- Sx/num+meanSx
	Score <- (Score-meanSx)/Sstdev
	return(Score)
}

