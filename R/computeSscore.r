computeSscore <- function(cel1, cel2, probeFile, bgp, method, 
                          infoKey, SF1 = NULL, SF2 = NULL, 
                          verbose = FALSE, trim = NULL, clean.chip = NULL) {
  
  # Define inverse operators:
  `%!in%` <- Negate(`%in%`)
  `%!like%` <- Negate(`%like%`)
   #  SET LIST OF XTA-STYLE CHIPS:
  XTA.chips <- c("mta10","rta10","hta20","clariomdhuman")
  
  
  #CREATE DATA TABLES FROM CEL FILES
	celdat1 <- data.table(X = cel1$x, Y = cel1$y,
	                      Intensities = cel1$intensities,
	                      STDVS = cel1$stdvs, nPixels = cel1$pixels)
	celdat2 <- data.table(X = cel2$x, Y = cel2$y,
	                      Intensities = cel2$intensities,
	                      STDVS = cel2$stdvs, nPixels = cel2$pixels)

	# OTHER INTENSITY TRANSFORMATION IDEA:
	    # Untransformed intensities does seem to be the best bet
  # Try using pre-normalized intensity values from oligo:
	# celdat1 <- data.table(X = cel1$x, Y = cel1$y,
	#                       # Intensities = cel1$intensities,
	#                       Intensities = (normData@assayData[["exprs"]])[,colnames(normData@assayData[["exprs"]])==celfile1],
	#                       STDVS = cel1$stdvs, nPixels = cel1$pixels)
	# 
	# celdat2 <- data.table(X = cel2$x, Y = cel2$y,
	#                       # Intensities = cel2$intensities,
	#                       Intensities = (normData@assayData[["exprs"]])[,colnames(normData@assayData[["exprs"]])==celfile2],
	#                       STDVS = cel2$stdvs, nPixels = cel2$pixels)
	# 
  # Maybe the problem the whole time has been the lack of the log2 transformation:
	# celdat1 <- data.table(X = cel1$x, Y = cel1$y,
	#                       Intensities = log2(cel1$intensities),
	#                       STDVS = log2(cel1$stdvs), nPixels = cel1$pixels)
	# celdat2 <- data.table(X = cel2$x, Y = cel2$y,
	#                       Intensities = log2(cel2$intensities),
	#                       STDVS = log2(cel2$stdvs), nPixels = cel2$pixels)

	
	
	# ############## Quantile Normalize the CEL intensities 
	#   # This method will use the whole chip to perform the normlization step:
	# celdat <- data.table(cF1 = celdat1$Intensities,
	#                      cf2 = celdat2$Intensities)
	# message("performing quantile normlization of CEL files (using whole chip)")
	# celdat.norm <- as.data.table(normalize.quantiles(as.matrix(celdat)))
	# # placing the normalized values back into celdat intensity values:
	# celdat1[,Intensities := celdat.norm$V1]
	# celdat2[,Intensities := celdat.norm$V2]
	# ##############
	
	#CREATE PROBETAB OBJECT FOR CONSISTENCY
	if (!any(identical(method, "pmmm"), identical(method, "gc"))) probeFile <- probeFile[GC.count < bgp[,min(GC.count)], GC.count := bgp[,min(GC.count)]]


	#SCALING FACTOR CEL 1
	message("calculating scaling factors (SF1 & SF2)")
	if (is.null(SF1)) {
		if (identical(method, "pmmm")) {
			# SF1 <- calcSF(celdat1[,Intensities][probeFile[,fid]] - celdat1[,Intensities][bgp[,fid]], probeFile, trim)
			SF1 <- calcSF(celdat1[,Intensities][probeFile[,fid]] - celdat1[,Intensities][probeFile[,MM_fid]], probeFile, trim, clean.chip)
		} else {
			SF1 <- celdat1[,calcSF(Intensities[probeFile[,fid]] - mismatch(probeFile, bgp, Intensities), probeFile, trim, clean.chip)]
		}
	} else if ((!is.null(SF1)) & SF1 > 0) {
	    SF1 <- SF1
	}
	
	#SCALING FACTOR CEL 2
	if (is.null(SF2)) {
		if (identical(method, "pmmm")) {
			SF2 <- calcSF(celdat2[,Intensities][probeFile[,fid]] - celdat2[,Intensities][probeFile[,MM_fid]], probeFile, trim, clean.chip)
		} else {
			SF2 <- celdat2[,calcSF(Intensities[probeFile[,fid]] - mismatch(probeFile, bgp, Intensities), probeFile, trim, clean.chip)]
		}
	} else if ((!is.null(SF2)) & SF2 > 0) {
	    SF2 <- SF2
	}
	if (verbose) {
		message(paste("SF1:", SF1))
		message(paste("SF2:", SF2))
	}
	
	#INTENSITY ADJUSTMENTS
	intense1 <- celdat1[,Intensities * SF1]
	intense2 <- celdat2[,Intensities * SF2]
	
	#RAWQ VALUES AFTER ZONE DIVISION
	message("calculating rawQ and SDT values")
	zonerq1 <- zoneRQ(celdat1, cel1, trim)
	zonerq2 <- zoneRQ(celdat2, cel2, trim)
	if (verbose) {
		message(paste("rawQ1:", zonerq1))
		message(paste("rawQ2:", zonerq2))
	}
	
	#SDT CALCULATIONS USING AVERAGE RAWQ
	SDT1 <- 4*zonerq1*SF1
	SDT2 <- 4*zonerq2*SF2
	if (verbose) {
		message(paste("SDT1:", SDT1))
	 	message(paste("SDT2:", SDT2))
	}
	
	# EDIT 03.08.20:
	# ALL OF THE ABOVE SHOULD BE THE SAME AND IT SHOULD INCLUDE ALL OF THE PROBESETIDS (>= 4 PROBES), 
	# SINCE THAT ACCOUNTS FOR THE ALMOST ALL OF PROBES ON THE CHIP:
	
	# IF CHIP-TYPE IS XTA-STYLE AND METHOD = TCID,
	# THEN REMOVE THE JUNCTION PROBES FROM THE PROBEFILE, SINCE THERE ARE PURELY FOR EXON-LEVEL ANALYSIS:
	if (clean.chip %in% XTA.chips & key(probeFile) == "transcriptclusterid"){
	  # remove probes with a 'probesetid' with a 'JUC' designation:
	  probeFile <- probeFile[probesetid %!like% "JUC"]
	}
	#CALCULATE PM/MM VALUES
	if (method[1] == "pmmm") {
		#SUBSET INTENSITY VALUES BASED ON PM/MM LISTS
		diff1 <- intense1[probeFile[,fid]] - intense1[probeFile[,MM_fid]]
		diff2 <- intense2[probeFile[,fid]] - intense2[probeFile[,MM_fid]]
	} else {
		diff1 <- probeFile[,intense1[fid] - mismatch(probeFile, bgp, intense1)]
		diff2 <- probeFile[,intense2[fid] - mismatch(probeFile, bgp, intense2)]
	}
	
	#CALCULATE SSCORE
	message("calculating and normalizing GCS-score results")
	probeFile[,rawS := rawScore(diff1, diff2, SDT1, SDT2)]
	# Score <- probeFile[,.(Score = sum(rawS/sqrt(.N))), keyby = key(probeFile)]
	# Score <- probeFile[,.(Score = sum(rawS/sqrt(.N))), keyby = key(info)]
	Score <- probeFile[,.(Score = sum(rawS/sqrt(.N))), keyby = infoKey]
	Score[,Score := normalize(Score)]
	
	if (verbose) {
	  	celName1 <- strsplit(cel1$header$filename, "/")[[1]]; celName1 <- celName1[length(celName1)]
	  	celName2 <- strsplit(cel2$header$filename, "/")[[1]]; celName2 <- celName2[length(celName2)]
	  	  # hist(Score$Score, breaks = seq(-max(abs(Score$Score)),max(abs(Score$Score)),0.1), xlim = c(-5,5), 
	  	  # main = paste(celName1, "vs", celName2),xlab = "S-scores \n (normalized)")
	}

	probeFile[,rawS := NULL]
	return(Score)
}
