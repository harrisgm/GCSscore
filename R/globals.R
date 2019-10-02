# globals.R script 

# used to solve following issues when running 'R CMD check':
  # Undefined global functions or variables:

utils::globalVariables(c("GC.count","MM_fid","fid","fsetid","man_fsetid","meta_fsetid",
                         "netaffxTranscript","netaffxProbeset","locustype","probesetType",
                         "pmSequence","probe_id","probesetid","transcriptclusterid",
                         "type","transcript_cluster_id", "probeDiff","xZone","X","yZone",
                         "Y","Intensities","STDVS","nPixels","rawS","."))




# function with all previously used variables -----------------------------
# utils::globalVariables(c(".lgExtraParanoia","GC.count","MM_fid","fid",
#                          "fsetid","man_fsetid","meta_fsetid","netaffxTranscript",
#                          "netaffxProbeset","locustype","fs1_man_fsetid","probeset_id",
#                          "category","transcriptID","geneName","probe_type","probesetType",
#                          "probesetID","pmSequence","probe_id","probesetid","transcriptclusterid",
#                          "type","transcript_cluster_id", "probeDiff","xZone","X","yZone","Y","Intensities",
#                          "STDVS","nPixels","rawS","."))


# console NOTES if no globals are defined: --------------------------------

# ClariomDXTApFBuilder: no visible binding for global variable
# ‘transcript_cluster_id’
# ClariomDXTApFBuilder: no visible binding for global variable ‘fsetid’
# ClariomDXTApFBuilder: no visible binding for global variable ‘fid’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘pmSequence’
# ClariomDXTApFBuilder: no visible binding for global variable ‘GC.count’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘netaffxTranscript’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘netaffxProbeset’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘transcriptclusterid’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘locustype’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘probesetid’
# ClariomDXTApFBuilder: no visible binding for global variable
# ‘man_fsetid’
# ClariomSpFBuilder: no visible binding for global variable ‘fsetid’
# ClariomSpFBuilder: no visible binding for global variable ‘fid’
# ClariomSpFBuilder: no visible binding for global variable ‘pmSequence’
# ClariomSpFBuilder: no visible binding for global variable ‘GC.count’
# ClariomSpFBuilder: no visible binding for global variable
# ‘netaffxTranscript’
# ClariomSpFBuilder: no visible binding for global variable
# ‘transcriptclusterid’
# ClariomSpFBuilder: no visible binding for global variable ‘man_fsetid’
# GCSscore: no visible binding for global variable ‘probesetid’
# GCSscore: no visible binding for global variable ‘probe_id’
# GCSscore: no visible binding for global variable ‘transcriptclusterid’
# GCSscore: no visible binding for global variable ‘meta_fsetid’
# GCSscore: no visible global function definition for ‘.’
# GCSscore: no visible binding for global variable ‘type’
# GCSscore: no visible binding for global variable ‘fsetid’
# GCSscore: no visible binding for global variable ‘MM_fid’
# GCSscore: no visible binding for global variable ‘fid’
# GCSscore: no visible binding for global variable ‘man_fsetid’
# IVT3primepFBuilder: no visible binding for global variable ‘fsetid’
# IVT3primepFBuilder: no visible binding for global variable ‘fid’
# IVT3primepFBuilder: no visible binding for global variable ‘pmSequence’
# IVT3primepFBuilder: no visible binding for global variable ‘GC.count’
# WTgeneexonpFBuilder: no visible binding for global variable ‘fsetid’
# WTgeneexonpFBuilder: no visible binding for global variable ‘fid’
# WTgeneexonpFBuilder: no visible binding for global variable
# ‘pmSequence’
# WTgeneexonpFBuilder: no visible binding for global variable ‘GC.count’
# calcSF: no visible binding for global variable ‘probeDiff’
# calcSF: no visible global function definition for ‘.’
# computeSscore: no visible binding for global variable ‘GC.count’
# computeSscore: no visible binding for global variable ‘Intensities’
# computeSscore: no visible binding for global variable ‘fid’
# computeSscore: no visible binding for global variable ‘MM_fid’
# computeSscore: no visible binding for global variable ‘rawS’
# computeSscore: no visible global function definition for ‘.’
# mismatch: no visible global function definition for ‘.’
# mismatch: no visible binding for global variable ‘fid’
# mismatch: no visible binding for global variable ‘GC.count’
# zoneRQ: no visible binding for global variable ‘xZone’
# zoneRQ: no visible binding for global variable ‘X’
# zoneRQ: no visible binding for global variable ‘yZone’
# zoneRQ: no visible binding for global variable ‘Y’
# zoneRQ: no visible binding for global variable ‘Intensities’
# zoneRQ: no visible global function definition for ‘.’
# zoneRQ: no visible binding for global variable ‘STDVS’
# zoneRQ: no visible binding for global variable ‘nPixels’
# Undefined global functions or variables:
#   . GC.count Intensities MM_fid STDVS X Y fid fsetid locustype
# man_fsetid meta_fsetid nPixels netaffxProbeset netaffxTranscript
# pmSequence probeDiff probe_id probesetid rawS transcript_cluster_id
# transcriptclusterid type xZone yZone


# List of previous NULLED variables ---------------------------------------

# 
# # ClariomDpFbuilder:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# transcript_cluster_id <- fsetid <- fid <- pmSequence <- GC.count <- netaffxTranscript <- NULL 
# netaffxProbeset <- transcriptclusterid <- locustype <- probesetid <- man_fsetid <- NULL
# .lgExtraParanoia <- NULL
# 
# # main GCSscore:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# fs1_man_fsetid <- transcript_cluster_id <- probeset_id <- transcriptClusterID <- category <- probesetType <- transcriptID <- '.' <- geneName <- probe_type <- probesetID <- NULL
# 
# # calcSF:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# probeDiff <- '.' <- NULL
# 
# # ZoneRQ:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# xZone <- X <- yZone <- Y <- Intensities <- '.' <- STDVS <- nPixels <- NULL
# 
# # mismatch:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# GC.count <- fid <- '.' <- NULL
# 
# # computeSscore:
# # Insert fake-NULL fix to prevent 'no visible binding for global variable' in R CMD check:
# GC.count <- Intensities <- fid <- rawS <- '.' <- NULL
# 









