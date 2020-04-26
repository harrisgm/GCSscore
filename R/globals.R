# globals.R script 

# used to solve following issues when running 'R CMD check':
  # Undefined global functions or variables:

utils::globalVariables(c("GC.count","MM_fid","fid","fsetid","man_fsetid",
                         "meta_fsetid","netaffxTranscript","netaffxProbeset",
                         "locustype","probesetType","pmSequence","probe_id",
                         "probesetid","transcriptclusterid","type",
                         "transcript_cluster_id", "probeDiff","xZone","X",
                         "yZone","Y","Intensities","STDVS",
                         "nPixels","rawS",".","name","chr","symbol","ref_id",
                         "nProbes","i.nProbes","category","probesettype"))
