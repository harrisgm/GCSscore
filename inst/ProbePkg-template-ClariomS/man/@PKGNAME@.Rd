\name{@PKGNAME@}
\docType{data}
\alias{@PKGNAME@}
\title{Parsed Affymetrix probe-level information for chip: @ARRAYTYPE@}
\description{This package was created from probe-level data in the platform design package: @CHIP.PD@.  This package utilizes a customized version of the makeProbePackage() function and a custom 'ProbePkg-template', both of which are published the AnnotationForge package (version 1.28.0).}
\usage{@PKGNAME@::@PKGNAME@}
\format{A data.table with @NROW@ rows and @NCOL@ columns, as follows.
  \tabular{lll}{
    transcriptclusterid \tab character \tab Affymetrix transcriptclusterids for chip-type: @ARRAYTYPE@ \cr
    category  \tab factor \tab transcriptclusterid category parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    locustype  \tab factor \tab transcriptclusterid locus-type parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    fid  \tab integer  \tab probe index on array: @ARRAYTYPE@ \cr
    fsetid  \tab integer  \tab featureset number for probe on array: @ARRAYTYPE@ \cr
    type  \tab factor \tab probe type identification number for probe on array: @ARRAYTYPE@ \cr
    x  \tab integer  \tab x-coordinate on array: @ARRAYTYPE@ \cr
    y  \tab integer  \tab y-coordinate on array: @ARRAYTYPE@ \cr
    GC.count  \tab integer  \tab GC-content of probe on array: @ARRAYTYPE@ \cr
}}

\source{@DATASOURCE@}

\keyword{Annotation datasets}
\examples{
   @ARRAYTYPE@.probeFile <- @PKGNAME@::@PKGNAME@
}
