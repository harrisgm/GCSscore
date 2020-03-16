\name{@PKGNAME@}
\docType{data}
\alias{@PKGNAME@}
\title{Parsed gene-level netaffx annotations for chip: @ARRAYTYPE@}
\description{This package contains a data.table object with parsed annotations netaffx annotations from the platform design (pd) package: @CHIP.PD@, and the annotation package: @PKG.DB@. This package utilizes a customized version of the makeProbePackage() function and a custom 'ProbePkg-template', both of which are sourced from the AnnotationForge package (version 1.28.0).}
\usage{@PKGNAME@::@PKGNAME@}
\format{A data.table with @NROW@ rows and @NCOL@ columns, as follows.
  \tabular{lll}{
    transcriptclusterid \tab character \tab Affymetrix transcriptclusterids for given chip-type \cr
    symbol  \tab factor  \tab all gene symbols parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    name  \tab factor  \tab all gene names parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    ref_id  \tab factor  \tab main reference id parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    db.symbol  \tab factor \tab well-annotated gene symbols parsed from annotation (.db) package: @PKG.DB@ \cr
    db.name  \tab factor \tab well-annotated gene names parsed from annotation (.db) package: @PKG.DB@ \cr
    chr  \tab factor \tab chromosome number for transcriptclusterid \cr
    start  \tab integer \tab start location on chromosome for transcriptclusterid \cr
    stop  \tab integer \tab stop location on chromosome for transcriptclusterid \cr
    nProbes  \tab integer \tab number of probes contained within the transcriptclusterid \cr
    locustype  \tab factor \tab transcriptclusterid locus-type parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
    category  \tab factor \tab transcriptclusterid category parsed from netaffxTranscript.rda data in package: @CHIP.PD@ \cr
}}

\source{@DATASOURCE@}

\keyword{Annotation datasets}
\examples{
   @ARRAYTYPE@.TC.netaffx.annot <- @PKGNAME@::@PKGNAME@
}
