\name{@PKGNAME@}
\docType{data}
\alias{@PKGNAME@}
\title{Parsed exon-level netaffx annotations for chip: @ARRAYTYPE@}
\description{This package contains a data.table object with parsed netaffx annotations from the following Bioconductor annotation packages: (1) platform design (pd) package: @CHIP.PD@, (2) annotation (db) package: @PKG.DB@. This package utilizes a customized version of the makeProbePackage() function and a customized 'ProbePkg-template', both of which are sourced from the AnnotationForge package (version 1.28.0).}
\usage{@PKGNAME@::@PKGNAME@}
\format{A data.table with @NROW@ rows and @NCOL@ columns, as follows.
  \tabular{lll}{
    probesetid \tab character \tab Affymetrix probesetids for given chip-type \cr
    symbol  \tab factor  \tab all gene symbols parsed from netaffxProbeset.rda data in package: @CHIP.PD@ \cr
    ref_id  \tab factor  \tab main reference id parsed from netaffxProbeset.rda data in package: @CHIP.PD@ \cr
    db.symbol  \tab factor \tab well-annotated gene symbols parsed from annotation (.db) package: @PKG.DB@ \cr
    db.name  \tab factor \tab well-annotated gene names parsed from annotation (.db) package: @PKG.DB@ \cr
    chr  \tab factor \tab chromosome number for probesetid \cr
    start  \tab integer \tab start location on chromosome for probesetid \cr
    stop  \tab integer \tab stop location on chromosome for probesetid \cr
    nProbes  \tab integer \tab number of probes contained within the probesetid \cr
    locustype  \tab factor \tab probesetid locus-type information parsed from netaffxProbeset.rda data in package: @CHIP.PD@ \cr
    probesettype  \tab factor \tab probesetid probeset-type information parsed from netaffxProbeset.rda data in package: @CHIP.PD@ \cr
}}

\source{@DATASOURCE@}

\keyword{Annotation datasets}
\examples{
   @ARRAYTYPE@.PSR.netaffx.annot <- @PKGNAME@::@PKGNAME@
}
