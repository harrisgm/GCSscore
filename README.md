
THE `GCSscore` PACKAGE
===================

The `GCSscore` package is used for the differential expression analysis of Affymetrix microarrays. It is designed to be able to support all Affymetrix chips.  It functions by reading information directly from array files in their native format (.CEL files), and it performs calculations on individual probe comparisons between arrays without the need for prior probe grouping summarization and normalizations steps.

Installation
------------

`GCSscore` is a github package that is in the process being uploaded to the BioConductor repository. Until it is uploaded to Bioconductor, the recommended way to install it is to load `R` and ensure that all dependencies are install prior to install the `GCSscore` package from source, as shown below:

Dependencies from CRAN (run commands in Rstudio/R.app):

```r
install.packages("data.table")
install.packages("R.utils")
install.packages("dplR")
```

Dependencies from Bioconductor (run commands in Rstudio/R.app):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("affxparser")
```

Once the dependencies are installed, install the `GCSscore` package from source.

NOTE: the user must have the `GCSscore_1.0.0.tar.gz` file in the current working directory or the user must provide the path to the file location.  Alternatively, the user can use the package install GUI.

GUI use in Rstudio (Tools -> Install Packages --> select 'Install from: Package Archive File (.zip; .tar.gz)' --> use 'Browse' to navigate to the source file location.

(OR, run command in Rstudio/R.app):

```r
install.packages("GCSscore_1.0.0.tar.gz", repos = NULL, type = "source")
```

What does the package offer?
-------------------
The GC-Sscore algorithm is a comparative method for gene expression data analysis that performs tests of hypotheses directly from probe level data. It is based on an error model in which the detected signal is assumed to be proportional to the probe signal for highly expressed genes, but assumed to approach a background level (rather than 0) for genes with low levels of expression. This error model is used to calculate relative changes in probe intensities that converts probe signals into multiple measurements with equalized errors, which are summed over a probe set to form the significance score (S-score).  The original S-score method required the mismatch (MM) probes to estimate non-specific binding (NSB) for each perfect-match (PM) probes, and the MM probes were removed the arrays beginning with the Affymetrix Whole Transcriptome (WT) style arrays. This new algorithm uses a gc-content based NSB, thus eliminating the original algorithm's dependence on MM probes.  The S-score2 algorithm is capable of working on all modern Affymetrix array types (3' IVT and up). Assuming no expression differences between chips, the S-score2 output follows a standard normal distribution. Thus, a separate step estimating the probe set expression summary values is not needed and p-values can be easily calculated from the S-score2 output.

Contributing
------------

1. Fork it.
2. Create a branch (`git checkout -b my_contrib`)
3. Commit your changes (`git commit -am "My Contributions"`)
4. Push to the branch (`git push origin my_contrib`)
5. Create an [Issue][1] with a link to your branch

Contact
-------

Guy Harris, M.S.
<harrisgm@vcu.edu>

Michael F. Miles, M.D., Ph.D.
<Michael.Miles@vcuhealth.org>

[1]: https://github.com/harrisgm/GCSscore
