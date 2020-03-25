
THE `GCSscore` PACKAGE
===================
`GCSscore` is an R package that is available on GitHub that is also available in the BioConductor repository (release 3.10 and higher).  The `GCSscore` package is used for the differential expression analysis of Affymetrix microarrays. In the current release, the GCSscore package fully supports all ClariomS and Clariom/XTA chips, and it has full support for the following two 3' IVT GeneChips: the Mouse Genome 430 2.0 Array, and the Human Genome U133 Plus 2.0 Array.  In general, the package reads probe level information directly from array files in their native format (.CEL files), and then performs calculations on individual probe comparisons between arrays without the need for prior probe grouping summarization and normalizations steps.

The `GCSscore` package builds and installs its own packages for probe-specific data and up-to-date annotations, using source information by using the platform design (.pd) and annotation (.db) packages found on the BioConductor repository.  Therefore, the user must have require the ability to compile packages from source, even if the package is installed from BioConductor.

Installation
------------
The `GCSscore` package builds and installs its own packages for probe-specific data and up-to-date annotations, using source information by using the platform design (.pd) and annotation (.db) packages found on the BioConductor repository.  Therefore, the user must have require the ability to compile packages from source, even if the package is installed from BioConductor.

The package depends on the "devtools" package to install it directly from github and build the parsed probe-level and annotation packages and to install this package directly from github.  More information regarding the installation and updating of the devtools package can be found in the "[Updating to the latest version of devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)" section of the package documentation.
 
GCSscore requries the ability to compile R packages, and so it de
       
For macOS users, command line tools needs to be installed.  In terminal, run the following command: 
		
	xcode-select --install
        
For windows users, install [Rtoools35.exe](https://cran.r-project.org/bin/windows/Rtools/) and ensure that it is added to the PATH               

Install GCSscore package directly from github using 'devtools': will use "devtools" to build all of the packages from source on the fly.


```r
devtools::install_github("harrisgm/GCSscore")
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
