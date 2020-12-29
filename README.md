## CFTK

CFTK is a toolkit for the analysis of co-fractionation mass spectrometry (CF-MS) data. It implements many common tasks for CF-MS data analysis, including preprocessing of CF-MS datasets, evaluating data quality, and inference of protein-protein interaction networks. To learn more about the package, download and open the "Introduction to CFTK" vignette at [https://github.com/fosterlab/CFTK/tree/master/vignettes](https://github.com/fosterlab/CFTK/tree/master/vignettes). 

### System requirements

CFTK relies on functions from the following R packages:

```
	dplyr (>= 0.7.4),
	tidyr (>= 1.1.2),
	readr (>= 1.3.1),
	purrr (>= 0.3.4),
	stringr (>= 1.4.0),
	magrittr (>= 1.5),
	arules (>= 1.6-1),
	gtools (>= 3.5.0),
	lsa (>= 0.73.1),
	pbapply (>= 1.3.4),
	propr (>= 3.1.8),
	treeClust (>= 1.1-7),
	vegan (>= 2.5-6),
	wccsom (>= 1.2.11),
	GENIE3 (>= 1.6.0),
	Pigengene (>= 1.12.0),
	WGCNA (>= 1.51),
	pcaPP (>= 1.9-73),
	AUC (>= 0.3.0),
	reshape2 (>= 1.4.3),
	naivebayes (>= 0.9.1),
	ranger (>= 0.8.0),
	LiblineaR (>= 2.10-8),
	speedglm (>= 0.3-2),
	tester (>= 0.1.7)
```

CFTK has been tested with R version 3.6.0 and higher.

### Installation

To install CFTK, first install the devtools package, if it is not already installed: 

```r
> install.packages("devtools") 
```

Then, install the required dependencies from [Bioconductor](https://www.bioconductor.org/):

```r
> if (!requireNamespace("BiocManager", quietly = TRUE))
>     install.packages("BiocManager")
> BiocManager::install(c("AnnotationDbi", "GO.db", "preprocessCore", "impute", "Pigengene", "GENIE3"))
```

Also install the `kohonen` package, which is required by `wccsom`: 

```r
> install.packages("kohonen")
```

Next, the `wccsom` package, which is no longer distributed via CRAN, must be manually installed. To do so, download and extract the `wccsom` package from [https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz](https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz). For example, on Unix systems, this can be done using `wget` by entering the following commands into a terminal: 

```
wget https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz
tar -xzvf wccsom_1.2.11.tar.gz
```

Then, open R and install `wccsom` using devtools:

```r
> devtools::install("wccsom")
```

On a Windows machine, installing `wccsom` is slightly more complicated. First, navigate to [https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz](https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz) using a web browser in order to download the `wccsom` package. Then, download [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) and follow the installation instructions. Last, from R, install `wccsom` by entering the following command, replacing the string `MyUserName` with your own username, or updating the path to the `wccsom` package as needed. 

```r
> install.packages("C:\\Users\\MyUserName\\Downloads\\wccsom_1.2.11.tar.gz", repos = NULL, type = "source")
```

Finally, install CFTK from GitHub using devtools:

```r
> devtools::install_github("fosterlab/CFTK")
```

This should take no more than a few minutes.

## Usage

For an in-depth guide to the functionality of CFTK, download the "Introduction to CFTK" vignette (file `intro-to-CFTK.html`) at [https://github.com/fosterlab/CFTK/tree/master/vignettes](https://github.com/fosterlab/CFTK/tree/master/vignettes), and open this HTML file in a web browser. The vignette demonstrates the application of CFTK to carry out several exemplary analyses in a published CF-MS dataset. This data is available from the `data-raw` directory (at [https://github.com/fosterlab/CFTK/tree/master/data-raw](https://github.com/fosterlab/CFTK/tree/master/data-raw)). 
