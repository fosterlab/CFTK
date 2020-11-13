## CFTK

CFTK is a toolkit for the analysis of co-fractionation mass spectrometry (CF-MS) data. It implements many common tasks for CF-MS data analysis, including preprocessing of CF-MS datasets, evaluating data quality, and inference of protein-protein interaction networks. To learn more about the package, open the "Introduction to CFTK" vignette at [https://github.com/fosterlab/CFTK/tree/master/vignettes](https://github.com/fosterlab/CFTK/tree/master/vignettes). 

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
	GENIE3 (>= 1.8.0),
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

Additionally, the `preprocessCore` package must be installed from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html) in order to perform quantile normalization.

### Installation

To install CFTK, first install the devtools package, if it is not already installed: 

```r
> install.packages("devtools") 
```

Last, install CFTK from GitHub using devtools:

```r
> devtools::install_github("fosterlab/CFTK")
```

This should take no more than a few minutes.

## Usage

For an in-depth guide to the functionality of CFTK, see the "Introduction to CFTK" vignette at [https://github.com/fosterlab/CFTK/tree/master/vignettes](https://github.com/fosterlab/CFTK/tree/master/vignettes). 
