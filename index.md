# <img src="man/figures/SPARROW_logo.png" width="40" style="vertical-align: middle; margin-right: 10px;" /> SPARROW

**The sparrow may be small, but it is fully equipped with all vital organs...**

SPARROW is a method for subsampling large scale data to obtain a lightweight presentation such that most of the statistical power is preserved with respect to downstream analysis.  

The name "sparrow" comes from the Chinese idiom, ["the sparrow may be small, but it is fully equiped with all of the vital organs"](https://en.wiktionary.org/wiki/%E9%BA%BB%E9%9B%80%E9%9B%96%E5%B0%8F%EF%BC%8C%E4%BA%94%E8%87%9F%E4%BF%B1%E5%85%A8#Chinese).  The letters in the acronym stand for: **S**ubmodular **P**ower-**A**daptive **R**eduction for **R**epresentative D**o**wnstream **W**orkflows.

## Installation

To install SPARROW from GitHub:


```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("kaishumason/SPARROW")
```

SPARROW can be paired with [spotGLM](http://kaishumason.github.io/SpotGLM/) for fast and scalable spatial omic data analysis.  
We also recommend for you to install spotGLM and check out the [tutorials](https://kaishumason.github.io/SpotGLM/index.html) from there:

```r
# Install from GitHub
devtools::install_github("kaishumason/spotGLM")
```


## Getting Started

Please follow these tutorials to get started on your data:

[Demonstration on a MERFISH dataset](articles/Subsampling_Data_Using_SPARROW.html)

[Demonstration on a simulated spatial-barcoding data set that requires deconvolution](articles/Subsampling_Data_Using_SPARROW.html)

[Demonstration on a Visium HD data set](https://kaishumason.github.io/SpotGLM/articles/Vignette_VisiumHD_Mouse_Kidney_analysis.html)


SPARROW was originally designed for large scale spatial transcriptomics data analysis, to be used in conjunction with [the spotGLM package](https://kaishumason.github.io/spotGLM/) for identifying spatial signals.  However, it can be used as a general purpose data selection tool.   

