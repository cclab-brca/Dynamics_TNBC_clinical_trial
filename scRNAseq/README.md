# scRNA-seq analysis

This folder contains the code and processed data required to reproduce the scRNA-seq analysis done in this study. 

The code begins with aggregating individual count matrices from each sample into a single matrix, processing it and generate the plots. Intermediate processed files are supplied for user convenience. In case the user wishes to reproduce the processing from scratch these intermediate files can be overwritten by setting the *rebuild* parameter in **PAR1006_scRNAseq_analysis.R** to TRUE. 

## Requirements
The data was processed and analysed with the [*metacell*](https://github.com/tanaylab/metacell) R package. It is required to reproduce the processing and for plot generation.

Additional required R packages:

* data.table
* dplyr
* ggplot2
* ggpubr
* glue
* MASS
* Matrix
* officer
* pheatmap
* R.utils
* RColorBrewer
* readr
* rstatix
* tgconfig
* tgstat
* tgutil

We use [*scrublet*](https://github.com/AllonKleinLab/scrublet) to detect and remove doublets. If you wish to reproduce that step (by setting the rebuild variable) you should also install:

* reticulate
* scrublet under a virtual environment named *r-scrublet*

## Supplied data

The *data* folder contains the following:

* Auxiliary files used for the analysis (sample information, list of TFs, gene information) 
* Processed count matrices per sample (output of 10X cellranger after further QC filtering). These are under the scdb folder and their names begins with mat.#_PARA1006...Rda
* Intermediate/final processed files are also found under the scdb folder. All files under the scdb folder will be copied to a new scdb folder under the output directory, which will be the working directory. 
