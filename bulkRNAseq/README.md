# Bulk RNA-seq Analysis

This folder contains the code to reproduce the figures related to bulk RNA-seq data from [Modeling drug responses and evolutionary dynamics using triple negative breast cancer patient-derived xenografts](https://doi.org/10.1101/2023.01.10.523259).

## Setup
Instruction to setup are included in the R notebook, we recommend using [`mamba`](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to install all required packages. We provided a yaml file to reproduce our installation:

```{bash}
# Create environment
mamba env create --file yaml/bulkrnaseqpdx.yml
# Activate environment
mamba activate bulkrnapdx
# Re-run this R notebook to obtain HTML or open in Rstudio to work on it interactively
Rscript -e "rmarkdown::render('BulkRNASeq_Analysis.Rmd')"
```

## Data

All necessary data can be folder in the `data` folder
