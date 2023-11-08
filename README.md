![Logo](inst/extdata/scprep_Logo.png)

***

# **scprep**
An R package for aggregating single-cell RNA-Seq data and metadata in an ExpressionSet object, along with biomaRt gene annotation and basic cell filtering

***

## Installation

Clone the 'scprep' repository
```
git clone https://github.com/g-duclos/scprep.git
```

Run R and install the package
```
install.packages("scprep", repos=NULL, type="source")
```

***

## Getting Started

#### Annotation
Specify the sample metadata - view the annotation file here: [scprep_annotation.csv](inst/extdata/scprep_annotation.csv) (click "Raw file" at the top right to download)
* Column 1: "Sample_ID" corresponds to the name of each sample
* Column 2: "Index" corresponds to the name of the PCR index used during library preparation
* Column 3: "Sample_Project" corresponds to the name of the project affiliated with all samples
* Column 4: "Reference" corresponds to the name of reference used when running Cell Ranger
* Additional columns can be added corresponding to useful metadata for a particular experiment

#### Parameters
Specify the pipeline parameters - view the parameters file here: [scprep_parameters.csv](inst/extdata/scprep_parameters.csv) (click "Raw file" at the top right to download)

**Critical Parameters**
* *dir_input* (path to input directory)
* *dir_output* (path to output directory)

#### Input Directory
Define the input directory (*dir_input*), which must contain a subdirectory named after each sample. The input directory (*dir_input*) must also contain the following files produced by 10X Genomics' *Cell Ranger* pipeline:

* **filtered_feature_matrix_bc.h5** (the filtered gene counts matrix for each sample)

***

## Overview

Template function to aggregate gene counts matrices from multiple samples, store aggregated counts matrices and metadata (for samples and genes) in an "ExpressionSet" S4 object, perform cell filtering, and calculate select QC metrics.
```
library(Biobase)
dataset <- scprep::template_scprep(dir_output=dir_output)
```

Core functions include:

* Read filtered_feature_matrix_bc.h5 file for each sample listed in scprep_annotation.csv into an ExpressionSet object. Add sample metadata from scprep_annotation.csv to the "pData" slot of the ExpressionSet object. Calculate transcripts ("UMIs") per cell and genes ("Genes") per cell (>=1 transcript detected) and add to the "pData" slot of the ExpressionSet object.

* Add gene-level metadata to the "fData" slot of the ExpressionSet object (Ensembl ID, gene ID, chromosome #, chromosome start, chromosome stop, biotype). Calculate transcripts per cell derived from each biotype and chromosome and add to the "pData" slot of the ExpressionSet object.

* Store parameters specified in the scprep_parameters.csv file in a list labeled "Parameters" in the "assayData" slot labeled "Params" in the ExpressionSet object

* Generate random seeds utilized for this analysis and store in a list labeled "Seeds" in the "assayData" slot labeled "Params" in the ExpressionSet object

* Assign barcodes the status of "Dead" if a high percentage of total transcripts are derived from mitochondrial genes ("max_mito" fraction specified in parameters.csv). Assign barcodes the status of "Debris" if low numbers of total transcripts or genes were detected per cell ("min_umi" and "min_gene" specified in parameters.csv). Assign barcodes the status of "Cell" if percent mitochondrial transcripts is less than "max_mito" and if total transcripts is greater than "min_umi".

* Select genes with at least 3 transcript counts in 0.1% of cells and label as "Expressed" in "fData" slot of ExpressionSet object. All other genes are labeled as "Not_Expressed" in "fData" slot of ExpressionSet object.

* Save ExpressionSet object in *dir_output*

***
