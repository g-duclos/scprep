![Logo](inst/extdata/scprep_Logo.png)

***

# **scprep**
An R package for aggregating single-cell RNA-Seq data and metadata in an ExpressionSet object, along with biomaRt gene annotation, basic cell filtering and QC metric calculation.

***

## Installation

The **scprep** R package can be installed from Github using devtools:
```
devtools::install_github("g-duclos/scprep")
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
Define the input directory (*dir_input*), which must contain a subdirectory named after each sample. The input directory (*dir_input*) must also contain the following file(s) produced by 10X Genomics' *Cell Ranger* pipeline:

* **filtered_feature_matrix_bc.h5** (the filtered gene counts matrix for each sample)
* If working with 10x Genomics Immune Profiling assay that includes 5' RNA-Seq with TCR or Ig V(D)J data: **filtered_contig_annotations.csv**

#### Dependency Note:
The ["Seurat" R package](https://satijalab.org/seurat/) must be installed (v3, v4, or v5 is acceptable, only the 'Read10X_h5' function is required) in order to use **scprep**. However, Seurat is NOT included as a formal package dependency due to common installation complications.

***

## Overview

Template function to aggregate gene counts matrices from multiple samples, store aggregated counts matrices and metadata (for samples and genes) in an "ExpressionSet" S4 object, add biomaRt gene annotation, perform cell filtering, and calculate select QC metrics.
```
library(Biobase)
dataset <- scprep::template_scprep(dir_output=dir_output)
```

**Core functions:**

#### Read filtered_feature_matrix_bc.h5 file for each sample listed in scprep_annotation.csv into an ExpressionSet object. Add sample metadata from scprep_annotation.csv to the "pData" slot of the ExpressionSet object. Calculate transcripts ("UMIs") per cell and genes ("Genes") per cell (>=1 transcript detected) and add to the "pData" slot of the ExpressionSet object. If working with multi-modal RNA/V(D)J, CITE, or RNA/ATAC data, this function will also store the V(D)J, CITE ADT surface protein, or ATAC information in the ExpressionSet object.</summary>

```
# Build ExpressionSet object with GEX counts and cell metadata
dataset <- scprep::scprep_eset_build(
	sample_paths=sample_paths,
	annotation=annotation,
	vdj=vdj,
	cite=cite,
	atac=atac)
```

<details>
	<summary>Add gene-level metadata to the "fData" slot of the ExpressionSet object (Ensembl ID, gene ID, chromosome #, chromosome start, chromosome stop, biotype) using the <a href="https://bioconductor.org/packages/release/bioc/html/biomaRt.html">biomaRt R package</a>. Calculate transcripts per cell derived from each biotype and chromosome and add to the "pData" slot of the ExpressionSet object.</summary>
<pre>
# Add gene-level metadata to the ExpressionSet object
dataset <- scprep::scprep_eset_biomart(
	dataset=dataset,
	ensembl_target=ensembl_target,
	reference=reference)
</pre>
</details>


<details>
	<summary>Assign barcodes the status of "Dead" if a high percentage of total transcripts are derived from mitochondrial genes ("max_mito" fraction specified in "scprep_parameters.csv"). Assign barcodes the status of "Debris" if low numbers of total transcripts or genes were detected per cell ("min_umi" and "min_gene" specified in "scprep_parameters.csv"). Assign barcodes the status of "Cell" if percent mitochondrial transcripts is less than "max_mito" and if total transcripts is greater than "min_umi".</summary>
<pre>
# Assign status of high quality "Cell", "Dead", or "Debris" to each barcode
dataset$Cell_Filter <- as.factor(scprep::scprep_cell_filter_multi(
	dataset=dataset,
	min_umi=min_umi,
	min_gene=min_gene,
	max_mito=max_mito))
</pre>
</details>


<details>
	<summary>Additional features</summary>
<ul><li>Store parameters specified in the scprep_parameters.csv file in a list labeled "Parameters" in the "assayData" slot labeled "Params" in the ExpressionSet object</li>

<li>Generate random seeds utilized for this analysis and store in a list labeled "Seeds" in the "assayData" slot labeled "Params" in the ExpressionSet object</li>

<li>Select genes with at least 3 transcript counts in a pre-specified (see "gene_filter" in scprep_parameters.csv) percentage (default = 0.1%) of cells and label as "Expressed" in "fData" slot of ExpressionSet object. All other genes are labeled as "Not_Expressed" in "fData" slot of ExpressionSet object</li>

<li>Save ExpressionSet RDS object in *dir_output*</li>
</ul>
</details>

***
