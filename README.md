![Logo](inst/extdata/scprep_Logo_v2.png)

***

# **scprep**
An R package for aggregating single-cell RNA-Seq data and metadata in ExpressionSet, Seurat, or SingleCellExperiment objects, along with biomaRt gene annotation, basic cell filtering and QC metric calculation.

Data that may be used for testing **scprep** may be accessed via the **[scdata](https://github.com/g-duclos/scdata)** package.

For post-processing & secondary analysis of data that has been ingested by **scprep**, this package may be coupled with **[scpost](https://github.com/g-duclos/scpost)**.

***

## Installation

### Option 1: Install from GitHub
The **scprep** R package can be installed from Github using devtools:
```r
devtools::install_github("g-duclos/scprep")
```

### Option 2: Docker Container (Recommended)
For a containerized environment with all dependencies pre-installed:

**Quick Start with Docker:**
```bash
# Build the Docker image
docker build -t scprep .

# Run interactive R session with scprep
docker run -it --rm scprep

# Run with data directory mounted
docker run -it --rm -v /path/to/your/data:/scprep/data scprep
```

**Using Docker Compose:**
```bash
# Start the scprep container
docker-compose up scprep

# For RStudio interface (optional)
docker-compose up scprep-rstudio
# Then navigate to http://localhost:8788
```

**Docker Image Features:**
- Pre-installed R 4.2.1+ with all required dependencies
- Seurat v5, Matrix, Biobase, biomaRt, and other dependencies configured
- Optimized for single-cell RNA-Seq workflows
- Volume mounting for data input/output

***

## Documentation

📖 **[Complete Usage Guide & Object Types Vignette](https://htmlpreview.github.io/?https://github.com/g-duclos/scprep/blob/main/vignettes/vignette_scprep_usage.html)** - Comprehensive guide demonstrating how to read 10X data, convert between object types, and explore data structure

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
* *file_type* (input format: "h5" or "mtx")
* *output_type* (object type: "eset", "seurat", or "sce")
* *gene_id* (gene identifier type: "ensembl" or "symbol")

#### Input Directory
Define the input directory (*dir_input*), which must contain a subdirectory named after each sample. The input directory (*dir_input*) must also contain the following file(s) produced by 10X Genomics' *Cell Ranger* pipeline:

**Option 1: H5 Format (Recommended)**
* **filtered_feature_bc_matrix.h5** (the filtered gene counts matrix for each sample)

**Option 2: MTX Format**
* **matrix.mtx** (gene expression count matrix in Matrix Market format)
* **barcodes.tsv** (cell barcodes)
* **genes.tsv** or **features.tsv** (gene identifiers and names)

**Additional Files:**
* If working with 10x Genomics Immune Profiling assay that includes 5' RNA-Seq with TCR or Ig V(D)J data: **filtered_contig_annotations.csv**

**Format Selection:**
Set the `file_type` parameter in `scprep_parameters.csv`:
- `file_type="h5"` for H5 format (supports multimodal data like CITE-seq/ATAC-seq)
- `file_type="mtx"` for MTX format (RNA expression data only)

#### Dependencies:
**scprep** requires the following R packages:
- **Seurat** (>= 5.0.0): For reading 10X Genomics data formats and Seurat object creation
- **Matrix** (>= 1.2-0): For sparse matrix operations  
- **Biobase**: For ExpressionSet data structures
- **SingleCellExperiment**: For SingleCellExperiment object creation
- **R** (>= 4.0.0): Minimum R version

**Note:** All dependencies are automatically installed when using `devtools::install_github()` or the Docker container.

***

## Output Object Types

**scprep** now supports three different output object types:

### **ExpressionSet** (Default)
- Traditional Bioconductor S4 object for storing expression data
- Metadata stored in `pData()` and `fData()` slots
- Compatible with Bioconductor workflows
- Maintains full backward compatibility

### **Seurat** (v5)
- Popular single-cell analysis framework object
- Metadata stored in `meta.data` slot
- Ready for Seurat downstream analysis workflows
- Supports multi-modal data (protein, VDJ)

### **SingleCellExperiment**
- Modern Bioconductor S4 object for single-cell data
- Metadata stored in `colData()` slot
- Compatible with Bioconductor/scater workflows
- Supports alternative experiments (protein data)

***

## Overview

Template function to aggregate gene counts matrices from multiple samples, store aggregated counts matrices and metadata (for samples and genes) in an ExpressionSet, Seurat, or SingleCellExperiment object, add biomaRt gene annotation, perform cell filtering, and calculate select QC metrics.

**Output type is specified in the `scprep_parameters.csv` file:**
```r
library(Biobase)

# Output type determined by output_type parameter in scprep_parameters.csv
dataset <- scprep::template_scprep(dir_output=dir_output)

```

**Core functions:**

<details>
	<summary>Read filtered_feature_matrix_bc.h5 file for each sample listed in scprep_annotation.csv into an ExpressionSet, Seurat, or SingleCellExperiment object. Add sample metadata from scprep_annotation.csv to the object metadata. Calculate transcripts ("UMIs") per cell and genes ("Genes") per cell (>=1 transcript detected) and add to the object metadata. If working with multi-modal RNA/V(D)J, CITE, or RNA/ATAC data, this function will also store the V(D)J, CITE ADT surface protein, or ATAC information in the object.</summary>
<pre>
# Build Seurat object with RNA data (no additional modalities)
dataset <- scprep::scprep_build(
	sample_paths="path/to/sample",
	annotation=annotation,
	gene_id="symbol"
	output_type="seurat",
	vdj=FALSE,
	cite=FALSE,
	cite_ignore=TRUE,
	atac=FALSE,
	atac_ignore=TRUE)
</pre>
</details>


<details>
	<summary>Add gene-level metadata to the object (Ensembl ID, gene ID, chromosome #, chromosome start, chromosome stop, biotype) using the <a href="https://bioconductor.org/packages/release/bioc/html/biomaRt.html">biomaRt R package</a>. Calculate transcripts per cell derived from each biotype and chromosome and add to the object's cell-level metadata slot.</summary>
<pre>
# Add gene-level metadata to the object
dataset <- scprep::scprep_biomart(
	dataset=dataset,
	ensembl_target=ensembl_target,
	reference=reference)
</pre>
</details>


<details>
	<summary>Assign barcodes the status of "Dead" if a high percentage of total transcripts are derived from mitochondrial genes ("max_mito" fraction specified in "scprep_parameters.csv"). Assign barcodes the status of "Debris" if low numbers of total transcripts or genes were detected per cell ("min_umi" and "min_gene" specified in "scprep_parameters.csv"). Assign barcodes the status of "Cell" if percent mitochondrial transcripts is less than "max_mito" and if total transcripts is greater than "min_umi".</summary>
<pre>
# Assign status of high quality "Cell", "Dead", or "Debris" to each barcode and add to object metadata
dataset <- scprep::scprep_cell_filter_multi(
	dataset=dataset,
	min_umi=min_umi,
	min_gene=min_gene,
	max_mito=max_mito)
</pre>
</details>


<details>
	<summary>Additional features</summary>
<ul><li>Store parameters specified in the scprep_parameters.csv file in a list labeled "Parameters" in the "assayData" slot labeled "Params" in the ExpressionSet object</li>

<li>Generate random seeds utilized for this analysis and store in a list labeled "Seeds" in the "assayData" slot labeled "Params" in the ExpressionSet object</li>

<li>Select genes with at least 3 transcript counts in a pre-specified (see "gene_filter" in scprep_parameters.csv) percentage (default = 0.1%) of cells and label as "Expressed" in "fData" slot of ExpressionSet object. All other genes are labeled as "Not_Expressed" in "fData" slot of ExpressionSet object</li>

<li>Save ExpressionSet RDS file in *dir_output*</li>
</ul>
</details>

***
