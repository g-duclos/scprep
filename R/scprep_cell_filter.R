#' Function for discriminating "Cell" barcodes from "Dead" cells and "Debris" based on hard cut-offs for UMI transcript & gene content, and percent mitochondrial UMIs
#'
#' This function designates high transcript/gene content & low percent mitochondrial content barcodes as "Cell", high percent mitochondrial content barcodes as "Dead", & low transcript/gene content barcodes as "Debris"
#' @param dataset ExpressionSet, Seurat, or SingleCellExperiment object with scRNA data
#' @param sample Sample ID stored under "Sample" in metadata slot
#' @param min_umi UMI transcript threshold
#' @param min_gene Genes detected threshold
#' @param max_mito Percent mitochondrial transcript threshold (0 = No mito transcripts allowed, 1 = No limit for mito transcripts)
#' @return Sample-specific Cell, Dead, & Debris assignment
#' @export
#' @examples
#' scprep_cell_filter(dataset=object, sample=sample, min_umi=100, min_gene=0, max_mito=0.2)

#
scprep_cell_filter <- function(
	dataset=dataset,
	sample=sample,
	min_umi=100,
	min_gene=0,
	max_mito=0.2) {

	# Detect object type
	object_type <- class(dataset)[1]
	if (object_type == "Seurat") {
		object_type <- "Seurat"
	} else if ("SingleCellExperiment" %in% class(dataset)) {
		object_type <- "SingleCellExperiment" 
	} else if ("ExpressionSet" %in% class(dataset)) {
		object_type <- "ExpressionSet"
	} else {
		stop("Unsupported object type. Please provide ExpressionSet, Seurat, or SingleCellExperiment object.")
	}

	# Extract metadata based on object type
	if (object_type == "ExpressionSet") {
		cell_data <- pData(dataset)
		cell_ids <- cell_data$ID
	} else if (object_type == "Seurat") {
		cell_data <- dataset@meta.data
		cell_ids <- rownames(cell_data)
	} else if (object_type == "SingleCellExperiment") {
		cell_data <- as.data.frame(SummarizedExperiment::colData(dataset))
		cell_ids <- rownames(cell_data)
	}

	# Check required columns exist
	required_cols <- c("Sample", "Mitochondrial", "UMIs", "Genes")
	missing_cols <- setdiff(required_cols, colnames(cell_data))
	if (length(missing_cols) > 0) {
		stop(paste("Missing required columns in metadata:", paste(missing_cols, collapse = ", ")))
	}

	# Filter for specified sample
	sample_mask <- cell_data$Sample == sample
	if (sum(sample_mask) == 0) {
		stop(paste("Sample", sample, "not found in dataset"))
	}

	sample_cells <- cell_ids[sample_mask]
	sample_mito <- cell_data$Mitochondrial[sample_mask]
	sample_umis <- cell_data$UMIs[sample_mask]
	sample_genes <- cell_data$Genes[sample_mask]

	# Calculate % mito transcripts
	percent.mito <- 100 * (sample_mito / sample_umis)
	names(percent.mito) <- sample_cells

	# Dead = Above max %mito threshold 
	dead <- names(percent.mito)[which(percent.mito > 100*max_mito)]

	# Debris = Below min transcript threshold
	debris_mask <- (sample_umis < min_umi) | (sample_genes < min_gene)
	debris <- setdiff(sample_cells[debris_mask], dead)

	# Cell = Below max %mito and Above min transcript threshold
	cell <- setdiff(sample_cells, c(dead, debris))

	# Create assignment vector
	assignment <- c(rep("Cell", length(cell)), rep("Dead", length(dead)), rep("Debris", length(debris)))
	names(assignment) <- c(cell, dead, debris)

	# Reorder to match cell order in dataset for the specific sample
	assignment <- assignment[sample_cells]
	
	return(assignment)
}
