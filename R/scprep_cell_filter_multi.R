#' Function for multi-sample discrimination of "Cell" barcodes from "Dead" cells and "Debris" based on hard cut-offs for UMI transcript & gene content, and percent mitochondrial UMIs
#'
#' Across all samples (Sample designation in object metadata slot), this function designates high transcript/gene content & low percent mitochondrial content barcodes as "Cell", high percent mitochondrial content barcodes as "Dead", & low transcript/gene content barcodes as "Debris"
#' @param dataset ExpressionSet, Seurat, or SingleCellExperiment object with scRNA data
#' @param min_umi UMI transcript threshold
#' @param min_gene Genes detected threshold
#' @param max_mito Percent mitochondrial transcript threshold (0 = No mito transcripts allowed, 1 = No limit for mito transcripts)
#' @return Multi-sample Cell, Dead, & Debris assignment
#' @export
#' @examples
#' scprep_cell_filter_multi(dataset=object, min_umi=100, min_gene=0, max_mito=0.2)
#
#
scprep_cell_filter_multi <- function(
	dataset=dataset,
	min_umi=100,
	min_gene=0,
	max_mito=0.2) {
	#
	# Detect object type and extract metadata
	object_type <- class(dataset)[1]
	if (object_type == "Seurat") {
		object_type <- "Seurat"
		cell_data <- dataset@meta.data
	} else if ("SingleCellExperiment" %in% class(dataset)) {
		object_type <- "SingleCellExperiment" 
		cell_data <- as.data.frame(SummarizedExperiment::colData(dataset))
	} else if ("ExpressionSet" %in% class(dataset)) {
		object_type <- "ExpressionSet"
		cell_data <- pData(dataset)
	} else {
		stop("Unsupported object type. Please provide ExpressionSet, Seurat, or SingleCellExperiment object.")
	}
	
	# Check if Sample column exists
	if (!"Sample" %in% colnames(cell_data)) {
		stop("Sample column not found in object metadata")
	}
	
	# Get unique samples
	unique_samples <- unique(cell_data$Sample)
	cat("Processing", length(unique_samples), "samples:", paste(unique_samples, collapse = ", "), "\n")
	
	#
	cell.filter.out <- unlist(lapply(unique_samples, function(sample) {
		#
		cell.filter.sub <- scprep_cell_filter(dataset=dataset,
			sample=sample,
			min_umi=min_umi,
			min_gene=min_gene,
			max_mito=max_mito)
		#
		cat(paste("Barcode Assignment:", sample), "\n")
		print(table(as.factor(cell.filter.sub)))
		#
		return(cell.filter.sub)
		#
	}))

	# Add results to object metadata
	if (object_type == "Seurat") {
		#
		dataset@meta.data$Cell_Filter <- as.factor(cell.filter.out)
		#
	} else if ("SingleCellExperiment" %in% class(dataset)) {
		#
		SummarizedExperiment::colData(dataset)$Cell_Filter <- as.factor(cell_filter_assignment)
		#
	} else if ("ExpressionSet" %in% class(dataset)) {
		#
		dataset$Cell_Filter <- as.factor(cell.filter.out)
		#
	}
	#
	return(dataset)
	#
}