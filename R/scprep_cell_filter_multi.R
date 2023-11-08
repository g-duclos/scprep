#' Function for multi-sample discrimination of "Cell" barcodes from "Dead" cells and "Debris" based on hard cut-offs for UMI transcript & gene content, and percent mitochondrial UMIs
#'
#' Across all samples (Sample designation in ExpressionSet pData slot), this function designates high transcript/gene content & low percent mitochondrial content barcodes as "Cell", high percent mitochondrial content barcodes as "Dead", & low transcript/gene content barcodes as "Debris"
#' @param dataset ExpressionSet object with scRNA data
#' @param min_umi UMI transcript threshold
#' @param min_gene Genes detected threshold
#' @param max_mito Percent mitochondrial transcript threshold (0 = No mito transcripts allowed, 1 = No limit for mito transcripts)
#' @return Multi-sample Cell, Dead, & Debris assignment
#' @export
#' @examples
#' scprep_cell_filter_multi(dataset=ExpressionSet, min_umi=100, min_gene=0, max_mito=0.2)
#

#
scprep_cell_filter_multi <- function(
	dataset=dataset,
	min_umi=100,
	min_gene=0,
	max_mito=0.2) {
	#
	cell.filter.out <- unlist(lapply(unique(dataset$Sample), function(sample) {
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
	#
	return(cell.filter.out)
	#
}