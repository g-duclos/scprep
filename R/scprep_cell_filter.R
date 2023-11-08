#' Function for discriminating "Cell" barcodes from "Dead" cells and "Debris" based on hard cut-offs for UMI transcript & gene content, and percent mitochondrial UMIs
#'
#' This function designates high transcript/gene content & low percent mitochondrial content barcodes as "Cell", high percent mitochondrial content barcodes as "Dead", & low transcript/gene content barcodes as "Debris"
#' @param dataset ExpressionSet object with scRNA data
#' @param sample Sample ID stored under "Sample" in pData slot
#' @param min_umi UMI transcript threshold
#' @param min_gene Genes detected threshold
#' @param max_mito Percent mitochondrial transcript threshold (0 = No mito transcripts allowed, 1 = No limit for mito transcripts)
#' @return Sample-specific Cell, Dead, & Debris assignment
#' @export
#' @examples
#' cell_filter(dataset=ExpressionSet, sample=sample, min_umi=100, min_gene=0, max_mito=0.2)
#

#
scprep_cell_filter <- function(
	dataset=dataset,
	sample=sample,
	min_umi=100,
	min_gene=0,
	max_mito=0.2) {

	# Calculate % mito transcripts
	percent.mito <- 100*(dataset$Mitochondrial[which(dataset$Sample == sample)]/dataset$UMIs[which(dataset$Sample == sample)])
	names(percent.mito) <- dataset$ID[which(dataset$Sample == sample)]

	# Dead = Above max %mito threshold 
	dead <- names(percent.mito)[which(percent.mito > 100*max_mito)]

	# Debris = Below min transcript threshold
	debris <- setdiff(dataset$ID[which(dataset$Sample == sample & dataset$UMIs < min_umi | dataset$Sample == sample & dataset$Genes < min_gene)], dead)

	# Cell = Below max %mito and Above min thranscript threshold
	cell <- setdiff(dataset$ID[which(dataset$Sample == sample)], c(dead, debris))

	#
	assignment <- c(rep("Cell", length(cell)), rep("Dead", length(dead)), rep("Debris", length(debris)))
	names(assignment) <- c(cell, dead, debris)

	# Reorder to match cell order in ESet
	assignment <- assignment[dataset$ID[which(dataset$Sample == sample)]]
	#
	return(assignment)
	#
}
