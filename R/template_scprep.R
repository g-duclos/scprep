#' Template function for building scRNA-Seq data objects
#'
#' This function builds an ExpressionSet, Seurat, or SingleCellExperiment object with scRNA data
#' @param dir_output Output directory for analysis results
#' @return ExpressionSet, Seurat, or SingleCellExperiment object (depending on output_type parameter)
#' @export
#' @examples
#' # Build ExpressionSet (default)
#' template_scprep(dir_output="/path/to/output")
#' # Output type can be specified in scprep_parameters.csv
#

#
template_scprep <- function(
	dir_output=dir_output) {
	#
	cat(paste(Sys.time()), "\n")
	cat("Initiate Preparation of scRNA-Seq Data", "\n")

	# Read Parameters
	parameters <- read.csv(file.path(dir_output, "scprep_parameters.csv"), stringsAsFactors=FALSE, row.names=1)
	
	# Read annotation
	annotation <- suppressWarnings(read.csv(file.path(dir_output, "scprep_annotation.csv"), stringsAsFactors=FALSE))
	rownames(annotation) <- annotation$Sample_ID

	# Pass parameters to variables & incorporate into ExpressionSet
	param.list <- lapply(rownames(parameters), function(param) {
		#
		cat(paste(parameters[param, "Brief"], parameters[param, "Selection"]), "\n")
		return(parameters[param, "Selection"])
		#
	})
	#
	names(param.list) <- rownames(parameters)
	#
	param.list["mem"] <- list(as.numeric(param.list[["mem"]]))
	param.list["min_umi"] <- list(as.numeric(param.list[["min_umi"]]))
	param.list["min_gene"] <- list(as.numeric(param.list[["min_gene"]]))
	param.list["max_mito"] <- list(as.numeric(param.list[["max_mito"]]))
	param.list["gene_filter"] <- list(as.numeric(param.list[["gene_filter"]]))

	# Parameters
	cat("Initiate ESet Build", "\n")
	#
	cat("Memory-per-sample:", paste(param.list[["mem"]], "G", sep=""), "\n")
	nsamples <- length(annotation$Sample_ID)
	cat("Sample Number:", nsamples, "\n")

	# Initiate data object
	cat("Initiate", param.list[["output_type"]], "object", "\n")
	dataset <- scprep::scprep_eset_build(
		sample_paths=file.path(param.list[["dir_input"]], annotation$Sample_ID),
		annotation=annotation,
		file_type=param.list[["file_type"]],
		gene_id=param.list[["gene_id"]],
		vdj=param.list[["vdj"]],
		cite=param.list[["cite"]],
		atac=param.list[["atac"]],
		output_type=param.list[["output_type"]])

	# Common processing for all object types
	
	# BiomaRt Feature Annotation (now supports all object types)
	cat("BiomaRt Feature Annotation of", param.list[["output_type"]], "object", "\n")
	dataset <- scprep::scprep_eset_biomart(
		dataset=dataset,
		ensembl_target=param.list[["ensembl"]],
		reference=unique(annotation[,"Reference"]))

	# Cell filtering (now supports all object types)
	cat("Barcode Filter:", "\n")
	cat("Cell: High Transcript Content & Low Mitochondrial Content", "\n")
	cat("Dead: High Mitochondrial Content", "\n")
	cat("Debris: Low Transcript Content", "\n")
	cell_filter_assignment <- scprep::scprep_cell_filter_multi(
		dataset=dataset,
		min_umi=param.list[["min_umi"]],
		min_gene=param.list[["min_gene"]],
		max_mito=param.list[["max_mito"]])
	
	# Add cell filter results to metadata based on object type
	if (param.list[["output_type"]] == "eset") {
		dataset$Cell_Filter <- as.factor(cell_filter_assignment)
		cells <- dataset$ID[which(dataset$Cell_Filter == "Cell")]
	} else if (param.list[["output_type"]] == "seurat") {
		dataset@meta.data$Cell_Filter <- as.factor(cell_filter_assignment)
		cells <- rownames(dataset@meta.data)[which(dataset@meta.data$Cell_Filter == "Cell")]
	} else if (param.list[["output_type"]] == "sce") {
		SummarizedExperiment::colData(dataset)$Cell_Filter <- as.factor(cell_filter_assignment)
		cells <- rownames(SummarizedExperiment::colData(dataset))[which(SummarizedExperiment::colData(dataset)$Cell_Filter == "Cell")]
	}
	
	cat("Total Barcode Assignment:", "\n")
	print(table(as.factor(cell_filter_assignment)))

	# Object-specific processing
	if (param.list[["output_type"]] == "eset") {
		#
		cat("Store Parameters & Random Seeds in ExpressionSet assayData Params", "\n")
		Biobase::assayData(dataset)$Params["scprep_Parameters"] <- list(param.list)
		Biobase::assayData(dataset)$Params["Seeds"] <- list(scprep::scprep_eset_seeds(n.seeds=1000))
		
		# Gene filtering for ExpressionSet
		cat("Select Cells", "\n")
		#
		if (is.na(param.list[["gene_filter"]])) {
			#
			cat(paste("Gene Filter = 'NA': >= 3 Transcripts Detected in Default of 0.1% of Cells", sep=""), "\n")
			cat(paste("NOTE: ALL genes will be used as input", sep=""), "\n")
			genes <- rownames(Biobase::exprs(dataset)[,cells])[which(rowSums(Biobase::exprs(dataset)[,cells] >= 3) >= ((0.1/100)*ncol(Biobase::exprs(dataset)[,cells])))]
			cat(paste("Expressed Genes:", length(genes), "/", nrow(Biobase::exprs(dataset))), "\n")
		} else {
			#
			cat(paste("Gene Filter: >= 3 Transcripts Detected in ", param.list[["gene_filter"]], "% of Cells", sep=""), "\n")
			genes <- rownames(Biobase::exprs(dataset)[,cells])[which(rowSums(Biobase::exprs(dataset)[,cells] >= 3) >= ((param.list[["gene_filter"]]/100)*ncol(Biobase::exprs(dataset)[,cells])))]
			cat(paste("Expressed Genes:", length(genes), "/", nrow(Biobase::exprs(dataset))), "\n")
		}
		#
		expressed <- c(rep("Expressed", length(genes)), rep("Not_Expressed", length(setdiff(rownames(Biobase::exprs(dataset)), genes))))
		names(expressed) <- c(genes, setdiff(rownames(Biobase::exprs(dataset)), genes))
		#
		Biobase::fData(dataset)$Gene_Filter <- as.factor(expressed[rownames(Biobase::exprs(dataset))])

		# save ExpressionSet as RDS
		cat("Save ExpressionSet", "\n")
		saveRDS(dataset, file.path(dir_output, "ExpressionSet.rds"));
		
	} else {
		# For Seurat and SingleCellExperiment objects
		cat("Gene filtering is currently only available for ExpressionSet objects", "\n")
		cat("BiomaRt annotation and cell filtering have been applied to your", param.list[["output_type"]], "object", "\n")
		cat("For gene filtering, please use object-specific methods or convert to ExpressionSet", "\n")
		
		# Save object with appropriate filename
		output_filename <- paste0(param.list[["output_type"]], ".rds")
		cat("Save", param.list[["output_type"]], "\n")
		saveRDS(dataset, file.path(dir_output, output_filename));
	}
	
	#
	return(dataset)
	#
}
