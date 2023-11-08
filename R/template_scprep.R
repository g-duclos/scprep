#' Template function for building ExpressionSet object
#'
#' This function builds an ExpressionSet object with scRNA data
#' @param dir_output Output directory for analysis results
#' @return ExpressionSet object
#' @export
#' @examples
#' template_scprep(dir_output="/new/directory")
#

#
template_scprep <- function(
	dir_output=dir_output) {

	# Load start up 'scprep' message
	#scprep_load_startup_message(load=TRUE);
	#
	cat(paste(Sys.time()), "\n")
	cat("Build ExpressionSet", "\n")
	
	# Read Parameters
	parameters <- read.csv(file.path(dir_output, "scprep_parameters.csv"), stringsAsFactors=FALSE, row.names=1)
	
	# Read annotation
	annotation <- read.csv(file.path(dir_output, "scprep_annotation.csv"), stringsAsFactors=FALSE)
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
	#

	# Parameters
	cat("Initiate ESet Build", "\n")
	#
	cat("Memory-per-sample:", paste(param.list[["mem"]], "G", sep=""), "\n")
	nsamples <- length(annotation$Sample_ID)
	cat("Sample Number:", nsamples, "\n")

	# Initiate ExpressionSet
	cat("Initiate ExpressionSet", "\n")
	dataset <- scprep_eset_build(
		sample_paths=file.path(annotation$Sample_Path, annotation$Sample_ID),
		annotation=annotation,
		vdj=param.list[["vdj"]],
		cite=param.list[["cite"]],
		atac=param.list[["atac"]])

	# Biomart
	cat("BiomaRt Feature Annotation of ExpressionSet", "\n")
	dataset <- scprep_eset_biomart(
		dataset=dataset,
		ensembl_target=param.list[["ensembl"]],
		reference=unique(annotation[,"Reference"]))
	#
	cat("Store Parameters & Random Seeds in ExpressionSet assayData Params", "\n")
	assayData(dataset)$Params["Parameters"] <- list(param.list)
	assayData(dataset)$Params["Seeds"] <- list(eset_seeds(n.seeds=1000))
	#
	cat("Barcode Filter:", "\n")
	cat("Cell: High Transcript Content & Low Mitochondrial Content", "\n")
	cat("Dead: High Mitochondrial Content", "\n")
	cat("Debris: Low Transcript Content", "\n")
	dataset$Cell_Filter <- as.factor(scprep_cell_filter_multi(
		dataset=dataset,
		min_umi=param.list[["min_umi"]],
		min_gene=param.list[["min_gene"]],
		max_mito=param.list[["max_mito"]]))
	#
	cat("Total Barcode Assignment:", "\n")
	print(table(dataset$Cell_Filter))
	
	#
	cat("Select Cells", "\n")
	cells <- dataset$ID[which(dataset$Cell_Filter == "Cell")]
	#
	if (is.na(param.list[["gene_filter"]])) {
		#
		cat(paste("Gene Filter = 'NA': >= 3 Transcripts Detected in Default of 0.1% of Cells", sep=""), "\n")
		cat(paste("NOTE: ALL genes will be used as Seurat input", sep=""), "\n")
		genes <- rownames(exprs(dataset)[,cells])[which(rowSums(exprs(dataset)[,cells] >= 3) >= ((0.1/100)*ncol(exprs(dataset)[,cells])))]
		cat(paste("Expressed Genes:", length(genes), "/", nrow(exprs(dataset))), "\n")
	} else {
		#
		cat(paste("Gene Filter: >= 3 Transcripts Detected in ", param.list[["gene_filter"]], "% of Cells", sep=""), "\n")
		genes <- rownames(exprs(dataset)[,cells])[which(rowSums(exprs(dataset)[,cells] >= 3) >= ((param.list[["gene_filter"]]/100)*ncol(exprs(dataset)[,cells])))]
		cat(paste("Expressed Genes:", length(genes), "/", nrow(exprs(dataset))), "\n")
	}
	#
	expressed = c(rep("Expressed", length(genes)), rep("Not_Expressed", length(setdiff(rownames(exprs(dataset)), genes))))
	names(expressed) <- c(genes, setdiff(rownames(exprs(dataset)), genes))
	#
	fData(dataset)$Gene_Filter <- as.factor(expressed[rownames(exprs(dataset))])

	# save ExpressionSet as RDS
	cat("Save ExpressionSet", "\n")
	saveRDS(dataset, file.path(dir_output, "ExpressionSet.rds"));
	
	# Calculate sample number from annotation
	nsamples <- length(unique(dataset$Sample))
	cat(paste("Sample number:", nsamples), "\n")

	#
	return(dataset)
	#
}
