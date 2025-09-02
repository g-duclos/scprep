#' Function to annotate feature data in ExpressionSet, Seurat, or SingleCellExperiment objects
#'
#' This function annotates feature data with scRNA data using biomaRt
#' @param dataset ExpressionSet, Seurat, or SingleCellExperiment object
#' @param ensembl_target Species input for biomaRt feature annotation
#' @param reference Reference genome (on annotation.csv & SampleSheet.csv): "refdata-cellranger-GRCh38-3.0.0" or "GRCh38_pre-mRNA" for Ensembl 93; "refdata-gex-GRCh38-2020-A" or "GRCh38_pre-mRNA-2020-A" for Ensembl 98
#' @return Object of same type as input containing annotated feature data
#' @export
#' @examples
#' scprep_eset_biomart(dataset=object, ensembl_target="hsapiens_gene_ensembl", reference="refdata-cellranger-GRCh38-3.0.0")
#

#
scprep_eset_biomart <- function(
	dataset=dataset,
	ensembl_target=ensembl_target,
	reference=reference) {
	#
	cat("BiomaRt Feature Annotation", "\n")
	
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
	
	cat("Processing", object_type, "object", "\n")

	# Extract expression matrix and row names based on object type
	if (object_type == "ExpressionSet") {
		expr_matrix <- exprs(dataset)
		gene_names <- rownames(dataset)
		cell_names <- colnames(dataset)
		cell_data <- pData(dataset)
	} else if (object_type == "Seurat") {
		expr_matrix <- Seurat::GetAssayData(dataset, slot = "counts")
		gene_names <- rownames(dataset)
		cell_names <- colnames(dataset)
		cell_data <- dataset@meta.data
	} else if (object_type == "SingleCellExperiment") {
		expr_matrix <- SingleCellExperiment::counts(dataset)
		gene_names <- rownames(dataset)
		cell_names <- colnames(dataset)
		cell_data <- as.data.frame(SummarizedExperiment::colData(dataset))
	}

	# Read stored feature annotation
	if (reference == "refdata-cellranger-GRCh38-3.0.0" | reference == "GRCh38_pre-mRNA") {
		#
		cat("Retrieve Annotation for hsapiens Ensembl 93", "\n")
		data(biomart_annotation_ens_93, package="scprep")
		biomart_annotation <- biomart_annotation_ens_93
		#
	} else if (reference == "refdata-gex-GRCh38-2020-A" | reference == "GRCh38_pre-mRNA-2020-A") {
		#
		cat("Retrieve Annotation for hsapiens Ensembl 98", "\n")			
		data(biomart_annotation_ens_98, package="scprep")
		biomart_annotation <- biomart_annotation_ens_98
		#
	}
	
	# Select only genes from annotation present in the dataset
	biomart_annotation <- biomart_annotation[intersect(gene_names, rownames(biomart_annotation)), ]
	
	# Add feature annotation to object based on type
	if (object_type == "ExpressionSet") {
		# Initialize featureData slot of ExpressionSet
		fData(dataset) <- data.frame(row.names=gene_names)
		for (fdata in colnames(biomart_annotation)) {
			fData(dataset)[[fdata]] <- biomart_annotation[intersect(rownames(biomart_annotation), gene_names), fdata]
		}
	} else if (object_type == "Seurat") {
		# Add feature metadata to Seurat object
		feature_df <- data.frame(row.names = gene_names)
		for (fdata in colnames(biomart_annotation)) {
			feature_df[[fdata]] <- biomart_annotation[intersect(rownames(biomart_annotation), gene_names), fdata]
		}
		dataset[["RNA"]]@meta.features <- feature_df
	} else if (object_type == "SingleCellExperiment") {
		# Add feature metadata to SCE object
		feature_df <- data.frame(row.names = gene_names)
		for (fdata in colnames(biomart_annotation)) {
			feature_df[[fdata]] <- biomart_annotation[intersect(rownames(biomart_annotation), gene_names), fdata]
		}
		SummarizedExperiment::rowData(dataset) <- feature_df
	}

	# Quantify Biotype information per cell
	cat("Quantify Biotype Annotation per Cell", "\n")
	#
	biotype.vars <- c("Biotype_Edit", "PC_Biotype")
	#
	for (biotype.var in biotype.vars) {
		# Get biotype info based on object type
		if (object_type == "ExpressionSet") {
			biotypes <- unique(fData(dataset)[[biotype.var]])
			feature_data <- fData(dataset)
		} else if (object_type == "Seurat") {
			if (biotype.var %in% colnames(dataset[["RNA"]]@meta.features)) {
				biotypes <- unique(dataset[["RNA"]]@meta.features[[biotype.var]])
				feature_data <- dataset[["RNA"]]@meta.features
			} else {
				next
			}
		} else if (object_type == "SingleCellExperiment") {
			if (biotype.var %in% colnames(SummarizedExperiment::rowData(dataset))) {
				biotypes <- unique(SummarizedExperiment::rowData(dataset)[[biotype.var]])
				feature_data <- as.data.frame(SummarizedExperiment::rowData(dataset))
			} else {
				next
			}
		}
		
		biotypes <- biotypes[!is.na(biotypes)]
		#
		biotype.quant <- sapply(biotypes, function (biotype) {
			cat(paste("Biotype:", biotype), "\n")
			gene_idx <- which(feature_data[[biotype.var]] == biotype)
			if (length(gene_idx) > 0) {
				sapply(1:ncol(expr_matrix), function(cell) {
					return(sum(expr_matrix[gene_idx, cell]))
				})
			} else {
				rep(0, ncol(expr_matrix))
			}
		})
		colnames(biotype.quant) <- biotypes
		rownames(biotype.quant) <- cell_names

		# Add biotype quantification to metadata
		for (biotype in biotypes) {
			if (object_type == "ExpressionSet") {
				if ("ID" %in% colnames(pData(dataset))) {
					dataset[[biotype]] <- biotype.quant[dataset$ID, biotype]
				} else {
					pData(dataset)[[biotype]] <- biotype.quant[, biotype]
				}
			} else if (object_type == "Seurat") {
				dataset@meta.data[[biotype]] <- biotype.quant[, biotype]
			} else if (object_type == "SingleCellExperiment") {
				SummarizedExperiment::colData(dataset)[[biotype]] <- biotype.quant[, biotype]
			}
		}
	}

	# Quantify Chromosomal Annotation per Cell
	cat("Quantify Chromosomal Annotation per Cell", "\n")
	#
	chroms <- c(1:22, "X", "Y")
	
	# Get chromosome info based on object type  
	if (object_type == "ExpressionSet") {
		if ("Chr" %in% colnames(fData(dataset))) {
			feature_chr <- fData(dataset)$Chr
		} else {
			cat("Warning: Chr column not found in feature data, skipping chromosome quantification", "\n")
			return(dataset)
		}
	} else if (object_type == "Seurat") {
		if ("Chr" %in% colnames(dataset[["RNA"]]@meta.features)) {
			feature_chr <- dataset[["RNA"]]@meta.features$Chr
		} else {
			cat("Warning: Chr column not found in feature data, skipping chromosome quantification", "\n")
			return(dataset)
		}
	} else if (object_type == "SingleCellExperiment") {
		if ("Chr" %in% colnames(SummarizedExperiment::rowData(dataset))) {
			feature_chr <- SummarizedExperiment::rowData(dataset)$Chr
		} else {
			cat("Warning: Chr column not found in feature data, skipping chromosome quantification", "\n")
			return(dataset)
		}
	}
	
	#
	chr.quant <- sapply(chroms, function (chr) {
		cat(paste("Chromosome:", chr), "\n")
		gene_idx <- which(feature_chr == chr)
		if (length(gene_idx) > 0) {
			sapply(1:ncol(expr_matrix), function(cell) {
				return(sum(expr_matrix[gene_idx, cell]))
			})
		} else {
			rep(0, ncol(expr_matrix))
		}
	})
	colnames(chr.quant) <- chroms
	rownames(chr.quant) <- cell_names

	# Add chromosome quantification to metadata
	for (chr in chroms) {
		chr_col_name <- paste("Chr_", chr, sep="")
		if (object_type == "ExpressionSet") {
			if ("ID" %in% colnames(pData(dataset))) {
				dataset[[chr_col_name]] <- chr.quant[dataset$ID, chr]
			} else {
				pData(dataset)[[chr_col_name]] <- chr.quant[, chr]
			}
		} else if (object_type == "Seurat") {
			dataset@meta.data[[chr_col_name]] <- chr.quant[, chr]
		} else if (object_type == "SingleCellExperiment") {
			SummarizedExperiment::colData(dataset)[[chr_col_name]] <- chr.quant[, chr]
		}
	}
	#
	return(dataset)
	#
}