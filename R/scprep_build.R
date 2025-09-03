#' Function to initiate single-cell data object
#'
#' This function builds an ExpressionSet, Seurat, or SingleCellExperiment object with scRNA data
#' @param sample_paths Character vector of paths to filtered_feature_bc_matrix.h5 for all samples in project
#' @param input_data Logical; Default is FALSE to ignore; if TRUE, input_data should be a list of pre-loaded count matrices (instead of file paths) - IMPORTANT: if you use this form of input, it will override 'sample_paths' argument! If the 'Sample_ID' column in annotation.csv does not match the names of the elements of this list, the 'annotation' input will be ignored and your metadata will only contain the names of the elements in this list (if the list elements are NOT named, they will be give numeric IDs corresponding to order of the list)
#' @param annotation Matrix of project-specific annotation.csv
#' @param file_type File type ("h5" or "mtx"); if "mtx", the following files must be in sample directory: matrix.mtx, barcodes.tsv, genes.tsv or features.tsv
#' @param gene_id Gene ID type ("ensembl" or "symbol")
#' @param vdj Data is multimodal RNA & VDJ
#' @param cite Data is multimodal RNA & Protein
#' @param cite_ignore FALSE to ignore; Default is TRUE, ignore the "Antibody Capture" slot in the h5 file input and only extract/store the "Gene Expression" slot; This is only relevant if your data contains both RNA and protein modalities and you want to ignore the protein modality
#' @param atac Data is multimodal RNA & ATAC
#' @param atac_ignore FALSE to ignore; Default is TRUE, ignore the "ATAC" slot in the h5 file input and only extract/store the "Gene Expression" slot; This is only relevant if your data contains both RNA and ATAC modalities and you want to ignore the ATAC modality
#' @param output_type Output object type ("eset" (ExpressionSet), "seurat", or "sce" (SingleCellExperiment))
#' @param verbose Logical; if TRUE (default), print progress messages; if FALSE, suppress output
#' @return ExpressionSet, Seurat, or SingleCellExperiment object containing expression data for all samples & annotation
#' @export
#' @examples
#' \dontrun{
#' # Build ExpressionSet (default)
#' dataset <- scprep_build(sample_paths=sample_paths, input_data=FALSE, annotation=annotation, 
#'                              file_type="h5", gene_id="symbol", vdj=FALSE, 
#'                              cite=FALSE, atac=FALSE)
#' }
#
scprep_build <- function(
	sample_paths=NULL,
	input_data=FALSE,
	annotation=NULL,
	file_type="h5",
	gene_id="symbol",
	vdj=FALSE,
	cite=FALSE,
	cite_ignore=TRUE,
	atac=FALSE,
	atac_ignore=TRUE,
	output_type="eset",
	verbose=TRUE) {
	#
	if (verbose) cat("Building", output_type, "object", "\n")
	
	# Handle input_data parameter when provided as list of matrices
	if (!isFALSE(input_data) && is.list(input_data)) {
		if (verbose) cat("Using pre-loaded count matrices from input_data", "\n")
		
		# Handle case where input_data list elements are named
		if (is.null(names(input_data))) {
			names(input_data) <- paste("Sample", 1:length(input_data), sep="_")
			if (verbose) cat("Input data not named, assigning numeric IDs:", paste(names(input_data), collapse=", "), "\n")
		}
		
		# Process input_data matrices
		for (i in 1:length(input_data)) {
			# Convert to matrix if needed
			counts <- as.matrix(input_data[[i]])
			
			# Add sample suffix to cell barcodes
			colnames(counts) <- paste(colnames(counts), "-", i, sep="")
			
			if (i == 1) {
				counts.all <- counts
				if (verbose) cat("Loaded matrix for", names(input_data)[i], ":", ncol(counts), "cells", "\n")
			} else {
				counts.all <- cbind(counts.all, counts)
				if (verbose) cat("Loaded matrix for", names(input_data)[i], ":", ncol(counts), "cells", "\n")
			}
		}
		
		# Create annotation based on input_data names
		sample_names <- names(input_data)
		n.cells.per.sample <- sapply(input_data, ncol)
		
		# Check if annotation Sample_ID matches input_data names
		if (!is.null(annotation) && all(sample_names %in% annotation$Sample_ID)) {
			if (verbose) cat("Using provided annotation matching input_data names", "\n")
			# Reorder annotation to match input_data order
			annotation <- annotation[match(sample_names, annotation$Sample_ID), ]
			
			annot.all <- do.call(rbind, lapply(1:length(sample_names), function(i) {
				return(sapply(1:ncol(annotation), function(x) {
					return(rep(annotation[i, x], n.cells.per.sample[i]))
				}))
			}))
			colnames(annot.all) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
		} else {
			if (verbose) cat("Annotation Sample_ID does not match input_data names or annotation not provided", "\n")
			if (verbose) cat("Creating minimal annotation from input_data names", "\n")
			# Create minimal annotation from input_data names
			annot.all <- do.call(rbind, lapply(1:length(sample_names), function(i) {
				return(data.frame(Sample_ID = rep(sample_names[i], n.cells.per.sample[i]),
								stringsAsFactors = FALSE))
			}))
			colnames(annot.all) <- "Sample"
		}
		rownames(annot.all) <- colnames(counts.all)
		
		if (verbose) cat("Total cells loaded:", ncol(counts.all), "\n")
		
		# Skip file reading logic and jump to metadata creation
		skip_file_reading <- TRUE
	} else {
		skip_file_reading <- FALSE
	}
	
	# Set use.names parameter based on gene_id preference
	use_names <- ifelse(gene_id == "symbol", TRUE, FALSE)
	#
	if (!skip_file_reading) {
		if (vdj == FALSE) {
		#
		counts.file.name <- "filtered_feature_bc_matrix.h5"
		#
	} else if (vdj == TRUE) {
		#
		counts.file.name <- "sample_filtered_feature_bc_matrix.h5"
		#
		vdj.file.name <- "filtered_contig_annotations.csv"
		# Initiate list of vdj outputs that will be loaded into Eset
		vdj.list <- list()
		#
	}
	#
	if (length(sample_paths) == 1) {
		#
		if (file_type == "mtx") {
			# If file type is mtx (vs h5), set dir_counts to file dir without filename
			dir_counts <- sample_paths
			# Set vdj/cite/atac to FALSE
			vdj <- FALSE
			atac <- FALSE
			cite <- FALSE
		} else {
			# Default file_type is "h5", which requires dir_counts to include dir + filename
			dir_counts <- file.path(sample_paths, counts.file.name)
		}
		#
		if (cite == FALSE & atac == FALSE) {
			#
			if (file_type == "mtx") {
				#
				counts.all <- suppressMessages(Seurat::Read10X(dir_counts, gene.column=1))
				#
			} else {
				#
				counts.all <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=use_names))
				#
			}
			#
			if (typeof(counts.all) == "S4") {
				#
				counts.all <- as.matrix(counts.all)
				#
			} else if (typeof(counts.all) == "list") {
				#
				counts.all <- as.matrix(counts.all[["Gene Expression"]])
				#
			}
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
		} else if (cite == TRUE && cite_ignore == FALSE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
			if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			#
			counts.cite.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
			if (verbose) cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			colnames(counts.cite.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
		} else if (atac == TRUE && atac_ignore == FALSE && cite == FALSE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
			if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
		} else {
			# When cite_ignore=TRUE or atac_ignore=TRUE, only extract Gene Expression
			if (file_type == "mtx") {
				#
				counts.all <- suppressMessages(Seurat::Read10X(dir_counts, gene.column=1))
				#
			} else {
				#
				counts.all <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=use_names))
				#
			}
			#
			if (typeof(counts.all) == "S4") {
				#
				counts.all <- as.matrix(counts.all)
				#
			} else if (typeof(counts.all) == "list") {
				#
				counts.all <- as.matrix(counts.all[["Gene Expression"]])
				#
			}
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			if (verbose) cat("Read Filtered GEX Matrix (ignoring multimodal data):", basename(sample_paths), "\n")
			if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
		}
		# Check if number of unique samples equals sample number specified in annotation
		n.samples <- as.numeric(unique(mat.barcodes[seq(2,length(mat.barcodes),2)]))
		#
		if (length(n.samples) == nrow(annotation)) {
			#
			n.barcodes <- unlist(lapply(unique(mat.barcodes[seq(2,length(mat.barcodes),2)]), function(i) {
				return(length(mat.barcodes[seq(2,length(mat.barcodes),2)][which(mat.barcodes[seq(2,length(mat.barcodes),2)] == i)]))
				}))
			#
			annot.all <- do.call(rbind, lapply(n.samples, function(i) {
				if (verbose) cat(n.barcodes[i], "Cells Detected in", annotation$Sample_ID[i], "\n")
				return(sapply(1:ncol(annotation), function(x) {
					return(rep(annotation[i, x], n.barcodes[i]))
				}))
			}))
			colnames(annot.all) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
			rownames(annot.all) <- colnames(counts.all)
			#
			} else {
				if (verbose) cat("Matrix Sample Number Does NOT Match Sample Number Specified in Annotation!", "\n")
			}
		#			
		} else {
			#
			for (i in 1:length(sample_paths)) {
				#
				if (file_type == "mtx") {
					# If file type is mtx (vs h5), set dir_counts to file dir without filename
					dir_counts <- sample_paths[i]
					# Set vdj/cite/atac to FALSE
					vdj <- FALSE
					atac <- FALSE
					cite <- FALSE
				} else {
					# Default file_type is "h5", which requires dir_counts to include dir + filename
					dir_counts <- file.path(sample_paths[i], counts.file.name)
				}
				#
				if (i == 1) {
					#
					if (cite == FALSE & atac == FALSE) {		
						#
						if (file_type == "mtx") {
							#
							counts.all <- suppressMessages(Seurat::Read10X(dir_counts, gene.column=1))
							#
						} else {
							#
							counts.all <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=use_names))
							#
						}
						#
						if (typeof(counts.all) == "S4") {
							#
							counts.all <- as.matrix(counts.all)
							#
						} else if (typeof(counts.all) == "list") {
							#
							counts.all <- as.matrix(counts.all[["Gene Expression"]])
							#
						}
						#
						mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
						colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						#
						if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
						#
					} else if (cite == TRUE) {
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
						if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						#
						counts.cite.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
						if (verbose) cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
						#
						mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
						colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						colnames(counts.cite.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						#
						if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
						#
					} else if (atac == TRUE & cite == FALSE) {
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
						if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						#
						mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
						colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						#
						if (verbose) cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
						#
					}
					#
					annot.all <- sapply(1:ncol(annotation), function(x) {
						return(rep(annotation[i, x], ncol(counts.all)))
					})
					colnames(annot.all) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
					rownames(annot.all) <- colnames(counts.all)
					#
					} else {
						#
						if (cite == FALSE & atac == FALSE) {
							#
							if (file_type == "mtx") {
								#
								counts <- suppressMessages(Seurat::Read10X(dir_counts, gene.column=1))
								#
							} else {
								#
								counts <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=use_names))
								#
							}
							#
							if (typeof(counts) == "S4") {
								#
								counts <- as.matrix(counts)
								#
							} else if (typeof(counts) == "list") {
								#
								counts <- as.matrix(counts[["Gene Expression"]])
								#
							}
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							if (verbose) cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						} else if (cite == TRUE) {
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
							if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							#
							counts.cite <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
							if (verbose) cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							colnames(counts.cite) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							if (verbose) cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						} else if (atac == TRUE & cite == FALSE) {
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
							if (verbose) cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							if (verbose) cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						}
						#
						annot <- sapply(1:ncol(annotation), function(x) {
							return(rep(annotation[i, x], ncol(counts)))
						})
						colnames(annot) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
						rownames(annot) <- colnames(counts)
						#
						counts.all <- cbind(counts.all, counts)
						#
						if (cite == TRUE && cite_ignore == FALSE) {
							counts.cite.all <- cbind(counts.cite.all, counts.cite)
						}
						#
						annot.all <- rbind(annot.all, annot)
						#
					}
				}
			}
	#
	if (vdj == TRUE) {
		#
		if (length(sample_paths) == 1) {
			#
			if (file.exists(file.path(sample_paths, vdj.file.name)) == TRUE) {
				if (verbose) cat("Read VDJ data:", vdj.file.name, "\n")
				vdj.list[as.character(basename(sample_paths))] <- list(read.csv(file.path(sample_paths, vdj.file.name)))
			} else {
				if (verbose) cat("VDJ data NOT found !!!", "\n")
			}
			#
		} else {
			#
			for (i in 1:length(sample_paths)) {
			#
				if (file.exists(file.path(sample_paths, vdj.file.name)) == TRUE) {
					if (verbose) cat("Read VDJ data:", vdj.file.name, "\n")
					#
					vdj.list[as.character(basename(sample_paths[i]))] <- list(read.csv(file.path(sample_paths[i], vdj.file.name)))
					#
				} else {
					if (verbose) cat("VDJ data NOT found !!!", "\n")
				}
			}
		}
	}
	}
	
	# Create common metadata data frame
	metadata_df <- data.frame(
		ID = colnames(counts.all),
		stringsAsFactors = FALSE
	)
	rownames(metadata_df) <- colnames(counts.all)
	
	# Add annotation variables to metadata
	for (i in 1:ncol(annot.all)) {
		metadata_df[[colnames(annot.all)[i]]] <- as.factor(annot.all[,i])
	}
	
	# Calculate common QC metrics
	if (verbose) cat("Calculate UMIs per cell", "\n")
	metadata_df$UMIs <- colSums(counts.all)

	if (verbose) cat("Calculate Genes per cell", "\n")
	metadata_df$Genes <- colSums(counts.all > 0)

	# Create object based on output type
	if (output_type == "seurat") {
		if (verbose) cat("Creating Seurat object directly", "\n")
		# Create Seurat object directly
		dataset <- Seurat::CreateSeuratObject(
			counts = counts.all,
			meta.data = metadata_df
		)
		# Add protein data if available
		if (cite == TRUE && cite_ignore == FALSE && exists("counts.cite.all")) {
			dataset[["Protein"]] <- Seurat::CreateAssayObject(counts = counts.cite.all)
		}
		# Add VDJ data if available (store in misc slot)
		if (vdj == TRUE && exists("vdj.list")) {
			dataset@misc$VDJ <- vdj.list
		}
		#
		if (verbose) cat("Seurat object created with", ncol(dataset), "cells and", nrow(dataset), "genes", "\n")
		return(dataset)
		#
	} else if (output_type == "sce") {
		if (verbose) cat("Creating SingleCellExperiment object directly", "\n")
		# Create SingleCellExperiment object directly
		dataset <- SingleCellExperiment::SingleCellExperiment(
			assays = list(counts = counts.all),
			colData = metadata_df
		)
		# Add protein data if available
		if (cite == TRUE && cite_ignore == FALSE && exists("counts.cite.all")) {
			SingleCellExperiment::altExp(dataset, "Protein") <- SingleCellExperiment::SingleCellExperiment(
				assays = list(counts = counts.cite.all)
			)
		}
		# Add VDJ data if available (store in metadata)
		if (vdj == TRUE && exists("vdj.list")) {
			S4Vectors::metadata(dataset)$VDJ <- vdj.list
		}
		#
		if (verbose) cat("SingleCellExperiment object created with", ncol(dataset), "cells and", nrow(dataset), "genes", "\n")
		return(dataset)
		#
	} else {
		# Create ExpressionSet object (when output_type == "eset")
		if (verbose) cat("Creating ExpressionSet object directly", "\n")
		dataset <- new("ExpressionSet")
		
		# Set up expression data
		assayData.list <- list()
		assayData.list$exprs <- counts.all
		Biobase::assayData(dataset) <- as.environment(assayData.list)
		
		# Set up phenotype data
		Biobase::pData(dataset) <- metadata_df
		
		# Add multimodal data to assayData slots
		if (cite == TRUE && cite_ignore == FALSE && exists("counts.cite.all")) {
			Biobase::assayData(dataset)$Protein <- list()
			Biobase::assayData(dataset)$Protein["Counts"] <- list(Matrix::Matrix(counts.cite.all, sparse=TRUE))
		}
		#
		if (vdj == TRUE && exists("vdj.list")) {
			Biobase::assayData(dataset)$VDJ <- list()
			Biobase::assayData(dataset)$VDJ["VDJ"] <- list(vdj.list)
		}
		#		
		if (verbose) cat("ExpressionSet object created with", ncol(dataset), "cells and", nrow(dataset), "genes", "\n")
		return(dataset)
		#
	}
	#
}
#
