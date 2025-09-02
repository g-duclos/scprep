#' Function to initiate single-cell data object
#'
#' This function builds an ExpressionSet, Seurat, or SingleCellExperiment object with scRNA data
#' @param sample_paths Character vector of paths to filtered_feature_bc_matrix.h5 for all samples in project
#' @param annotation Matrix of project-specific annotation.csv
#' @param vdj Data is multimodal RNA & VDJ
#' @param cite Data is multimodal RNA & Protein
#' @param cite_ignore FALSE to ignore; Default is TRUE, ignore the "Antibody Capture" slot in the h5 file input and only extract/store the "Gene Expression" slot; This is only relevant if your data contains both RNA and protein modalities and you want to ignore the protein modality
#' @param atac Data is multimodal RNA & ATAC
#' @return ExpressionSet object containing expression data for all samples & annotation in pData slot & UMIs/Genes per cell
#' @export
#' @examples
#' scprep_eset_build(sample_paths=sample_paths, annotation=annotation, file_type="h5", vdj=vdj, cite=cite, atac=atac)
#
scprep_eset_build <- function(
	sample_paths=sample_paths,
	annotation=annotation,
	file_type="h5",
	gene_id="symbol",
	vdj=vdj,
	cite=cite,
	atac=atac) {
	#
	cat("Build ExpressionSet", "\n")
	#
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
			cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
		} else if (cite == TRUE && cite_ignore == FALSE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
			cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			#
			counts.cite.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
			cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			colnames(counts.cite.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
		} else if (atac == TRUE && atac_ignore == FALSE && cite == FALSE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
			cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
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
			cat("Read Filtered GEX Matrix (ignoring multimodal data):", basename(sample_paths), "\n")
			cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
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
				cat(n.barcodes[i], "Cells Detected in", annotation$Sample_ID[i], "\n")
				return(sapply(1:ncol(annotation), function(x) {
					return(rep(annotation[i, x], n.barcodes[i]))
				}))
			}))
			colnames(annot.all) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
			rownames(annot.all) <- colnames(counts.all)
			#
			} else {
				cat("Matrix Sample Number Does NOT Match Sample Number Specified in Annotation!", "\n")
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
						cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
						#
					} else if (cite == TRUE) {
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
						cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						#
						counts.cite.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
						cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
						#
						mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
						colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						colnames(counts.cite.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						#
						cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
						#
					} else if (atac == TRUE & cite == FALSE) {
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
						cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
						#
						mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
						colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
						#
						cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
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
							cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						} else if (cite == TRUE) {
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
							cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							#
							counts.cite <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=TRUE)[["Antibody Capture"]]))
							cat("Read Filtered Surface Protein Matrix:", basename(sample_paths[i]), "\n")
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							colnames(counts.cite) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						} else if (atac == TRUE & cite == FALSE) {
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=use_names)[["Gene Expression"]]))
							cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
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
				cat("Read VDJ data:", vdj.file.name, "\n")
				vdj.list[as.character(basename(sample_paths))] <- list(read.csv(file.path(sample_paths, vdj.file.name)))
			} else {
				cat("VDJ data NOT found !!!", "\n")
			}
			#
		} else {
			#
			for (i in 1:length(sample_paths)) {
			#
				if (file.exists(file.path(sample_paths, vdj.file.name)) == TRUE) {
					cat("Read VDJ data:", vdj.file.name, "\n")
					#
					vdj.list[as.character(basename(sample_paths[i]))] <- list(read.csv(file.path(sample_paths[i], vdj.file.name)))
					#
				} else {
					cat("VDJ data NOT found !!!", "\n")
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
	cat("Calculate UMIs per cell", "\n")
	metadata_df$UMIs <- colSums(counts.all)

	cat("Calculate Genes per cell", "\n")
	dataset$Genes <- unlist(lapply(1:ncol(exprs(dataset)), function(x) {
		return(length(rownames(exprs(dataset))[which(exprs(dataset)[,x] >= 1)]))
	}))

	# Define additional assayData slots
	assayData(dataset)$Params <- list()
	assayData(dataset)$Seurat <- list()
	assayData(dataset)$QC <- list()
	assayData(dataset)$Secondary <- list()
	#
	return(dataset)
	#
}
#
