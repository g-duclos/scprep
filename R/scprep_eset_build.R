#' Function to initiate ExpressionSet object
#'
#' This function builds an ExpressionSet object with scRNA data
#' @param sample_paths Character vector of paths to filtered_feature_bc_matrix.h5 for all samples in project
#' @param annotation Matrix of project-specific annotation.csv
#' @param vdj Data is multimodal RNA & VDJ
#' @param cite Data is multimodal RNA & Protein
#' @param atac Data is multimodal RNA & ATAC
#' @return ExpressionSet object containing expression data for all samples & annotation in pData slot & UMIs/Genes per cell
#' @export
#' @examples
#' scprep_eset_build(sample_paths=sample_paths, annotation=annotation, cite=cite)
#
scprep_eset_build <- function(
	sample_paths=sample_paths,
	annotation=annotation,
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
		dir_counts <- file.path(sample_paths, counts.file.name)
		#
		if (cite == FALSE & atac == FALSE) {
			#
			counts.all <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=FALSE))
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
		} else if (cite == TRUE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
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
		} else if (atac == TRUE & cite == FALSE) {
			#
			counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
			cat("Read Filtered GEX Matrix:", basename(sample_paths), "\n")
			#
			mat.barcodes <- unlist(strsplit(colnames(counts.all),"-"))
			colnames(counts.all) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
			#
			cat(ncol(counts.all), "Cells Detected in", basename(sample_paths[i]), "\n")
			#
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
				dir_counts <- file.path(sample_paths[i], counts.file.name)
				#
				if (i == 1) {
					#
					if (cite == FALSE & atac == FALSE) {
						#
						counts.all <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=FALSE))
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
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
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
						counts.all <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
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
							counts <- suppressMessages(Seurat::Read10X_h5(dir_counts, use.names=FALSE))
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
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
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
							counts <- suppressMessages(as.matrix(Seurat::Read10X_h5(dir_counts, use.names=FALSE)[["Gene Expression"]]))
							cat("Read Filtered GEX Matrix:", basename(sample_paths[i]), "\n")
							#
							mat.barcodes <- unlist(strsplit(colnames(counts),"-"))
							colnames(counts) <- paste(mat.barcodes[seq(1,length(mat.barcodes),2)], "-", i, sep="")
							#
							cat(ncol(counts), "Cells Detected in", basename(sample_paths[i]), "\n")
							#
						}
						#
						annot = sapply(1:ncol(annotation), function(x) {
							return(rep(annotation[i, x], ncol(counts)))
							})
						colnames(annot) <- c("Sample", colnames(annotation)[2:ncol(annotation)])
						rownames(annot) <- colnames(counts)
						#
						counts.all = cbind(counts.all, counts)
						#
						if (cite == TRUE) {
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
	#
	# Initiate ExpressionSet
	dataset <- new("ExpressionSet")
	#
	pData(dataset) <- data.frame(
    	ID = colnames(counts.all),
		stringsAsFactors = FALSE)
	#
	assayData.list <- list()
	#
	assayData.list$exprs <- matrix(
		data = 0L, nrow = nrow(counts.all), ncol = ncol(counts.all),
    	dimnames = list(rownames(counts.all), colnames(counts.all)))
	#
    colnames(assayData.list$exprs) <- colnames(counts.all)
   	#
    assayData.list$exprs <- counts.all
	#
	dataset$ID <- colnames(assayData.list$exprs)
	rownames(pData(dataset)) <- dataset$ID
	#
	assayData(dataset) <- as.environment(assayData.list)
	
	#	
	if (cite == TRUE) {
		# Add CITE-Seq matrix to assayData slot
		assayData(dataset)$Protein <- list()
		assayData(dataset)$Protein["Counts"] <- list(Matrix::Matrix(counts.cite.all, sparse=TRUE))
	}

	#	
	if (vdj == TRUE) {
		# Add CITE-Seq matrix to assayData slot
		assayData(dataset)$VDJ <- list()
		assayData(dataset)$VDJ["VDJ"] <- list(vdj.list)
	}

	# Add annotation variables to pData slot
	for (i in 1:ncol(annot.all)) {
		dataset[[colnames(annot.all)[i]]] <- as.factor(annot.all[,i])
	}
	
	# Total counts per sampl
	cat("Calculate UMIs per cell", "\n")
	dataset$UMIs <- colSums(counts.all)

	# Total genes per sample
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
