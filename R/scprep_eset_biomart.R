#' Function to annotate feature data (fData) in ExpressionSet object
#'
#' This function annotate feature data (fData) in ExpressionSet object with scRNA data
#' @param dataset Existing ExpressionSet object created with eset_build()
#' @param ensembl_target Species input for biomaRt feature annotation
#' @param reference Reference genome (on annotation.csv & SampleSheet.csv): "refdata-cellranger-GRCh38-3.0.0" or "GRCh38_pre-mRNA" for Ensembl 93; "refdata-gex-GRCh38-2020-A" or "GRCh38_pre-mRNA-2020-A" for Ensembl 98
#' @return ExpressionSet object containing annotated feature data in fData slot
#' @export
#' @examples
#' scprep_eset_biomart(dataset=ExpressionSet, ensembl_target="hsapiens_gene_ensembl", reference="refdata-cellranger-GRCh38-3.0.0")
#

#
scprep_eset_biomart <- function(
	dataset=dataset,
	ensembl_target=ensembl_target,
	reference=reference) {
	#
	cat("BiomaRt Feature Annotation", "\n")

	# Initialize featureData slot of ExpressionSet
	fData(dataset) <- data.frame(row.names=rownames(exprs(dataset)))

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
	} else if (reference == "refdata-gex-mm10-2020-A") {
		#
		cat("Retrieve Annotation for mmusculus Ensembl 98", "\n")			
		data(biomart_annotation_ens_98_mm10, package="scprep")
		biomart_annotation <- biomart_annotation_ens_98_mm10
		#
	}
	#
 	for (fdata in colnames(biomart_annotation)) {
 		fData(dataset)[[fdata]] <- biomart_annotation[intersect(rownames(biomart_annotation), rownames(exprs(dataset))), fdata]
 	}

	# Quantify Biotype information per cell and include in assayData slot
	cat("Quantify Biotype Annotation per Cell", "\n")
	#
	biotype.vars <- c("Biotype_Edit", "PC_Biotype")
	#
	for (biotype.var in biotype.vars) {
		#
		biotypes <- unique(fData(dataset)[[biotype.var]])
		#
		biotype.quant <- sapply(biotypes, function (biotype) {
			cat(paste("Biotype:", biotype), "\n")
			sapply(1:ncol(exprs(dataset)), function(cell) {
				return(sum(exprs(dataset)[rownames(exprs(dataset))[which(fData(dataset)[[biotype.var]] == biotype)], cell]))
			})
		})
		colnames(biotype.quant) <- biotypes
		rownames(biotype.quant) <- colnames(exprs(dataset))

		# Add biotype quantification to ExpressionSet pData slot
 		for (biotype in biotypes) {
 			dataset[[biotype]] <- biotype.quant[dataset$ID, biotype]
 		}
 	}

	# Quantify Biotype information per cell and include in assayData slot
	cat("Quantify Chromosomal Annotation per Cell", "\n")
	#
	chroms <- c(1:22, "X", "Y")
	#
	chr.quant <- sapply(chroms, function (chr) {
		cat(paste("Chromosome:", chr), "\n")
		sapply(1:ncol(exprs(dataset)), function(cell) {
			return(sum(exprs(dataset)[rownames(exprs(dataset))[which(fData(dataset)$Chr == chr)], cell]))
		})
	})
	colnames(chr.quant) <- chroms
	rownames(chr.quant) <- colnames(exprs(dataset))

	# Add biotype quantification to ExpressionSet pData slot
 	for (chr in chroms) {
 		dataset[[paste("Chr_", chr, sep="")]] <- chr.quant[dataset$ID, chr]
 	}
 	#
	return(dataset)
	#
}