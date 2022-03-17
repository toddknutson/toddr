#' RNAseq normalization methods
#'
#' Perform various RNAseq normalization methods on a count table.
#'
#' @param data [Required] Integer based count table. Rows represent genes (features), columns represent samples (observations). 
#' @param gene_lengths [Default: NULL] The nucleotide count (length) of each gene. Required for TPM method.
#' @param method [Default: c("DESeq2_norm_counts", "edgeR_TMM", "TPM")] Method of RNAseq normalization.
#' @param gene_names [Default: rownames(data)] Character vector of gene names.
#' @param gene_names_col_id [Default: "ensembl_id"] String column name that should be used for output table.
#' @param sample_names [Default: colnames(data)] Character vector of sample names.
#'
#' 
#'
#' @return A tibble::tibble data frame object of normalized counts. The first column will be a character type of gene names. Remaining columns will represent sample names.
#'
#'
#'
#' Good info about normalization methods:
#' [https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)
#' [https://www.biostars.org/p/317701](https://www.biostars.org/p/317701)
#' [https://www.biostars.org/p/317417](https://www.biostars.org/p/317417)
#' [https://www.biostars.org/p/210447](https://www.biostars.org/p/210447)
#'
#' @examples
#' set.seed(1)
#' nrow = 10000
#' ncol = 6
#' min = 0
#' max = 1000
#' gene_lengths = sample(400:10000, nrow, TRUE)
#' my_counts <- matrix(sample(min:max, nrow*ncol, TRUE), 
#'     nrow = nrow, ncol = ncol, 
#'     dimnames = list(paste0("gene", 1:nrow), paste0("sample", 1:ncol)))
#' 
#' my_counts_norm <- rnaseq_norm(my_counts, method = "DESeq2_norm_counts")
#' my_counts_norm <- rnaseq_norm(my_counts, method = "edgeR_TMM")
#' my_counts_norm <- rnaseq_norm(my_counts, gene_lengths = gene_lengths, method = "TPM")
#'
#'
#' @export
rnaseq_norm <- function(data, gene_lengths = NULL, method = c("DESeq2_norm_counts", "edgeR_TMM", "TPM"), gene_names = rownames(data), gene_names_col_id = "ensembl_id", sample_names = colnames(data)) {
	stopifnot(!is.null(data))
	stopifnot(!is.null(method))
	# DESeq2_norm_counts
	if (method == "DESeq2_norm_counts") {
		counts <- data
		counts <- counts %>% as.matrix(.)
		rownames(counts) <- gene_names
		dds <- DESeqDataSetFromMatrix(countData = counts,
									  colData = tibble(sample_names = sample_names),
									  design = ~ 1)
		dds <- estimateSizeFactors(dds)
		counts_norm <- t(t(counts) / sizeFactors(dds))
		# Add a small value to all to make sure there are no zeros. This will allow for log transformation.
		#counts_norm <- counts_norm + 0.5
		# Convert to tibble
		counts_norm <- counts_norm %>%
							as.data.frame(.) %>%
							rownames_to_column(., var = gene_names_col_id) %>%
							remove_rownames(.) %>%
							as_tibble(.)
		return(counts_norm)
	}
	# edgeR_TMM
	if (method == "edgeR_TMM") {
		# https://www.biostars.org/p/317701/
		counts <- data
		dge <- DGEList(counts, group = NULL)
		dge <- calcNormFactors(dge, method = "TMM")
		tmm <- cpm(dge)
		rownames(tmm) <- gene_names
		# Convert to tibble
		tmm <- tmm %>%
					as.data.frame(.) %>%
					rownames_to_column(., var = gene_names_col_id) %>%
					remove_rownames(.) %>%
					as_tibble(.)
		return(tmm)
	}
	# TPM
	if (method == "TPM") {
	    stopifnot(!is.null(gene_lengths))
		# https://support.bioconductor.org/p/91218/
		counts <- data
		x <- (counts / gene_lengths)
		tpm <- t(t(x) * 1e6 / colSums(x))
		rownames(tpm) <- gene_names
		tpm <- tpm %>%
					as.data.frame(.) %>%
					rownames_to_column(., var = gene_names_col_id) %>%
					remove_rownames(.) %>%
					as_tibble(.)
		return(tpm)
	}
}


