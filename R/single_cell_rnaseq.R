#!/usr/bin/env Rscript



























# ---------------------------------------------------------------------
# Merge clusters based on DE analysis
# ---------------------------------------------------------------------


#' @export
tk_get_medoids <- function(m, cluster) {
      # m:
      # cluster = columns
      # PCs = rows
      # cluster: vector of class labels
      cluster_names <- unique(cluster)
      cluster_num <- length(unique(cluster))
      
      all_cluster_medoid <- matrix(0L, nrow = nrow(m), ncol = cluster_num)
      # For each cluster
      for (i in 1:cluster_num) {
            curr_m <- m[, cluster == cluster_names[i]]
            # For each row (PC)
            curr_cluster_medoid <- rep(0, dim(curr_m)[1])
            for (j in seq_along(curr_cluster_medoid)) {
                  d <- dist(curr_m[j, ], upper = TRUE)
                  d <- as.matrix(d)
                  s <- apply(d, 2, sum)
                  # Which cell in this cluster
                  w <- which.min(s)
                  curr_cluster_medoid[j] <- curr_m[j, w]
            }
            all_cluster_medoid[, i] <- curr_cluster_medoid
      }
      colnames(all_cluster_medoid) <- cluster_names
      rownames(all_cluster_medoid) <- rownames(m)
      return(all_cluster_medoid)
}


# https://stackoverflow.com/questions/20396582/order-a-mixed-vector-numbers-with-letters
#' @export
tk_multi_mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
        lapply(list(...), function(l){
            if(is.character(l)){
                factor(l, levels=gtools::mixedsort(unique(l)))
            } else {
                l
            }
        }),
        list(na.last = na.last, decreasing = decreasing)
    ))
}





#' @export
tk_parallel_de_pairwise <- function(k, seurat_object, curr_clusters, curr_number_of_clusters, combos, lfc_threshold, merge_data_list) {
	# Start running function
	curr_cluster_a <- combos[1, k]
	curr_cluster_b <- combos[2, k]
	name <- paste0("cluster_", curr_cluster_a, "_vs_", curr_cluster_b)
	# If the DE comparison was done previously, and is not changed here, don't repeat.
	if (any(names(merge_data_list) %in% name)) {
		# Get previous results stored in "merge_data_list"
		prev_merge_data_list_idx <- which(names(merge_data_list) == name)
		print(paste0(names(merge_data_list)[prev_merge_data_list_idx], " comparison already done, skipping."))
		results <- merge_data_list[[prev_merge_data_list_idx]]
		#results <- NULL
	} else {
		print(paste0(name, " running DE analysis."))
		# Find the DE genes for a new comparison
		#curr_merge_data_list_idx <- (first_avail_merge_data_list_idx + k) - 1

		# Run DE analysis
		de <- FindMarkers(object = seurat_object, slot = "data", ident.1 = combos[1, k], ident.2 = combos[2, k], min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.1)
		# Extract significant DE genes
		lfc_col <- which(grepl("log2FC", colnames(de)))
		adj_p_col <- which(grepl("p_val_adj", colnames(de)))
		de_signif <- de[abs(de[, lfc_col]) >= lfc_threshold & de[, adj_p_col] <= 0.01, ]
		de_signif <- de_signif[order(de_signif[, lfc_col], decreasing = TRUE), ]
		de <- rownames_to_column(as.data.frame(de), var = "symbol")
		de_signif <- rownames_to_column(as.data.frame(de_signif), var = "symbol")
		# How many genes are significant?
		number_of_signif_genes <- dim(de_signif)[1]
		results <- list(name = name, cluster_a = curr_cluster_a, cluster_b = curr_cluster_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
	}
	return(results)
}




#' @export
tk_parallel_de_one_vs_all <- function(m, seurat_object, lfc_threshold) {
	curr_cluster_a <- paste0("cluster_", levels(factor(seurat_object@active.ident))[m])
	name <- paste0(curr_cluster_a, "_vs_all")
	print(paste0(name, " DE analysis."))
	de_cluster_vs_all <- FindMarkers(object = seurat_object, slot = "data", ident.1 = levels(factor(seurat_object@active.ident))[m], min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.1)
	# Extract significant DE genes
	lfc_col <- which(grepl("log2FC", colnames(de_cluster_vs_all)))
	adj_p_col <- which(grepl("p_val_adj", colnames(de_cluster_vs_all)))
	de_cluster_vs_all_signif <- de_cluster_vs_all[abs(de_cluster_vs_all[, lfc_col]) >= lfc_threshold & de_cluster_vs_all[, adj_p_col] <= 0.01, ]
	de_cluster_vs_all_signif <- de_cluster_vs_all_signif[order(de_cluster_vs_all_signif[, lfc_col], decreasing = TRUE), ]
	de_cluster_vs_all <- rownames_to_column(as.data.frame(de_cluster_vs_all), var = "symbol")
	de_cluster_vs_all_signif <- rownames_to_column(as.data.frame(de_cluster_vs_all_signif), var = "symbol")
	de_cluster_vs_all_signif_number_of_signif_genes <- dim(de_cluster_vs_all_signif)[1]
	results <- list(name = name, cluster_a = curr_cluster_a, cluster_b = "all_other_cells", de = de_cluster_vs_all, de_signif = de_cluster_vs_all_signif, number_of_signif_genes = de_cluster_vs_all_signif_number_of_signif_genes)
	return(results)
}




#' @export
tk_parallel_de_plots_pairwise <- function(r, seurat_object = seurat_object, merge_data_list = merge_data_list, mapping = mapping, final_df = final_df, out_dir2 = out_dir2, umap_plot_list = umap_plot_list) {
	# DE LISTS PAIRWISE
	ugly_name <- paste0("cluster_", final_df$comparison_a[r], "_vs_", final_df$comparison_b[r])
	clean_name <- paste0("cluster_", mapping$de_merge_final[mapping$de_merge_orig == as.character(final_df$comparison_a[r])], "_vs_", mapping$de_merge_final[mapping$de_merge_orig == as.character(final_df$comparison_b[r])])
	write_tsv(merge_data_list[[ugly_name]]$de, paste0(out_dir2, "/de_pairwise/", clean_name, "_de.txt"))
	write_tsv(merge_data_list[[ugly_name]]$de_signif, paste0(out_dir2, "/de_pairwise/", clean_name, "_de_signif.txt"))
	# GENE PLOTS PAIRWISE
	orig_wd <- getwd()
	setwd(paste0(out_dir2, "/de_pairwise/umap_gene_expression"))
	geneset_up <- merge_data_list[[ugly_name]]$de_signif$symbol[merge_data_list[[ugly_name]]$de_signif$avg_log2FC > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(clean_name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot_list[[2]], type = "feature_plot")
	geneset_dn <- merge_data_list[[ugly_name]]$de_signif$symbol[merge_data_list[[ugly_name]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(clean_name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot_list[[2]], type = "feature_plot")
	setwd(orig_wd)
	# Make Excel file
	excel_list_of_df_de <- merge_data_list[[ugly_name]]$de
	excel_list_of_df_de_signif <- merge_data_list[[ugly_name]]$de_signif
	signif_genes_across_all_clusters <- merge_data_list[[ugly_name]]$de_signif$symbol
	results <- list(ugly_name = ugly_name, clean_name = clean_name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}



#' @export
tk_parallel_de_plots_one_vs_all <- function(s, seurat_object = seurat_object, merge_data_list = merge_data_list, mapping = mapping, out_dir2 = out_dir2, umap_plot_list = umap_plot_list) {
	# DE LISTS ALL_VS_ONE
	ugly_name <- paste0("cluster_", mapping$de_merge_orig[s], "_vs_all")
	clean_name <- paste0("cluster_", mapping$de_merge_final[s], "_vs_all")
	write_tsv(merge_data_list[[ugly_name]]$de, paste0(out_dir2, "/de_one_vs_all/", clean_name, "_de.txt"))
	write_tsv(merge_data_list[[ugly_name]]$de_signif, paste0(out_dir2, "/de_one_vs_all/", clean_name, "_de_signif.txt"))
	# GENE PLOTS ALL_VS_ONE
	orig_wd <- getwd()
	setwd(paste0(out_dir2, "/de_one_vs_all/umap_gene_expression"))
	geneset_up <- merge_data_list[[ugly_name]]$de_signif$symbol[merge_data_list[[ugly_name]]$de_signif$avg_log2FC > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(clean_name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot_list[[2]], type = "feature_plot")
	geneset_dn <- merge_data_list[[ugly_name]]$de_signif$symbol[merge_data_list[[ugly_name]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(clean_name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot_list[[2]], type = "feature_plot")
	setwd(orig_wd)
	# Make Excel file
	excel_list_of_df_de <- merge_data_list[[ugly_name]]$de
	excel_list_of_df_de_signif <- merge_data_list[[ugly_name]]$de_signif
	signif_genes_across_all_clusters <- merge_data_list[[ugly_name]]$de_signif$symbol
	results <- list(ugly_name = ugly_name, clean_name = clean_name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}




#' @export
tk_cluster_merge <- function(seurat_object, seurat_object_name, lfc_threshold, num_genes_diff_between_clusters_threshold) {
	print(seurat_object_name)
	# Get the matrix and labels from the clustering results
	# These class labels imported as 1-indexed because they were converted by different function (above) after clustering
	cluster_cell_class <- seurat_object@meta.data[, grepl("res.", colnames(seurat_object@meta.data)), drop = FALSE]
	cluster_res_name <- str_replace_all(colnames(cluster_cell_class), fixed("res."), "res_")
	cell_ids <- rownames(cluster_cell_class)
	return_list <- list()
	# For each clustering resolution
	for (j in seq_along(cluster_res_name))  {
		# Create dir for all results
		out_dir2 <- cluster_res_name[j]
		if (!dir.exists(out_dir2)) {dir.create(out_dir2, recursive = TRUE)}
		print(cluster_res_name[j])
		# Factor the current clustering resolution cell classes (should start with 1-index based)
		class_orig <- factor(cluster_cell_class[, j])
		names(class_orig) <- cell_ids
		# Re-assign the class names to @ident slot with current cell classes
		seurat_object@active.ident <- class_orig
		iteration_counter <- 0
		continue_merging_clusters <- TRUE
		# List to store all DE comparison data
		#merge_data_list <- vector("list", 1000000)
		merge_data_list <- list()
		iteration_data_list <- vector("list", 1000000)
		# While the variable above is true, continue to merge clusters, until unique clusters are found (based on DE analysis)
		while (continue_merging_clusters == TRUE) {
			iteration_counter <- iteration_counter + 1
			print(cluster_res_name[j])
			print(paste0("iteration: ", iteration_counter))
			#first_avail_merge_data_list_idx <- sum(!sapply(merge_data_list, is.null)) + 1
			curr_clusters <- gtools::mixedsort(levels(seurat_object@active.ident))
			curr_number_of_clusters <- length(curr_clusters)
			# Generate pairwise cluster numbering for all cluster comparisons
			combos <- combn(curr_clusters, 2)
			# Figure out the right assay to use
			# If this is an individual library (single sample) and the SCT has been run, the assay should still be the default SCT (leave it). The slot should be "data" by default, which are "SCT corrected" log normalized UMI counts.
			# Else, if integrated assay, switch to RNA:data assay:slot for FeaturePlots
			orig_assay <- DefaultAssay(object = seurat_object)
			if (DefaultAssay(object = seurat_object) == "integrated") {
				DefaultAssay(object = seurat_object) <- "RNA"
			}
			# Run all pairwise cluster comparisons, find DE genes
			k_vect <- 1:dim(combos)[2]
			# Help with parallel lapply:
			# https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
			out <- mclapply(X = k_vect, FUN = tk_parallel_de_pairwise, seurat_object = seurat_object, curr_clusters = curr_clusters, curr_number_of_clusters = curr_number_of_clusters, combos = combos, lfc_threshold = lfc_threshold, merge_data_list = merge_data_list,
				mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
			DefaultAssay(object = seurat_object) <- orig_assay
			# name the list elements
			names(out) <- sapply(out, function(x) x[["name"]])
			number_of_sig_genes_for_each_comparison <- unname(sapply(out, function(x) x[["number_of_signif_genes"]]))
			# Are there any pairwise comparisons that have less than "num_genes_diff_between_clusters_threshold" DE genes between them? If so, merge clusters.
			combos2 <- as.data.frame(t(combos))
			combos2 <- data.frame(combos2, number_of_sig_genes_for_each_comparison)
			colnames(combos2) <- c("comparison_a", "comparison_b", "number_of_sig_genes_for_each_comparison")
			# Update the merge_data_list
			# Find length of prev merge data list
			merge_data_list_length <- length(merge_data_list)
			# append lists
			for (i in seq_along(out)) {
				idx <- merge_data_list_length + i
				merge_data_list[[idx]] <- out[[i]]
				names(merge_data_list)[idx] <- out[[i]]$name
			}
			# Get rid of any duplicate entries in the appended merge_data_list
			merge_data_list <- merge_data_list[!duplicated(names(merge_data_list))]
			if (any(number_of_sig_genes_for_each_comparison < num_genes_diff_between_clusters_threshold)) {
				print(paste0("There are cluster comparisons with less than ", num_genes_diff_between_clusters_threshold, " DE genes between clusters -- starting to merge them."))
				# Which clusters to merge?
				# Which comparisons have too few DE genes, and will be merged?
				# Calculate the cluster medoids in PCA space, only for the cells that have "too few" DE genes between them
				if (DefaultAssay(seurat_object) == "RNA" & "pca_rna" %in% names(seurat_object@reductions)) {
					medoids <- tk_get_medoids(m = t(seurat_object@reductions$pca_rna@cell.embeddings), cluster = seurat_object@active.ident)
				} else {
					medoids <- tk_get_medoids(m = t(seurat_object@reductions$pca@cell.embeddings), cluster = seurat_object@active.ident)
				}
				medoids_dist <- as.matrix(dist(t(medoids), diag = TRUE, upper = FALSE))
				medoids_dist[upper.tri(medoids_dist, diag = TRUE)] <- NA
				medoids_dist_df <- as.data.frame(medoids_dist)
				medoids_dist_df <- rownames_to_column(medoids_dist_df, var = "cluster")
				medoids_dist_df <- gather(medoids_dist_df, cluster2, value, -cluster)
				medoids_dist_df <- medoids_dist_df[!is.na(medoids_dist_df$value), ]
				medoids_dist_df$new_order1 <- NA
				medoids_dist_df$new_order2 <- NA
				for (v in 1:dim(medoids_dist_df)[1]) {
					new_order <- gtools::mixedsort(levels(factor(c(medoids_dist_df$cluster[v], medoids_dist_df$cluster2[v]))))
					medoids_dist_df$new_order1[v] <- new_order[1]
					medoids_dist_df$new_order2[v] <- new_order[2]
				}
				medoids_dist_df <- medoids_dist_df[tk_multi_mixedorder(medoids_dist_df$new_order1, medoids_dist_df$new_order2), ]
				merged_medoids <- merge(combos2, medoids_dist_df, by.x = c("comparison_a", "comparison_b"), by.y = c("new_order1", "new_order2"))
				# Subset list to include only above clusters available for merging.
				merged_medoids <- merged_medoids[merged_medoids$number_of_sig_genes_for_each_comparison < num_genes_diff_between_clusters_threshold, ]
				merged_medoids <- merged_medoids[order(merged_medoids$value, decreasing = FALSE), ]
				# Get the top cluster comparison for merging
				merge_small <- as.character(merged_medoids$comparison_a[1])
				merge_large <- as.character(merged_medoids$comparison_b[1])
				# Re-create cell new merged classes
				curr_classes <- seurat_object@active.ident	
				new_classes <- as.character(curr_classes)
				new_classes[curr_classes == merge_small] <- paste0(merge_small, ".", merge_large)
				new_classes[curr_classes == merge_large] <- paste0(merge_small, ".", merge_large)
				names(new_classes) <- names(curr_classes)
				new_classes <- factor(new_classes)
				new_classes <- droplevels(new_classes)
				# Re-assign the class names to @ident slot
				seurat_object@active.ident <- new_classes
				# Repeat pairwise again
				continue_merging_clusters <- TRUE
				clusters_merged <- c(merge_small, merge_large)
			} else {
				curr_classes <- seurat_object@active.ident
				final_classes <- curr_classes
				print(paste0("All clusters have more than ", num_genes_diff_between_clusters_threshold, " DE genes between them. No merging."))
				continue_merging_clusters <- FALSE
				clusters_merged <- "none"
			}
			names(iteration_data_list)[iteration_counter] <- paste0("iteration_", iteration_counter)
			iteration_data_list[[iteration_counter]] <- list(class_orig = class_orig, pre_iteration_classes = curr_classes, post_iteration_classes = seurat_object@active.ident, pairwise_signif_genes = combos2, were_clusters_merged = continue_merging_clusters, clusters_merged = clusters_merged)
		}	
		# MERGING HAS STOPPED, run the following code
		# Redo DE testing for cluster "a" vs all other cells (all other clusters combined)
		# Re-assign the identity factor 
		seurat_object@active.ident <- final_classes		
		curr_number_of_clusters_final <- length(levels(factor(seurat_object@active.ident)))
		# Run all pairwise cluster comparisons, find DE genes
		# Set variables to global env, so they can be used in mclapply parallel function below
		first_avail_merge_data_list_idx <- sum(!sapply(merge_data_list, is.null)) + 1
		# Figure out the right assay to use
		# If this is an individual library (single sample) and the SCT has been run, the assay should still be the default SCT (leave it). The slot should be "data" by default, which are "SCT corrected" log normalized UMI counts.
		# Else, if integrated assay, switch to RNA:data assay:slot for FeaturePlots
		orig_assay <- DefaultAssay(object = seurat_object)
		if (DefaultAssay(object = seurat_object) == "integrated") {
			DefaultAssay(object = seurat_object) <- "RNA"
		}
		m_vect <- 1:curr_number_of_clusters_final
		out <- mclapply(X = m_vect, FUN = tk_parallel_de_one_vs_all, seurat_object = seurat_object, lfc_threshold = lfc_threshold,
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		DefaultAssay(object = seurat_object) <- orig_assay
		# name the list elements
		names(out) <- sapply(out, function(x) x[["name"]])
		# update merge_data_list	
		merge_data_list_length <- length(merge_data_list)
		# append lists
		for (i in seq_along(out)) {
			idx <- merge_data_list_length + i
			merge_data_list[[idx]] <- out[[i]]
			names(merge_data_list)[idx] <- out[[i]]$name
		}
		# Update object with final class names 
		# Rename the final clusters to be regular integers, for plotting and file output
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_orig")]] <- as.character(final_classes)
		# Rename the final clusters to be regular integers, for plotting and file output
		# Make them zero-based index
		cluster_cell_class_factor_final <- factor(as.numeric(final_classes) - 1)
		names(cluster_cell_class_factor_final) <- names(final_classes)
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final")]] <- as.character(cluster_cell_class_factor_final)	
		# Map the ugly names to final names
		mapping_full <- data.frame(de_merge_orig = final_classes, de_merge_final = cluster_cell_class_factor_final)
		mapping <- unique(mapping_full)
		mapping <- mapping %>% dplyr::arrange(de_merge_final)
		write_tsv(mapping, paste0(out_dir2, "/mapping_cluster_labels.txt"), col_names = TRUE)
		rownames(mapping) <- NULL
		# Export data here
		# Get rid of extra NA values
		#merge_data_list <- merge_data_list[!sapply(merge_data_list, is.null)]
		iteration_data_list <- iteration_data_list[!sapply(iteration_data_list, is.null)]
		final_cluster_labels <- data.frame(cell_id = names(cluster_cell_class_factor_final), clust_orig = class_orig, clust_de_merge_orig = final_classes, clust_de_merge_final = cluster_cell_class_factor_final)
		write_tsv(final_cluster_labels, paste0(out_dir2, "/final_cluster_labels.txt"), col_names = TRUE)
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# UMAP PLOTS
		# Plot original UMAP
		umap_plot_list <- list()
		# Re-assign the identity factor
		# Make this a nicely sorted factor of integers, zero-based
		class_orig_factor <- factor(as.numeric(class_orig) - 1)
		names(class_orig_factor) <- names(class_orig)
		seurat_object@active.ident <- class_orig_factor
		# get centers for labels
		if (DefaultAssay(seurat_object) == "RNA" & "umap_rna" %in% names(seurat_object@reductions)) {
			plot_data <- as.data.frame(seurat_object@reductions$umap_rna@cell.embeddings) %>%
				rename(UMAP_1 = rnaUMAP_1, UMAP_2 = rnaUMAP_2)
		} else {
			plot_data <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
		}
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[1]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nInitial pre-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[1] <- paste0(cluster_res_name[j], "_before_de_merge")
		# Plot final UMAP
		# Re-assign the identity factor
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# get centers for labels
		if (DefaultAssay(seurat_object) == "RNA" & "umap_rna" %in% names(seurat_object@reductions)) {
			plot_data <- as.data.frame(seurat_object@reductions$umap_rna@cell.embeddings) %>%
				rename(UMAP_1 = rnaUMAP_1, UMAP_2 = rnaUMAP_2)
		} else {
			plot_data <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
		}
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[2]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nFinal post-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[2] <- paste0(cluster_res_name[j], "_final")
		pdf(paste0(out_dir2, "/umap_final_de_merge.pdf"), width = 6, height = 6)
		print(umap_plot_list[[2]])
		dev.off()
		pdf(paste0(out_dir2, "/umap_original.pdf"), width = 6, height = 6)
		print(umap_plot_list[[1]])
		dev.off()
# 		saveRDS(seurat_object, paste0(out_dir2, "/seurat_object_after_cluster_merging.rds"))
# 		saveRDS(merge_data_list, paste0(out_dir2, "/merge_data_list.rds"))
# 		saveRDS(iteration_data_list, paste0(out_dir2, "/iteration_data_list.rds"))
# 		saveRDS(mapping, paste0(out_dir2, "/mapping.rds"))
# 		saveRDS(mapping_full, paste0(out_dir2, "/mapping_full.rds"))
# 		saveRDS(final_cluster_labels, paste0(out_dir2, "/final_cluster_labels.rds"))
# 		saveRDS(cluster_cell_class_factor_final, paste0(out_dir2, "/cluster_cell_class_factor_final.rds"))
# 		saveRDS(final_classes, paste0(out_dir2, "/final_classes.rds"))
# 		saveRDS(cluster_cell_class_factor_final, paste0(out_dir2, "/cluster_cell_class_factor_final.rds"))
# 		saveRDS(class_orig, paste0(out_dir2, "/class_orig.rds"))
# 		saveRDS(umap_plot_list, paste0(out_dir2, "/umap_plot_list.rds"))
		# GET PAIRWISE PLOTS AND LISTS
		if (!dir.exists(paste0(out_dir2, "/de_pairwise/umap_gene_expression"))) {dir.create(paste0(out_dir2, "/de_pairwise/umap_gene_expression"), recursive = TRUE)}
		# Figure out the right assay to use
		# If this is an individual library (single sample) and the SCT has been run, the assay should still be the default SCT (leave it). The slot should be "data" by default, which are "SCT corrected" log normalized UMI counts.
		# Else, if integrated assay, switch to RNA:data assay:slot for FeaturePlots
		orig_assay <- DefaultAssay(object = seurat_object)
		if (DefaultAssay(object = seurat_object) == "integrated") {
			DefaultAssay(object = seurat_object) <- "RNA"
		}
		# get the last iteration
		final_df <- iteration_data_list[[length(iteration_data_list)]]$pairwise_signif_genes
		r_vect <- 1:dim(final_df)[1]
		out <- mclapply(X = r_vect, FUN = tk_parallel_de_plots_pairwise, seurat_object = seurat_object, merge_data_list = merge_data_list, mapping = mapping, final_df = final_df, out_dir2 = out_dir2, umap_plot_list = umap_plot_list,
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		# Switch default back to original
		DefaultAssay(object = seurat_object) <- orig_assay
		names(out) <- sapply(out, function(x) x[["clean_name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# The openxlsx package is using an old zip method, so you might get warnings (but should still work)
		# Note: zip::zip() is deprecated, please use zip::zipr() instead
		# https://github.com/awalker89/openxlsx/issues/454
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")		
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = paste0(out_dir2, "/de_pairwise.xlsx"))
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = paste0(out_dir2, "/de_pairwise_signif.xlsx"))
		
		# HEATMAP PAIRWISE		
		if (!dir.exists(paste0(out_dir2, "/heatmaps"))) {dir.create(paste0(out_dir2, "/heatmaps"), recursive = TRUE)}
		orig_wd <- getwd()
		setwd(paste0(out_dir2, "/heatmaps"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_pairwise_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(orig_wd)
		
		# GET ONE_VS_ALL PLOTS AND LISTS
		if (!dir.exists(paste0(out_dir2, "/de_one_vs_all/umap_gene_expression"))) {dir.create(paste0(out_dir2, "/de_one_vs_all/umap_gene_expression"), recursive = TRUE)}
		# Figure out the right assay to use
		# If this is an individual library (single sample) and the SCT has been run, the assay should still be the default SCT (leave it). The slot should be "data" by default, which are "SCT corrected" log normalized UMI counts.
		# Else, if integrated assay, switch to RNA:data assay:slot for FeaturePlots
		orig_assay <- DefaultAssay(object = seurat_object)
		if (DefaultAssay(object = seurat_object) == "integrated") {
			DefaultAssay(object = seurat_object) <- "RNA"
		}		
		s_vect <- 1:curr_number_of_clusters_final
		out <- mclapply(X = s_vect, FUN = tk_parallel_de_plots_one_vs_all, seurat_object = seurat_object, merge_data_list = merge_data_list, mapping = mapping, out_dir2 = out_dir2, umap_plot_list = umap_plot_list,
			mc.cores = min(16, as.integer(system("echo $THREADS", intern = TRUE))))
		# Switch default back to original
		DefaultAssay(object = seurat_object) <- orig_assay
		names(out) <- sapply(out, function(x) x[["clean_name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")		
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = paste0(out_dir2, "/de_one_vs_all.xlsx"))
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = paste0(out_dir2, "/de_one_vs_all_signif.xlsx"))
		# HEATMAP ALL_VS_ONE
		if (!dir.exists(paste0(out_dir2, "/heatmaps"))) {dir.create(paste0(out_dir2, "/heatmaps"), recursive = TRUE)}
		orig_wd <- getwd()
		setwd(paste0(out_dir2, "/heatmaps"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_one_vs_all_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(orig_wd)

		# Get the final iteration counter for plotting
		iteration_counter_final <- iteration_counter - 1
		pdf(paste0(out_dir2, "/umap_original_and_final_de_merge.pdf"), width = 10, height = 6)
		grid.arrange(
			grobs = umap_plot_list,
			nrow = 1,
			top = paste0("UMAP plots, Pre- and Post-DE merging.\n", iteration_counter_final , " cluster merging steps required."))
		dev.off()
		# Export results to list of list
		curr_return_list <- list(seurat_object = seurat_object, 
			final_cluster_labels = final_cluster_labels, 
			iteration_data_list = iteration_data_list, 
			merge_data_list = merge_data_list, 
			mapping = mapping, 
			mapping_full = mapping_full, 
			lfc_threshold = lfc_threshold, 
			num_genes_diff_between_clusters_threshold = num_genes_diff_between_clusters_threshold,
			umap_plot_list = umap_plot_list)
		return_list[[j]] <- curr_return_list
		names(return_list)[j] <- cluster_res_name[j]
	}
	return(return_list)
}







# ---------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------


#' @export
tk_cluster_heatmap <- function(seurat_object, gene_names, filename_prefix, active.ident, heatmap_title = NULL) {
	# NOTE: the seurat_object@active.ident must have the correct cluster assignment labels!
	# The DefaultAssay(seurat_object) is used, so set before running this function.
	# Get the unique set of genes
	gene_names <- unique(gene_names)
	# Use the scaled data from seurat object
	# matrix: rows = genes, cols = cells
	
	# This will work for "individual" libraries
	matrix <- GetAssayData(object = seurat_object, slot = "scale.data", assay = DefaultAssay(seurat_object))
	if (all(dim(matrix) == c(0, 0))) {
		# scale.data is empty -- create it
		print("Scaling data for heatmap")
		seurat_object_temp2 <- seurat_object
		seurat_object_temp2 <- ScaleData(seurat_object_temp2, verbose = FALSE)
		matrix <- GetAssayData(object = seurat_object_temp2, slot = "scale.data", assay = DefaultAssay(seurat_object_temp2))
	}
	matrix <- matrix[rownames(matrix) %in% gene_names, , drop = FALSE]	
	


	if (dim(matrix)[1] <= 1) {
		# There are no genes to print.
		write_tsv(data.frame("Only 1 or fewer genes. Too few to plot as heatmap."), paste0(filename_prefix, ".txt"), col_names = FALSE)
	    continue <- FALSE
	} else if (dim(matrix)[1] > 1 & dim(matrix)[1] <= 399) {
		# Do no filtering based on number of genes
		print("heatmap matrix has less than 400 genes, continue with original genes.")
		continue <- TRUE
	} else {
		# Gene filtering (only if the number of genes is too large)
		# Use the genefilter function to find most variant genes
		# This function requires an ExpressionSet as input
		# It also requires a valid annotationDBI database for merging entrezIDs
		max_number_of_genes_allowed <- 400
		curr_total_number_of_genes <- dim(matrix)[1]
		if (curr_total_number_of_genes > max_number_of_genes_allowed) {
			# Subset the matrix to keep only max 400 genes
			eset <- ExpressionSet(assayData = matrix)
			fractional_cutoff <- 1 - (max_number_of_genes_allowed / curr_total_number_of_genes)
			fractional_cutoff <- round(fractional_cutoff, 1)
			# Cannot round to zero or one, force it to be a small or large fraction
			if (fractional_cutoff == 0) {
			    fractional_cutoff <- 0.001
			} else if (fractional_cutoff == 1) {
			    fractional_cutoff <- 0.95
			}
			# keep the most variant genes across all samples
			var <- genefilter::varFilter(eset, var.func = IQR, var.cutoff = fractional_cutoff, filterByQuantile = TRUE)
			# Get a logical, specifying which of the original genes to keep
			keep <- rownames(matrix) %in% rownames(var@assayData$exprs)
			# Subset the original dataset
			matrix <- matrix[keep, ]
		}
		continue <- TRUE
	}	
	
	if (continue == TRUE) {
        # Cell filtering
        # If there are too many cells to display well, reduce number of cells randomly
        max_number_of_cells_allowed <- 1500
        curr_total_number_of_cells <- dim(matrix)[2]
        if (curr_total_number_of_cells > max_number_of_cells_allowed) {
            # Randomly remove cells
            print("heatmap matrix has more than 1500 cells, randomly downsampling cells number to 1500 total.")
            set.seed(20)
            keep_index <- sample(seq_along(colnames(matrix)), size = max_number_of_cells_allowed, replace = FALSE)
            matrix <- matrix[, keep_index]
            active.ident <- active.ident[keep_index]
        }
        # Order the cells by cluster 
        cells_order <- order(active.ident, decreasing = FALSE)
        matrix <- matrix[, cells_order]
        # Create df for cluster assignments
        annotation <- data.frame(Cluster = factor(active.ident)[cells_order])
        # Get a bunch of discrete colors
        gg_color_hue <- function(n, alpha = 1) {
          hues <- seq(15, 375, length = n + 1)
          hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
        }
        my_heatmap_colors <- gg_color_hue(length(levels(annotation$Cluster)), alpha = 1)
        annotation_colors_col <- list(Cluster = my_heatmap_colors)
    
    
        # Finally, check to see if there is zero standard deviation for each gene, across all cells.
        # If there is, the clustering will fail.
        # Which rows (genes) do NOT have sd=0? Keep those rows.
        if (any(apply(matrix, 1, sd) == 0)) {
            drop_genes <- rownames(matrix)[which(apply(matrix, 1, sd) == 0)]
            warning(paste("Genes that have zero standard deviation across cells will be removed from the heatmap:", paste(drop_genes, collapse = " ")))
            write_tsv(data.frame(paste("Genes that have zero standard deviation across cells will be removed from the heatmap:", paste(drop_genes, collapse = " "))), paste0(filename_prefix, "_genes_removed.txt"), col_names = FALSE)
            # Keep all other genes
            matrix <- matrix[which(!apply(matrix, 1, sd) == 0), ]
        }
    
        # ward clustering different between versions 0.21.0 and 0.20.6. The newer version is correct.
        # https://github.com/renozao/NMF/issues/117
        # Thus, do the ward1 clustering manually, so it's correct, then feed that to aheatmap.
        # https://stats.stackexchange.com/questions/109949/what-algorithm-does-ward-d-in-hclust-implement-if-it-is-not-wards-criterion
        # http://girke.bioinformatics.ucr.edu/GEN242/pages/mydoc/Rclustering.html
        #matrix_hclust_cols <- hclust((as.dist(1-cor(matrix, method="pearson"))^2), method="ward.D")
        matrix_hclust_rows <- hclust((as.dist(1-cor(t(matrix), method="pearson"))^2), method="ward.D")

        
        # Plot heatmap
        if (length(rownames(matrix)) < 15 ) {
            pdf(paste0(filename_prefix, ".pdf"), width = 6, height = 3, onefile = FALSE)
        } else if (length(rownames(matrix)) >= 15 & length(rownames(matrix)) < 50) {
            pdf(paste0(filename_prefix, ".pdf"), width = 8, height = 4, onefile = FALSE)
        } else if (length(rownames(matrix)) >= 50 & length(rownames(matrix)) < 100) {
            pdf(paste0(filename_prefix, ".pdf"), width = 10, height = 6, onefile = FALSE)
        } else {
            pdf(paste0(filename_prefix, ".pdf"), width = 12, height = 8, onefile = FALSE)
        }
        map <- aheatmap(matrix,
                        labRow = rownames(matrix),
                        labCol = NA,
                        scale = "none",
                        cexRow = 0.75,
                        Colv = NA,
                        Rowv = matrix_hclust_rows,
                        annCol = annotation,
                        annColors = annotation_colors_col,
                        #distfun = "pearson",
                        #hclustfun = "ward", 
                        #Colv = NA,
                        #Rowv = TRUE,
                        color = colorRampPalette(c("forestgreen", "forestgreen", "forestgreen", "royalblue1", "white", "firebrick1", "yellow", "yellow", "yellow"))(300),
                        breaks = seq(-10, 10, length.out = 301),
                        main = heatmap_title
        )
        dev.off()
        write_tsv(data.frame(heatmap_rownames = rev(rownames(matrix)[map$rowInd])), paste0(filename_prefix, "_rownames.txt"), col_names = TRUE)
        # colnames are not clustered, so use original cells_order
        write_tsv(cbind(data.frame(heatmap_colnames = colnames(matrix)), annotation), paste0(filename_prefix, "_colnames_annot.txt"), col_names = TRUE)
        matrix_out <- matrix %>% as.data.frame() %>% rownames_to_column(var = "symbol") %>% as_tibble()
        write_tsv(matrix_out, paste0(filename_prefix, "_matrix.txt"), col_names = TRUE)
    }
}








# ---------------------------------------------------------------------
# Gene Plots
# ---------------------------------------------------------------------

#' @export
tk_gene_plots <- function(seurat_object, genes_list, genes_list_name, all_together = FALSE, number_of_ind_plots = 16, umap_ggplot_object = NULL, type = NULL, facet_by = NULL) {
	# Get only genes present in the dataset
	genes_list <- intersect(genes_list, rownames(GetAssayData(seurat_object, slot = "data")))
	out_dir <- paste0(genes_list_name, "/")
	if (length(genes_list) == 0) {
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		write_tsv(data.frame("No genes to plot"), paste0(genes_list_name, ".txt"), col_names = FALSE)
	} else {
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		# Print all genes individually
		if (length(genes_list) >= number_of_ind_plots) {
			genes_list_short <- genes_list[1:number_of_ind_plots]
		} else {
			genes_list_short <- genes_list
		}		
		ind_plots_list <- list()
		for (i in seq_along(genes_list_short)) {
			if (type == "feature_plot") {
			    pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 7, height = 7)
			    print(
			    ind_plots_list[[i]] <- tk_feature_plot(object = seurat_object, slot = "data", feature = genes_list_short[i], subtitle = NULL, facet_by = facet_by)
			    )
			    dev.off()
			} else if (type == "violin_plot") {
			    pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 7, height = 7)
			    print(
			    ind_plots_list[[i]] <- tk_violin_feature_plot(object = seurat_object, slot = "data", feature = genes_list_short[i], subtitle = NULL, facet_by = facet_by)
			    )
			    dev.off()
			} else if (type == "box_plot") {
			    pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 7, height = 7)
			    print(
			    ind_plots_list[[i]] <- tk_box_feature_plot(object = seurat_object, slot = "data", feature = genes_list_short[i], subtitle = NULL, facet_by = facet_by)
			    )
			    dev.off()
			} else if (type == "ridge_plot") {
			    pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 7, height = 7)
			    print(
			    ind_plots_list[[i]] <- tk_ridge_feature_plot(object = seurat_object, slot = "data", feature = genes_list_short[i], subtitle = NULL, facet_by = facet_by)
			    )
			    dev.off()
			} else if (type == "dot_plot") {
			    pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 4, height = 7)
			    print(
			    ind_plots_list[[i]] <- tk_dot_feature_plot(object = seurat_object, slot = "data", feature = genes_list_short[i], subtitle = NULL, facet_by = facet_by)
			    )
			    dev.off()
			} else {
			    stop("no plot type")
			}
		}
		if (all_together) {
		    if (type == "dot_plot") {
		        # Do not patchwork together previous individual plots. Do not add umap_ggplot_object
                seurat_object@active.ident <- factor(seurat_object@active.ident, levels = gtools::mixedsort(levels(seurat_object@active.ident)), labels = gtools::mixedsort(levels(seurat_object@active.ident)))
                # NOTE: regardless of data slot, the dataset is "scale()" across clusters
                # https://github.com/satijalab/seurat/blob/9843b843ed0c3429d86203011fda00badeb29c2e/R/visualization.R#L3464
                # The size of the dots ("percent expressed") will vary depending on which genes are included in the plot
                if (!is.null(facet_by)) {
                    plot <- DotPlot(seurat_object, features = genes_list_short, dot.scale = 8, group.by = facet_by) + RotatedAxis() + 
                        scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
                        labs(title = "", x = "Gene", y = "Cluster", fill = "Cluster") +
                            theme(plot.title = element_text(hjust = 0.5))
                    pdf(paste0(genes_list_name, ".pdf"), width = 6, height = 9)
                    print(plot)
                    dev.off()
                } else {
                    plot <- DotPlot(seurat_object, features = genes_list_short, dot.scale = 8) + RotatedAxis() + 
                        scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
                        labs(title = "", x = "Gene", y = "Cluster", fill = "Cluster") +
                            theme(plot.title = element_text(hjust = 0.5))
                    pdf(paste0(genes_list_name, ".pdf"), width = 6, height = 9)
                    print(plot)
                    dev.off()                
                }
		    } else {
                # Print all genes together in one plot, but limit to 16 genes
                if ("gg" %in% class(umap_ggplot_object)) {
                    # Include the UMAP cluster plot with the expression plots
                    if (length(ind_plots_list) >= 15) {
                        ind_plots_list_short <- ind_plots_list[c(1:15)]
                    } else {
                        ind_plots_list_short <- ind_plots_list
                    }
                    # Append the lists
                    umap_and_exprs_plotlist <- list()
                    umap_and_exprs_plotlist[[1]] <- umap_ggplot_object
                    for (p in seq_len(length(ind_plots_list_short))) {
                        umap_and_exprs_plotlist[[p + 1]] <- ind_plots_list_short[[p]]
                    }
                    # Figure out best pdf dimensions
                    n_plots <- length(umap_and_exprs_plotlist)
                    if (n_plots <= 4) {
                        pdf(paste0(genes_list_name, ".pdf"), width = 20, height = 4.5)
                        print(
                        patchwork::wrap_plots(umap_and_exprs_plotlist, guides = "keep", widths = 1, heights = 1, nrow = 1, ncol = 4)
                        )
                        dev.off()
                    } else if (n_plots > 4 & n_plots <= 8) {
                        pdf(paste0(genes_list_name, ".pdf"), width = 20, height = 9)
                        print(
                        patchwork::wrap_plots(umap_and_exprs_plotlist, guides = "keep", widths = 1, heights = 1, nrow = 2, ncol = 4)
                        )
                        dev.off()
                    } else if (n_plots > 8 & n_plots <= 12) {
                        pdf(paste0(genes_list_name, ".pdf"), width = 20, height = 13.5)
                        print(
                        patchwork::wrap_plots(umap_and_exprs_plotlist, guides = "keep", widths = 1, heights = 1, nrow = 3, ncol = 4)
                        )
                        dev.off()
                    } else {
                        pdf(paste0(genes_list_name, ".pdf"), width = 20, height = 18)
                        print(
                        patchwork::wrap_plots(umap_and_exprs_plotlist, guides = "keep", widths = 1, heights = 1, nrow = 4, ncol = 4)
                        )
                        dev.off()
                    }
                } else {
                    # Do not include the UMAP cluster plot with the expression plots
                    if (length(ind_plots_list) >= 16) {
                        ind_plots_list_short <- ind_plots_list[c(1:16)]
                    } else {
                        ind_plots_list_short <- ind_plots_list
                    }
                    umap_and_exprs_plotlist <- ind_plots_list_short
                    pdf(paste0(genes_list_name, ".pdf"), width = 14, height = 14)
                    print(
                    patchwork::wrap_plots(umap_and_exprs_plotlist, guides = "keep", widths = 1, heights = 1)
                    )
                    dev.off()
                }
			}
		}
	}
}















#' @export
tk_merge_clusters_by_de <- function(seurat_object, cluster_res_name = cluster_res_name[j], lfc_threshold, num_genes_diff_between_clusters_threshold) {
	iteration_counter <- 0
	continue_merging_clusters <- TRUE
	# List to store all DE comparison data
	#merge_data_list <- vector("list", 1000000)
	merge_data_list <- list()
	iteration_data_list <- vector("list", 1000000)
	# While the variable above is true, continue to merge clusters, until unique clusters are found (based on DE analysis)
	while (continue_merging_clusters == TRUE) {
		iteration_counter <- iteration_counter + 1
		print(cluster_res_name)
		print(paste0("iteration: ", iteration_counter))
		curr_clusters <- gtools::mixedsort(levels(seurat_object@active.ident))
		curr_number_of_clusters <- length(curr_clusters)
		# Generate pairwise cluster numbering for all cluster comparisons
		combos <- combn(curr_clusters, 2)
		# Run all pairwise cluster comparisons, find DE genes
		k_vect <- 1:dim(combos)[2]
		# Help with parallel lapply:
		# https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
		out <- mclapply(X = k_vect, FUN = tk_parallel_de_pairwise, seurat_object = seurat_object, curr_clusters = curr_clusters, curr_number_of_clusters = curr_number_of_clusters, combos = combos, lfc_threshold = lfc_threshold, merge_data_list = merge_data_list,
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		# name the list elements
		names(out) <- sapply(out, function(x) x[["name"]])
		number_of_sig_genes_for_each_comparison <- unname(sapply(out, function(x) x[["number_of_signif_genes"]]))
		# Are there any pairwise comparisons that have less than "num_genes_diff_between_clusters_threshold" DE genes between them? If so, merge clusters.
		combos2 <- as.data.frame(t(combos))
		combos2 <- data.frame(combos2, number_of_sig_genes_for_each_comparison)
		colnames(combos2) <- c("comparison_a", "comparison_b", "number_of_sig_genes_for_each_comparison")
		# Update the merge_data_list
		# Find length of prev merge data list
		merge_data_list_length <- length(merge_data_list)
		# append lists
		for (i in seq_along(out)) {
			idx <- merge_data_list_length + i
			merge_data_list[[idx]] <- out[[i]]
			names(merge_data_list)[idx] <- out[[i]]$name
		}
		# Get rid of any duplicate entries in the appended merge_data_list
		merge_data_list <- merge_data_list[!duplicated(names(merge_data_list))]
		if (any(number_of_sig_genes_for_each_comparison < num_genes_diff_between_clusters_threshold)) {
			print(paste0("There are cluster comparisons with less than ", num_genes_diff_between_clusters_threshold, " DE genes between clusters -- starting to merge them."))
			# Which clusters to merge?
			# Which comparisons have too few DE genes, and will be merged?
			# Calculate the cluster medoids in PCA space, only for the cells that have "too few" DE genes between them
			if (DefaultAssay(seurat_object) == "RNA" & "pca_rna" %in% names(seurat_object@reductions)) {
				medoids <- tk_get_medoids(m = t(seurat_object@reductions$pca_rna@cell.embeddings), cluster = seurat_object@active.ident)
			} else {
				medoids <- tk_get_medoids(m = t(seurat_object@reductions$pca@cell.embeddings), cluster = seurat_object@active.ident)
			}
			medoids_dist <- as.matrix(dist(t(medoids), diag = TRUE, upper = FALSE))
			medoids_dist[upper.tri(medoids_dist, diag = TRUE)] <- NA
			medoids_dist_df <- as.data.frame(medoids_dist)
			medoids_dist_df <- rownames_to_column(medoids_dist_df, var = "cluster")
			medoids_dist_df <- gather(medoids_dist_df, cluster2, value, -cluster)
			medoids_dist_df <- medoids_dist_df[!is.na(medoids_dist_df$value), ]
			medoids_dist_df$new_order1 <- NA
			medoids_dist_df$new_order2 <- NA
			for (v in 1:dim(medoids_dist_df)[1]) {
				new_order <- gtools::mixedsort(levels(factor(c(medoids_dist_df$cluster[v], medoids_dist_df$cluster2[v]))))
				medoids_dist_df$new_order1[v] <- new_order[1]
				medoids_dist_df$new_order2[v] <- new_order[2]
			}
			medoids_dist_df <- medoids_dist_df[tk_multi_mixedorder(medoids_dist_df$new_order1, medoids_dist_df$new_order2), ]
			merged_medoids <- merge(combos2, medoids_dist_df, by.x = c("comparison_a", "comparison_b"), by.y = c("new_order1", "new_order2"))
			# Subset list to include only above clusters available for merging.
			merged_medoids <- merged_medoids[merged_medoids$number_of_sig_genes_for_each_comparison < num_genes_diff_between_clusters_threshold, ]
			merged_medoids <- merged_medoids[order(merged_medoids$value, decreasing = FALSE), ]
			# Get the top cluster comparison for merging
			merge_small <- as.character(merged_medoids$comparison_a[1])
			merge_large <- as.character(merged_medoids$comparison_b[1])
			# Re-create cell new merged classes
			curr_classes <- seurat_object@active.ident	
			new_classes <- as.character(curr_classes)
			new_classes[curr_classes == merge_small] <- paste0(merge_small, ".", merge_large)
			new_classes[curr_classes == merge_large] <- paste0(merge_small, ".", merge_large)
			names(new_classes) <- names(curr_classes)
			new_classes <- factor(new_classes)
			new_classes <- droplevels(new_classes)
			# Re-assign the class names to @ident slot
			seurat_object@active.ident <- new_classes
			# Repeat pairwise again
			continue_merging_clusters <- TRUE
			clusters_merged <- c(merge_small, merge_large)
		} else {
			curr_classes <- seurat_object@active.ident
			final_classes <- curr_classes
			print(paste0("All clusters have more than ", num_genes_diff_between_clusters_threshold, " DE genes between them. No merging."))
			continue_merging_clusters <- FALSE
			clusters_merged <- "none"
		}
		names(iteration_data_list)[iteration_counter] <- paste0("iteration_", iteration_counter)
		iteration_data_list[[iteration_counter]] <- list(pre_iteration_classes = curr_classes, post_iteration_classes = seurat_object@active.ident, pairwise_signif_genes = combos2, were_clusters_merged = continue_merging_clusters, clusters_merged = clusters_merged)
	}
	results_list <- list(final_classes = final_classes, merge_data_list = merge_data_list, iteration_data_list = iteration_data_list)
	return(results_list)
}






#' @export
tk_int_parallel_de_pairwise <- function(k, seurat_object, curr_clusters, curr_number_of_clusters, combos, lfc_threshold, de_data_list, meta_data_col) {
	# Start running function
	curr_cluster_a <- combos[1, k]
	curr_cluster_b <- combos[2, k]
	curr_cluster_vect <- c(curr_cluster_a, curr_cluster_b)
	name <- paste0("cluster_", curr_cluster_a, "_vs_", curr_cluster_b, "_conserved_across_samples")
	# If the DE comparison was done previously, and is not changed here, don't repeat.
	if (any(names(de_data_list) %in% name)) {
		# Get previous results stored in "de_data_list"
		prev_de_data_list_idx <- which(names(de_data_list) == name)
		print(paste0(names(de_data_list)[prev_de_data_list_idx], " comparison already done, skipping."))
		results <- de_data_list[[prev_de_data_list_idx]]
	} else {
		# Initialize to "not" skip this comparison
		curr_de_skip <- FALSE
		unique_orig_ident <- unique(seurat_object@meta.data$orig.ident)
		for (v in seq_along(unique_orig_ident)) {
			for (w in seq_along(curr_cluster_vect)) {
				sample_specific_cells <- seurat_object@meta.data$orig.ident == unique_orig_ident[v] & seurat_object@meta.data[[meta_data_col]] == curr_cluster_vect[w]
				if (dim(GetAssayData(seurat_object, slot = "data")[, sample_specific_cells, drop = FALSE])[2] < 3) {
					# Too few files for DE analysis, skip DE analysis
					print(paste0("Too few cells for DE analysis, skipping DE analysis for orig.ident: ", unique_orig_ident[v], ", cluster: ", curr_cluster_vect[w]))
					curr_de_skip <- TRUE
					break
				}
			if (curr_de_skip == TRUE) {break}	
			}	
		}
		if (curr_de_skip == TRUE) {
			colnames_list <- list()
			for (v in seq_along(unique_orig_ident)) {
				colnames_list[[v]] <- c(
				paste0(unique_orig_ident[v], "_p_val"), 
				paste0(unique_orig_ident[v], "_avg_log2FC"), 
				paste0(unique_orig_ident[v], "_pct.1"), 
				paste0(unique_orig_ident[v], "_pct.2"), 
				paste0(unique_orig_ident[v], "_p_val_adj"))
			}
			colnames_vect <- c("symbol", unlist(list(colnames_list)), "max_pval", "minimump_p_val")
			# Create empty tibble
			de_blank <- matrix(0, 0, ncol = length(colnames_vect))
			colnames(de_blank) <- colnames_vect
			de_blank <- as_tibble(de_blank)
			de_blank <- de_blank %>%
  				mutate_all(as.character)
			de <- de_blank
			de_signif <- de_blank
			number_of_signif_genes <- 0
			results <- list(name = name, cluster_a = curr_cluster_a, cluster_b = curr_cluster_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
		} else {
			print(paste0(name, " running DE analysis."))
			# Find the DE genes for a new comparison
			# Run DE analysis
			de <- FindConservedMarkers(object = seurat_object, grouping.var = "orig.ident", slot = "data", ident.1 = combos[1, k], ident.2 = combos[2, k], min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.1)
			de <- de %>%
				tibble::rownames_to_column(var = "symbol") %>%
				as_tibble()
			# Extract significant DE genes
			# If multiple groups were integrated, then each "group" will have a set of columns
			lfc_col <- colnames(de)[(grepl("log2FC", colnames(de)))]
			# This should only be one column
			adj_p_col <- colnames(de)[(grepl("max_pval", colnames(de)))]
		
			# Pass thresholds
			# Make sure all samples have same the direction of log fold change (both pos, or both neg)
			# Sort the table by the first log2FC column group
			de_signif <- de %>%
				dplyr::filter(apply(abs(.[, lfc_col]) >= lfc_threshold, 1, all) & .[, adj_p_col] <= 0.01) %>%
				dplyr::filter(apply(sign(.[, lfc_col]), 1, function(x) length(unique(as.numeric(x))) == 1)) %>%
				dplyr::arrange(desc(!!sym(lfc_col[1])))
			# How many genes are significant?
			number_of_signif_genes <- dim(de_signif)[1]
			results <- list(name = name, cluster_a = curr_cluster_a, cluster_b = curr_cluster_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
		}
	}
	return(results)
}




#' @export
tk_int_parallel_de_one_vs_all <- function(m, seurat_object, lfc_threshold, meta_data_col) {
	curr_cluster <- levels(factor(seurat_object@active.ident))[m]
	curr_cluster_a <- paste0("cluster_", curr_cluster)
	name <- paste0(curr_cluster_a, "_vs_all_conserved_across_samples")
	# Initialize to "not" skip 
	curr_de_skip <- FALSE
	unique_orig_ident <- unique(seurat_object@meta.data$orig.ident)
    print(paste0(name, " DE analysis."))
    results <- tryCatch(
        expr = {
            de <- FindConservedMarkers(object = seurat_object, grouping.var = "orig.ident", slot = "data", ident.1 = curr_cluster, min.pct = 0.1, test.use = "wilcox", logfc.threshold = 0.1) %>%
                tibble::rownames_to_column(var = "symbol") %>%
                as_tibble()
            # Extract significant DE genes
            # If multiple groups were integrated, then each "group" will have a set of columns
            lfc_col <- colnames(de)[(grepl("log2FC", colnames(de)))]
            # This should only be one column
            adj_p_col <- colnames(de)[(grepl("max_pval", colnames(de)))]
            # Pass thresholds
            # Make sure all samples have same the direction of log fold change (both pos, or both neg)
            # Sort the table by the first log2FC column group
            de_signif <- de %>%
                dplyr::filter(apply(abs(.[, lfc_col]) >= lfc_threshold, 1, all) & .[, adj_p_col] <= 0.01) %>%
                dplyr::filter(apply(sign(.[, lfc_col]), 1, function(x) length(unique(as.numeric(x))) == 1)) %>%
                dplyr::arrange(desc(!!sym(lfc_col[1])))
            # How many genes are significant?
            number_of_signif_genes <- dim(de_signif)[1]
            list(name = name, cluster_a = curr_cluster_a, cluster_b = "all_other_cells", de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
        },
        error = function(e) {
            # (Optional)
            # Do this if an error is caught...
            colnames_list <- list()
            for (v in seq_along(unique_orig_ident)) {
                colnames_list[[v]] <- c(
                paste0(unique_orig_ident[v], "_p_val"), 
                paste0(unique_orig_ident[v], "_avg_log2FC"), 
                paste0(unique_orig_ident[v], "_pct.1"), 
                paste0(unique_orig_ident[v], "_pct.2"), 
                paste0(unique_orig_ident[v], "_p_val_adj"))
            }
            colnames_vect <- c("symbol", unlist(list(colnames_list)), "max_pval", "minimump_p_val")
            # Create empty tibble
            de_blank <- matrix(0, 0, ncol = length(colnames_vect))
            colnames(de_blank) <- colnames_vect
            de_blank <- as_tibble(de_blank)
            de_blank <- de_blank %>%
                mutate_all(as.character)
            de <- de_blank
    		de_signif <- de_blank
    		number_of_signif_genes <- 0
    		list(name = name, cluster_a = curr_cluster_a, cluster_b = "all_other_cells", de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
        },
        warning = function(w){
            # (Optional)
            # Do this if an warning is caught...
            colnames_list <- list()
            for (v in seq_along(unique_orig_ident)) {
                colnames_list[[v]] <- c(
                paste0(unique_orig_ident[v], "_p_val"), 
                paste0(unique_orig_ident[v], "_avg_log2FC"), 
                paste0(unique_orig_ident[v], "_pct.1"), 
                paste0(unique_orig_ident[v], "_pct.2"), 
                paste0(unique_orig_ident[v], "_p_val_adj"))
            }
            colnames_vect <- c("symbol", unlist(list(colnames_list)), "max_pval", "minimump_p_val")
            # Create empty tibble
            de_blank <- matrix(0, 0, ncol = length(colnames_vect))
            colnames(de_blank) <- colnames_vect
            de_blank <- as_tibble(de_blank)
            de_blank <- de_blank %>%
                mutate_all(as.character)
            de <- de_blank
    		de_signif <- de_blank
    		number_of_signif_genes <- 0
    		list(name = name, cluster_a = curr_cluster_a, cluster_b = "all_other_cells", de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)

        },
        finally = {
            # (Optional)
            # Do this at the end before quitting the tryCatch structure...
        }
    )
    return(results)
}	



#' @export
tk_int_parallel_de_pairwise_conserved_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	# set3: DE LISTS PAIRWISE
	if (!dir.exists("de_pairwise_cons/umap_gene_expression")) {dir.create("de_pairwise_cons/umap_gene_expression", recursive = TRUE)}
	if (!dir.exists("de_pairwise_cons/de_text_files")) {dir.create("de_pairwise_cons/de_text_files", recursive = TRUE)}
	write_tsv(de_data_list[[r]]$de, paste0("de_pairwise_cons/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_pairwise_cons/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS PAIRWISE
	starting_wd <- getwd()
	setwd("de_pairwise_cons/umap_gene_expression")
	# "de_signif" lists only contain log2FC values that match between each sample. Therefore, pick either sample to find up/dn.
	first_lfc_col <- which(str_detect(colnames(de_data_list[[r]]$de_signif), fixed("avg_log2FC")))[1]
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif[[first_lfc_col]] > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif[[first_lfc_col]] < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(starting_wd)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}




#' @export
tk_int_parallel_de_one_vs_all_conserved_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	# DE LISTS ALL_VS_ONE
	if (!dir.exists("de_one_vs_all_cons/de_text_files")) {dir.create("de_one_vs_all_cons/de_text_files", recursive = TRUE)}
	if (!dir.exists("de_one_vs_all_cons/umap_gene_expression")) {dir.create("de_one_vs_all_cons/umap_gene_expression", recursive = TRUE)}
	write_tsv(de_data_list[[r]]$de, paste0("de_one_vs_all_cons/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_one_vs_all_cons/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS ALL_VS_ONE
	starting_dir <- getwd()
	setwd("de_one_vs_all_cons/umap_gene_expression")
	# "de_signif" lists only contain log2FC values that match between each sample. Therefore, pick either sample to find up/dn.
	first_lfc_col <- which(str_detect(colnames(de_data_list[[r]]$de_signif), fixed("avg_log2FC")))[1]
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif[[first_lfc_col]] > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif[[first_lfc_col]] < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(starting_dir)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}



#' @export
tk_int_dotplots <- function(seurat_object, results_list) {
	# DOTPLOTS	
	# For each cluster, get 2 genes of interest 
	dotplot_genes <- vector()
	for (i in seq_along(results_list)) {
		geneset_up <- results_list[[i]]$geneset_up
		if (length(geneset_up) > 0) {
			# Get top 2 genes or less
			num_genes_to_keep <- min(2, length(geneset_up))
			dotplot_genes <- c(dotplot_genes, geneset_up[1:num_genes_to_keep])
		}
	}
	dotplot_genes <- unique(dotplot_genes)
	if (length(dotplot_genes) > 0) {
		# I think this is pulling the integrated assay (or whatever is default), scale.data slot
		# https://github.com/satijalab/seurat/issues/1620#issuecomment-504177635
		seurat_object_temp1 <- seurat_object
		new_ident <- factor(str_replace_all(seurat_object_temp1@active.ident, "_", "."))
		names(new_ident) <- rownames(seurat_object_temp1@meta.data)
		seurat_object_temp1@active.ident <- new_ident
		seurat_object_temp1@meta.data$sample <- str_replace_all(seurat_object_temp1@meta.data$orig.ident, "_", ".")
		num_of_samples <- length(levels(factor(seurat_object_temp1@meta.data$sample)))
		dotplot_sample_colors <- c("#e5323d",  "#4c43dc", "#46d53c", "#cac027", "#db48e9", "#7847c2", "#007100", "#4622a2", "#b4a1f3", "#b87000", "#a52fb1")[seq_len(num_of_samples)]
		pdf("dotplot_genes_split_by_sample_up.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8, split.by = "sample") + RotatedAxis())
		dev.off()
		pdf("dotplot_genes_all_samples_up.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8) + RotatedAxis())
		dev.off()
		rm(seurat_object_temp1)
	} else {
		write_tsv(data.frame("No genes to plot"), "dotplot_genes_all_samples_up.txt", col_names = FALSE)
	}
	# For each cluster, get 2 genes of interest 
	dotplot_genes <- vector()
	for (i in seq_along(results_list)) {
		geneset_dn <- results_list[[i]]$geneset_dn
		if (length(geneset_dn) > 0) {
			# Get top 2 genes or less
			num_genes_to_keep <- min(2, length(geneset_dn))
			dotplot_genes <- c(dotplot_genes, geneset_dn[1:num_genes_to_keep])
		} 
	}
	dotplot_genes <- unique(dotplot_genes)
	if (length(dotplot_genes) > 0 ) {
		# I think this is pulling the integrated assay (or whatever is default), scale.data slot
		# https://github.com/satijalab/seurat/issues/1620#issuecomment-504177635
		seurat_object_temp1 <- seurat_object
		new_ident <- factor(str_replace_all(seurat_object_temp1@active.ident, "_", "."))
		names(new_ident) <- rownames(seurat_object_temp1@meta.data)
		seurat_object_temp1@active.ident <- new_ident
		seurat_object_temp1@meta.data$sample <- str_replace_all(seurat_object_temp1@meta.data$orig.ident, "_", ".")
		num_of_samples <- length(levels(factor(seurat_object_temp1@meta.data$sample)))
		dotplot_sample_colors <- c("#e5323d",  "#4c43dc", "#46d53c", "#cac027", "#db48e9", "#7847c2", "#007100", "#4622a2", "#b4a1f3", "#b87000", "#a52fb1")[seq_len(num_of_samples)]
		pdf("dotplot_genes_split_by_sample_dn.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8, split.by = "sample") + RotatedAxis())
		dev.off()
		pdf("dotplot_genes_all_samples_dn.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8) + RotatedAxis())
		dev.off()
		rm(seurat_object_temp1)
	} else {
		write_tsv(data.frame("No genes to plot"), "dotplot_genes_all_samples_dn.txt", col_names = FALSE)
	}
	# For each cluster, get 2 genes of interest 
	dotplot_genes <- vector()
	for (i in seq_along(results_list)) {
		geneset_up <- results_list[[i]]$geneset_up
		geneset_dn <- results_list[[i]]$geneset_dn
		if (length(geneset_up) > 0) {
			keep_up <- geneset_up[1]
			dotplot_genes <- c(dotplot_genes, keep_up)
		}
		if (length(geneset_dn) > 0) {
			keep_dn <- geneset_dn[1]
			dotplot_genes <- c(dotplot_genes, keep_dn)
		}
	}
	dotplot_genes <- unique(dotplot_genes)
	if (length(dotplot_genes) > 0 ) {
		# I think this is pulling the integrated assay (or whatever is default), scale.data slot
		# https://github.com/satijalab/seurat/issues/1620#issuecomment-504177635
		seurat_object_temp1 <- seurat_object
		new_ident <- factor(str_replace_all(seurat_object_temp1@active.ident, "_", "."))
		names(new_ident) <- rownames(seurat_object_temp1@meta.data)
		seurat_object_temp1@active.ident <- new_ident
		seurat_object_temp1@meta.data$sample <- str_replace_all(seurat_object_temp1@meta.data$orig.ident, "_", ".")
		num_of_samples <- length(levels(factor(seurat_object_temp1@meta.data$sample)))
		dotplot_sample_colors <- c("#e5323d",  "#4c43dc", "#46d53c", "#cac027", "#db48e9", "#7847c2", "#007100", "#4622a2", "#b4a1f3", "#b87000", "#a52fb1")[seq_len(num_of_samples)]
		pdf("dotplot_genes_split_by_sample_updn.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8, split.by = "sample") + RotatedAxis())
		dev.off()
		pdf("dotplot_genes_all_samples_updn.pdf", width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = dotplot_sample_colors, dot.scale = 8) + RotatedAxis())
		dev.off()
		rm(seurat_object_temp1)
	} else {
		write_tsv(data.frame("No genes to plot"), "dotplot_genes_all_samples_updn.txt", col_names = FALSE)
	}		
}



#' @export
tk_int_gene_plots <- function(seurat_object, genes_list, genes_list_name, all_together = FALSE, number_of_ind_plots = 8, umap_ggplot_object = NULL) {
	# Get only genes present in the dataset
	genes_list <- intersect(genes_list, rownames(GetAssayData(seurat_object, slot = "data")))
	out_dir <- paste0(genes_list_name, "/")
	if (length(genes_list) == 0) {
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		write_tsv(data.frame("No genes to plot"), paste0(genes_list_name, ".txt"), col_names = FALSE)
	} else {
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		# Print all genes individually
		if (length(genes_list) >= number_of_ind_plots) {
			genes_list_short <- genes_list[1:number_of_ind_plots]
		} else {
			genes_list_short <- genes_list
		}		
		ind_plots_list <- list()
		for (i in seq_along(genes_list_short)) {
			pdf(paste0(out_dir, genes_list_short[i], ".pdf"), width = 12, height = 7)
			print(
			ind_plots_list[[i]] <- FeaturePlot(object = seurat_object, slot = "data", features = genes_list_short[i], order = TRUE, cols = colorRampPalette(c("white", "red"))(100), reduction = "umap", combine = TRUE, split.by = "orig.ident")
			)
			dev.off()
		}
		if (all_together) {
			# Print all genes together in one plot, but limit to 7 genes
			if ("gg" %in% class(umap_ggplot_object)) {
				# Include the UMAP cluster plot with the expression plots
				if (length(ind_plots_list) >= 7) {
					ind_plots_list_short <- ind_plots_list[c(1:7)]
				} else {
					ind_plots_list_short <- ind_plots_list
				}
				# Create blank plot
				blank_plot <- ggplot()+geom_blank(aes(1,1))+
				theme(plot.background = element_blank(), 
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_blank(), 
				axis.text.y = element_blank(),
				axis.ticks = element_blank(),
				axis.line = element_blank())
     			# Append all the plot lists
				umap_and_exprs_plotlist <- list()
				# Create a 2-plot list of the dimplot and a blank plot placeholder
				umap_and_exprs_plotlist[[1]] <- cowplot::plot_grid(umap_ggplot_object, blank_plot, align = "hv", ncol = 2)
				for (p in seq_len(length(ind_plots_list_short))) {
					umap_and_exprs_plotlist[[p + 1]] <- ind_plots_list_short[[p]]
				}
				pdf(paste0(genes_list_name, ".pdf"), width = 20, height = 15)
				print(
				cowplot::plot_grid(plotlist = umap_and_exprs_plotlist, align = "hv", ncol = 2)
				)
				dev.off()
			} else {
				# Include the UMAP cluster plot with the expression plots
				if (length(ind_plots_list) >= 8) {
					ind_plots_list_short <- ind_plots_list[c(1:8)]
				} else {
					ind_plots_list_short <- ind_plots_list
				}
				umap_and_exprs_plotlist <- ind_plots_list_short
				pdf(paste0(genes_list_name, ".pdf"), width = 14, height = 12)
				print(
				cowplot::plot_grid(plotlist = umap_and_exprs_plotlist, align = "hv", ncol = 4)
				)
				dev.off()
			}
		}
	}
}





#' @export
tk_int_parallel_de_pairwise_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	if (!dir.exists("de_pairwise/de_text_files")) {dir.create("de_pairwise/de_text_files", recursive = TRUE)}
	if (!dir.exists("de_pairwise/umap_gene_expression")) {dir.create("de_pairwise/umap_gene_expression", recursive = TRUE)}
	# DE LISTS PAIRWISE
	write_tsv(de_data_list[[r]]$de, paste0("de_pairwise/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_pairwise/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS PAIRWISE
	staring_dir <- getwd()
	setwd("de_pairwise/umap_gene_expression")
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(staring_dir)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}



#' @export
tk_int_parallel_de_one_vs_all_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	if (!dir.exists("de_one_vs_all/de_text_files")) {dir.create("de_one_vs_all/de_text_files", recursive = TRUE)}
	if (!dir.exists("de_one_vs_all/umap_gene_expression")) {dir.create("de_one_vs_all/umap_gene_expression", recursive = TRUE)}
	# DE LISTS ALL_VS_ONE
	write_tsv(de_data_list[[r]]$de, paste0("de_one_vs_all/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_one_vs_all/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS ALL_VS_ONE
	starting_dir <- getwd()
	setwd("de_one_vs_all/umap_gene_expression")
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC > 0]
	tk_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot, type = "feature_plot")
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(starting_dir)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}





#' @export
tk_int_parallel_de_across_samples_plots <- function(r, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot) {
	if (!dir.exists("de_across_samples/de_text_files")) {dir.create("de_across_samples/de_text_files", recursive = TRUE)}
	if (!dir.exists("de_across_samples/umap_gene_expression")) {dir.create("de_across_samples/umap_gene_expression", recursive = TRUE)}
	# DE LISTS ALL_VS_ONE
	write_tsv(de_data_list[[r]]$de, paste0("de_across_samples/de_text_files/", de_data_list[[r]]$name, "_de.txt"))
	write_tsv(de_data_list[[r]]$de_signif, paste0("de_across_samples/de_text_files/", de_data_list[[r]]$name, "_de_signif.txt"))
	# GENE PLOTS ALL_VS_ONE
	starting_dir <- getwd()
	setwd("de_across_samples/umap_gene_expression")
	geneset_up <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC > 0]
	tk_int_gene_plots(seurat_object, geneset_up, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), TRUE, umap_ggplot_object = umap_plot)
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_up"), new_dir = "../umap_gene_expression_flattened")
	geneset_dn <- de_data_list[[r]]$de_signif$symbol[de_data_list[[r]]$de_signif$avg_log2FC < 0]
	# Reverse the order of the downregulated genes, so genes with largest fold change are first
	geneset_dn <- rev(geneset_dn)
	tk_int_gene_plots(seurat_object, geneset_dn, paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), TRUE, umap_ggplot_object = umap_plot)
	#tk_move_and_symlink(original_dir = paste0(de_data_list[[r]]$name, "_de_results_signif_geneplots_dn"), new_dir = "../umap_gene_expression_flattened")
	setwd(starting_dir)
	# Make Excel file
	excel_list_of_df_de <- de_data_list[[r]]$de
	excel_list_of_df_de_signif <- de_data_list[[r]]$de_signif
	signif_genes_across_all_clusters <- de_data_list[[r]]$de_signif$symbol
	results <- list(name = de_data_list[[r]]$name, geneset_up = geneset_up, geneset_dn = geneset_dn, excel_list_of_df_de = excel_list_of_df_de, excel_list_of_df_de_signif = excel_list_of_df_de_signif, signif_genes_across_all_clusters = signif_genes_across_all_clusters)
	return(results)
}







# ---------------------------------------------------------------------
# custom geneset plotting
# ---------------------------------------------------------------------


#' @export
tk_dotplots_geneset_only <- function(seurat_object, geneset, file_name_prefix = "") {
	# DOTPLOTS
	geneset <- intersect(geneset, rownames(GetAssayData(seurat_object, slot = "data")))
	dotplot_genes <- unique(geneset)
	if (length(dotplot_genes) > 0) {
		# I think this is pulling the integrated assay (or whatever is default), scale.data slot
		# https://github.com/satijalab/seurat/issues/1620#issuecomment-504177635
		# Using the default active.ident as provided to this script.
		pdf(glue("{file_name_prefix}dotplot.pdf"), width = 12, height = 9)
		print(DotPlot(seurat_object, features = dotplot_genes, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis())
		dev.off()
	} else {
		write_tsv(data.frame("No genes to plot"), "dotplot.txt", col_names = FALSE)
	}	
}


#' @export
tk_int_dotplots_geneset_only <- function(seurat_object, geneset, file_name_prefix = "") {
	# DOTPLOTS
	geneset <- intersect(geneset, rownames(GetAssayData(seurat_object, slot = "data")))
	dotplot_genes <- unique(geneset)
	if (length(dotplot_genes) > 0) {
		# I think this is pulling the integrated assay (or whatever is default), scale.data slot
		# https://github.com/satijalab/seurat/issues/1620#issuecomment-504177635
		seurat_object_temp1 <- seurat_object
		new_ident <- factor(str_replace_all(seurat_object_temp1@active.ident, "_", "."))
		names(new_ident) <- rownames(seurat_object_temp1@meta.data)
		seurat_object_temp1@active.ident <- new_ident
		seurat_object_temp1@meta.data$sample <- str_replace_all(seurat_object_temp1@meta.data$orig.ident, "_", ".")
		pdf(glue("{file_name_prefix}dotplot_split_by_sample.pdf"), width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = c("blue", "red"), dot.scale = 8, split.by = "sample") + RotatedAxis())
		dev.off()
		pdf(glue("{file_name_prefix}dotplot_all_samples.pdf"), width = 12, height = 9)
		print(DotPlot(seurat_object_temp1, features = dotplot_genes, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis())
		dev.off()
		rm(seurat_object_temp1)
	} else {
		write_tsv(data.frame("No genes to plot"), "dotplot.txt", col_names = FALSE)
	}	
}





#' @export
tk_int_violinplot_geneset_only <- function(seurat_object, geneset, file_name_prefix = "", split_by, group_by) {
	# VIOLIN PLOTS	
	geneset <- intersect(geneset, rownames(GetAssayData(seurat_object, slot = "data")))
	plot_genes <- unique(geneset)
	if (length(plot_genes) > 0) {
		plot_list <- list()
		for (i in seq_along(plot_genes)) {
			# I think this is pulling the integrated assay (or whatever is default), scale.data slot
			plots <- VlnPlot(seurat_object, features = plot_genes[i], cols = c(brewer.pal(n = 8, name = "Set1"), brewer.pal(n = 8, name = "Set1")), split.by = split_by, group.by = group_by, pt.size = 0, combine = FALSE)
			plot_list[[i]] <- CombinePlots(plots = plots, ncol = 1)
		}
		if (length(plot_genes) <= 12) {
			pdf(glue("{file_name_prefix}violinplot_split_by_sample.pdf"), width = 16, height = 12)
			print(
			cowplot::plot_grid(plotlist = plot_list, align = "hv", ncol = 3)
			)
			dev.off()
		} else {
			pdf(glue("{file_name_prefix}violinplot_split_by_sample.pdf"), width = 32, height = 24)
			print(
			cowplot::plot_grid(plotlist = plot_list, align = "hv", ncol = 4)
			)
			dev.off()		
		}
	} else {
		write_tsv(data.frame("No genes to plot"), "violinlot.txt", col_names = FALSE)
	}	
}



#' @export
tk_violinplot_geneset_only <- function(seurat_object, geneset, file_name_prefix = "") {
	# VIOLIN PLOTS	
	# Uses the active.ident for cluster groups
	geneset <- intersect(geneset, rownames(GetAssayData(seurat_object, slot = "data")))
	plot_genes <- unique(geneset)
	if (length(plot_genes) > 0) {
		plot_list <- list()
		for (i in seq_along(plot_genes)) {
			# I think this is pulling the integrated assay (or whatever is default), scale.data slot
			plots <- VlnPlot(seurat_object, features = plot_genes[i], cols = c(brewer.pal(n = 8, name = "Set1"), brewer.pal(n = 8, name = "Set1"), brewer.pal(n = 8, name = "Set1")), pt.size = 0, combine = FALSE)
			plot_list[[i]] <- CombinePlots(plots = plots, ncol = 1)
		}
		if (length(plot_genes) <= 12) {
			pdf(glue("{file_name_prefix}violinplot.pdf"), width = 16, height = 12)
			print(
			cowplot::plot_grid(plotlist = plot_list, align = "hv", ncol = 3)
			)
			dev.off()
		} else {
			pdf(glue("{file_name_prefix}violinplot.pdf"), width = 32, height = 24)
			print(
			cowplot::plot_grid(plotlist = plot_list, align = "hv", ncol = 4)
			)
			dev.off()		
		}
	} else {
		write_tsv(data.frame("No genes to plot"), "violinlot.txt", col_names = FALSE)
	}	
}






























# ---------------------------------------------------------------------
# CITEseq Integrated merge
# ---------------------------------------------------------------------

# CITE-SEQ DEMULTIPLEXED DATA Expression data (not "integrated", but cells are classified by sample)
# Use the "RNA" assay for FeaturePlot(), dotplot, heatmap(scale.data), violinplot, ridgeplot when using the lognorm data
# Use the "SCT" assay for FeaturePlot(), dotplot, heatmap(scale.data), violinplot, ridgeplot when using the SCT data
# Choose either the pca or umap (SCT) or the pca_rna or umap_rna (RNA) slots




#' @export
tk_citeseq_cluster_merge <- function(seurat_object, seurat_object_name, lfc_threshold, num_genes_diff_between_clusters_threshold) {
	print(seurat_object_name)
	# Get the matrix and labels from the clustering results
	cluster_cell_class <- seurat_object@meta.data[, grepl("SCT_snn_res", colnames(seurat_object@meta.data)), drop = FALSE]
	cluster_res_name <- str_replace_all(colnames(cluster_cell_class), fixed("res."), "res_")
	cell_ids <- rownames(cluster_cell_class)
	orig_assay <- DefaultAssay(object = seurat_object)
	orig_wd <- getwd()
	return_list <- list()
	# For each clustering resolution
	for (j in seq_along(cluster_res_name))  {
	    numb_of_samples <- length(levels(factor(seurat_object@meta.data$orig.ident)))
		# Create dir for all results
		out_dir <- cluster_res_name[j]
		if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}"))
		print(cluster_res_name[j])
		# Factor the current clustering resolution cell classes
		class_orig <- factor(cluster_cell_class[, j])
		names(class_orig) <- cell_ids
		# Re-assign the class names to @ident slot with current cell classes
		seurat_object@active.ident <- class_orig
		
		
		# Figure out the right assay to use
		DefaultAssay(seurat_object) <- "SCT"
		print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
		
		
		
		
		print("Starting DE set1, PAIRWISE DE and MERGE CLUSTERS BASED ON SIGNIF DE GENES")
		# set1: PAIRWISE DE and MERGE CLUSTERS BASED ON SIGNIF DE GENES
		# Use all cells in a cluster, regardless of their orig.ident (i.e. sample)
		tk_merge_clusters_by_de_results <- tk_merge_clusters_by_de(seurat_object, cluster_res_name = cluster_res_name[j], lfc_threshold, num_genes_diff_between_clusters_threshold)
		# Extract results
		merge_data_list <- tk_merge_clusters_by_de_results[["merge_data_list"]]
		iteration_data_list <- tk_merge_clusters_by_de_results[["iteration_data_list"]]
		final_classes <- tk_merge_clusters_by_de_results[["final_classes"]]
	
	
		# MERGING HAS STOPPED, The clusters are now finalized
		# Re-assign the identity factor
		seurat_object@active.ident <- final_classes		
		curr_number_of_clusters_final <- length(levels(factor(seurat_object@active.ident)))		



		# RE-WRITE CLUSTER NAMES (e.g. convert the 1.2.3.4 names to 1) 
		# Rename the final clusters to be regular integers, for plotting and file output
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_orig")]] <- as.character(final_classes)
		# Rename the final clusters to be regular integers, for plotting and file output
		# Make them zero-based index and sorted by size
		nclasses <- length(levels(final_classes))
		cluster_cell_class_factor_final <- factor(forcats::fct_infreq(final_classes), labels = c(1:nclasses - 1))
		names(cluster_cell_class_factor_final) <- names(final_classes)
		seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final")]] <- cluster_cell_class_factor_final	
		# Map the ugly names to final names
		mapping_full <- data.frame(de_merge_orig = final_classes, de_merge_final = cluster_cell_class_factor_final)
		mapping <- unique(mapping_full)
		mapping <- mapping %>% dplyr::arrange(de_merge_final)
		write_tsv(mapping, glue("{orig_wd}/{out_dir}/mapping_cluster_labels.txt"), col_names = TRUE)
		rownames(mapping) <- NULL
		# Export data here
		iteration_data_list <- iteration_data_list[!sapply(iteration_data_list, is.null)]
		final_cluster_labels <- data.frame(cell_id = names(cluster_cell_class_factor_final), clust_orig = class_orig, clust_de_merge_orig = final_classes, clust_de_merge_final = cluster_cell_class_factor_final)
		write_tsv(final_cluster_labels, glue("{orig_wd}/{out_dir}/final_cluster_labels.txt"), col_names = TRUE)
        seurat_object@meta.data$orig.ident <- factor(as.character(seurat_object@meta.data$orig.ident), levels = gtools::mixedsort(as.character(seurat_object@meta.data$orig.ident)), labels = gtools::mixedsort(as.character(seurat_object@meta.data$orig.ident)))

		
		# UPDATE MERGE DATA LIST
		# Get only the final clusters info
		de_merge_orig_names <- names(merge_data_list)
		de_merge_orig_names2 <- str_replace(de_merge_orig_names, "cluster_", "")
		de_merge_orig_names3 <- str_split(de_merge_orig_names2, "_vs_", simplify = FALSE)
		de_data_list <- list()
		for (i in seq_along(de_merge_orig_names3)) {
			if (all(de_merge_orig_names3[[i]] %in% mapping$de_merge_orig)) {
				curr_idx <- length(de_data_list) + 1
				de_data_list[[curr_idx]] <- merge_data_list[[i]]
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "name"] <- "name_orig"
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "cluster_a"] <- "cluster_a_orig"
				names(de_data_list[[curr_idx]])[names(de_data_list[[curr_idx]]) == "cluster_b"] <- "cluster_b_orig"
				de_data_list[[curr_idx]]$cluster_a <- as.character(mapping$de_merge_final[mapping$de_merge_orig %in% merge_data_list[[i]]$cluster_a])
				de_data_list[[curr_idx]]$cluster_b <- as.character(mapping$de_merge_final[mapping$de_merge_orig %in% merge_data_list[[i]]$cluster_b])
				de_data_list[[curr_idx]]$name <- paste0("cluster_", de_data_list[[curr_idx]]$cluster_a, "_vs_", de_data_list[[curr_idx]]$cluster_b)
				de_data_list[[curr_idx]]$de_type <- "de_pairwise"
				names(de_data_list)[curr_idx] <- de_data_list[[curr_idx]]$name
			}
		}
		#lapply(de_data_list, function(x) {print(paste0("final: ", x$name)); print(paste0("orig: ", x$name_orig))})
		
		
		
		
		print("Starting DE set2, DE ONE VS ALL")
		# set2: DE ONE VS ALL
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# testing for cluster "a" vs all other cells (all other clusters combined)
		# Use all cells in a cluster, regardless of their orig.ident (i.e. sample)
		# Add these results to the de_data_list
		print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
		m_vect <- 1:curr_number_of_clusters_final
		out <- mclapply(X = m_vect, FUN = tk_parallel_de_one_vs_all, seurat_object = seurat_object, lfc_threshold = lfc_threshold,
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		# name the list elements
		names(out) <- base::sapply(out, function(x) x[["name"]])
		# update de_data_list	
		# append lists
		length_de_data_list <- length(de_data_list)
		for (i in seq_along(out)) {
			out[[i]]$de_type <- "de_one_vs_all"
			idx <- length_de_data_list + i
			de_data_list[[idx]] <- out[[i]]
			names(de_data_list)[idx] <- out[[i]]$name
		}
		
		
		
		
		
        if (numb_of_samples > 1) {
            print("Multiple samples in analysis. Starting across sample analyses.")

            print("Starting DE set3, PAIRWISE CONSERVED DE LISTS")
            # set3: PAIRWISE CONSERVED DE LISTS
            seurat_object@active.ident <- cluster_cell_class_factor_final
            # Generate pairwise cluster numbering for all cluster comparisons
            combos <- combn(levels(cluster_cell_class_factor_final), 2)
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Run all pairwise cluster comparisons, find DE genes
            k_vect <- 1:dim(combos)[2]
            # Help with parallel lapply:
            # https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
            curr_cluster_res_name <- glue("{cluster_res_name[j]}_de_merge_final")
            out <- mclapply(X = k_vect, FUN = tk_int_parallel_de_pairwise, seurat_object = seurat_object, curr_clusters = cluster_cell_class_factor_final, curr_number_of_clusters = curr_number_of_clusters_final, combos = combos, lfc_threshold = lfc_threshold, de_data_list = de_data_list, meta_data_col = curr_cluster_res_name,
                mc.cores = min(16, as.integer(system("echo $THREADS", intern = TRUE))))
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_pairwise_conserved"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		
        
            print("Starting DE set4, ONE VS ALL CONSERVED DE LISTS")
            # set4: ONE VS ALL CONSERVED DE LISTS
            seurat_object@active.ident <- cluster_cell_class_factor_final
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Run all pairwise cluster comparisons, find DE genes
            m_vect <- 1:curr_number_of_clusters_final
            # Help with parallel lapply:
            # https://stackoverflow.com/questions/15852482/mclapply-additional-arguments
            curr_cluster_res_name <- glue("{cluster_res_name[j]}_de_merge_final")
            out <- mclapply(X = m_vect, FUN = tk_int_parallel_de_one_vs_all, seurat_object = seurat_object, lfc_threshold = lfc_threshold, meta_data_col = curr_cluster_res_name,
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_one_vs_all_conserved"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		

        


            print("Starting DE set5, DE LISTS ACROSS SAMPLES -- WITHIN SAME CLUSTER")
            # set5: DE LISTS ACROSS SAMPLES -- WITHIN SAME CLUSTER
            print(glue("Using Seurat Object Assay: {DefaultAssay(seurat_object)}"))
            # Create a new identity labeling scheme that combines the original sample name and final cluster name
            col_cluster <- colnames(seurat_object@meta.data)[str_detect(colnames(seurat_object@meta.data), "_de_merge_final$")]
            col_sample <- "orig.ident"
            new_levels <- tidyr::crossing(seurat_object@meta.data[[col_cluster]], seurat_object@meta.data[[col_sample]]) %>%
                dplyr::rename(col_cluster = `seurat_object@meta.data[[col_cluster]]`, col_sample = `seurat_object@meta.data[[col_sample]]`) %>%
                dplyr::mutate(new_levels = as.character(glue("{col_cluster}.{col_sample}")))
            cluster_sample <- paste0(seurat_object@meta.data[[col_cluster]], ".", seurat_object@meta.data[[col_sample]])
            cluster_sample <- factor(cluster_sample, levels = new_levels$new_levels)
            cluster_sample <- fct_relevel(cluster_sample, gtools::mixedsort)
            seurat_object@meta.data[[paste0(cluster_res_name[j], "_de_merge_final_cluster_sample")]] <- cluster_sample
            names(cluster_sample) <- rownames(seurat_object@meta.data)
            seurat_object@active.ident <- cluster_sample
            cluster_sample_tbl <- tibble(cluster = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[1]),
                sample = sapply(str_split(levels(cluster_sample), fixed("."), simplify = FALSE), function(x) x[2]))
            uniq_clust <- unique(cluster_sample_tbl$cluster)
            cluster_sample_list <- list()
            for (i in seq_along(uniq_clust)) {
                curr <- cluster_sample_tbl %>% filter(cluster == uniq_clust[i])
                cluster_sample_list[[i]] <- list(cluster = uniq_clust[i], sample = curr$sample)
            }
            t_vect <- 1:curr_number_of_clusters_final
            out1 <- mclapply(X = t_vect, FUN = tk_int_parallel_de_across_samples, seurat_object = seurat_object, cluster_sample_list = cluster_sample_list, lfc_threshold = lfc_threshold,
                mc.cores = min(16, as.integer(system("echo $THREADS", intern = TRUE))))
            # Flatten one level of the lists
            out <- unlist(out1, recursive = FALSE)
            # name the list elements
            names(out) <- base::sapply(out, function(x) x[["name"]])
            # Find length of prev merge data list
            length_de_data_list <- length(de_data_list)
            # append lists
            for (i in seq_along(out)) {
                out[[i]]$de_type <- "de_across_samples"
                idx <- length_de_data_list + i
                de_data_list[[idx]] <- out[[i]]
                names(de_data_list)[idx] <- out[[i]]$name
            }		
        
        
        } else if (numb_of_samples == 1) {
            print("One sample in analysis. Skipping across sample analyses.")
        }





		print("Making UMAP plots")
		# UMAP PLOTS
		# Plot original UMAP
		umap_plot_list <- list()
		# Re-assign the identity factor
		# Make this a nicely sorted factor of integers, zero-based
		class_orig_factor <- factor(as.numeric(class_orig) - 1)
		names(class_orig_factor) <- names(class_orig)
		seurat_object@active.ident <- class_orig_factor
		plot_data <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
		print("Using seurat_object@reductions$umap@cell.embeddings for UMAP")
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[1]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nInitial pre-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[1] <- paste0(cluster_res_name[j], "_before_de_merge")
		# Plot final UMAP
		# Re-assign the identity factor
		seurat_object@active.ident <- cluster_cell_class_factor_final
		# get centers for labels
		plot_data$active.ident <- as.factor(seurat_object@active.ident)
		plot_data %>%
			dplyr::group_by(active.ident) %>%
			summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> centers
		umap_plot_list[[2]] <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = active.ident)) +
			geom_point(size = 0.5) +
			labs(color = "Cluster", title = paste0(seurat_object_name, "\nCluster resolution: ", cluster_res_name[j], "\nFinal post-DE merging")) +
			geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0) +
			geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
			theme_cowplot() +
			theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5))
		names(umap_plot_list)[2] <- paste0(cluster_res_name[j], "_final")
		pdf("umap_final_de_merge.pdf", width = 6, height = 6)
		print(umap_plot_list[[2]])
		dev.off()
		pdf("umap_original.pdf", width = 6, height = 6)
		print(umap_plot_list[[1]])
		dev.off()
		# Get the final iteration counter for plotting
		iteration_counter <- length(names(iteration_data_list))
		iteration_counter_final <- iteration_counter - 1
		pdf("umap_original_and_final_de_merge.pdf", width = 10, height = 6)
		grid.arrange(
			grobs = umap_plot_list,
			nrow = 1,
			top = paste0("UMAP plots, Pre- and Post-DE merging.\n", iteration_counter_final , " cluster merging steps required."))
		dev.off()
		
		
		
# 		# Print the names of de_data_list
# 		lapply(de_data_list, function(x) {print(x$name)})
# 		lapply(de_data_list, function(x) {names(x)})
# 		lapply(de_data_list, function(x) {print(head(x$de_signif))})
# 		# Not tibbles
# 		cluster_5_vs_all
# 		cluster_4_vs_5
		
		

        print("Starting most SCT variable gene plots")
        if (!is.null(seurat_object@assays$SCT)) {
            # Requires SCT assay
            most_variable_features_sorted_n50 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 50) %>%
                base::rownames(.) 
            most_variable_features_sorted_n100 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 100) %>%
                base::rownames(.) 
            most_variable_features_sorted_n200 <- seurat_object@assays$SCT@SCTModel.list$model1@feature.attributes %>%
                dplyr::arrange(desc(residual_variance)) %>%
                dplyr::slice_head(n = 200) %>%
                base::rownames(.)
            # Heatmap (top most variable genes across all cells)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"))) {dir.create(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/heatmap_most_var_genes"))
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n50, filename_prefix = "heatmap_top50_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 50 most variable genes across all cells\nEach row independently scaled")
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n100, filename_prefix = "heatmap_top100_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 100 most variable genes across all cells\nEach row independently scaled")
            tk_cluster_heatmap(seurat_object, gene_names = most_variable_features_sorted_n200, filename_prefix = "heatmap_top200_var", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 200 most variable genes across all cells\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
        }





		print("Starting set1 plots: de_pairwise")
		# Set1: GET PAIRWISE PLOTS
		r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_pairwise")
		out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_pairwise_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		names(out) <- base::sapply(out, function(x) x[["name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])       
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		tk_int_parallel_de_pairwise_plots_results <- out
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# The openxlsx package is using an old zip method, so you might get warnings (but should still work)
		# Note: zip::zip() is deprecated, please use zip::zipr() instead
		# https://github.com/awalker89/openxlsx/issues/454
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_pairwise.xlsx")
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_pairwise_signif.xlsx")
        # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
        comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
        comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_pairwise/umap_gene_expression_flattened"))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_pairwise/umap_gene_expression_flattened"))
		# Heatmap
		# Get only top DE genes from each comparison
		top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
        top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
        # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
        top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
        top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
		# Heatmap (top up and and down 5 genes)
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_and_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))    
		# Heatmap (top up or down genes, farthest from the mean)       
 		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap_top_up_or_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))             
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/heatmap"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_pairwise_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))
		# Dotplots
		seurat_object@active.ident <- cluster_cell_class_factor_final
		if (!file.exists(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))}
		setwd(glue("{orig_wd}/{out_dir}/de_pairwise/dotplot_gene_expression"))
		tk_int_dotplots(seurat_object, tk_int_parallel_de_pairwise_plots_results)
		setwd(glue("{orig_wd}/{out_dir}"))
		





		print("Starting set2 plots: de_one_vs_all")
		# Set2: GET ONE_VS_ALL PLOTS AND LISTS
		r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_one_vs_all")
		out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_one_vs_all_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
			mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
		names(out) <- base::sapply(out, function(x) x[["name"]])
		excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
		excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
		signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
		tk_int_parallel_de_one_vs_all_plots_results <- out
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
		excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
		ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
		excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
		# The openxlsx package is using an old zip method, so you might get warnings (but should still work)
		# Note: zip::zip() is deprecated, please use zip::zipr() instead
		# https://github.com/awalker89/openxlsx/issues/454
		# xlsx tab names cannot be longer than 31 characters, clip them
		excel_list_of_df_de_XLSX <- excel_list_of_df_de
		names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
		names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
		openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_one_vs_all.xlsx")
		openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_one_vs_all_signif.xlsx")
        # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
        comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
        comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_one_vs_all/umap_gene_expression_flattened"))
        base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_one_vs_all/umap_gene_expression_flattened"))		
		# Heatmap
		# Get only top DE genes from each comparison
		top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
        top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
        # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
        top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
        top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
		# Heatmap (top up and and down 5 genes)
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_and_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))    
		# Heatmap (top up or down genes, farthest from the mean)       
 		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap_top_up_or_down_reg"))
		tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))             
		if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"), recursive = TRUE)}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/heatmap"))
		tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_one_vs_all_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
		setwd(glue("{orig_wd}/{out_dir}"))
		# Dotplots
		seurat_object@active.ident <- cluster_cell_class_factor_final
		if (!file.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))}
		setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all/dotplot_gene_expression"))
		tk_int_dotplots(seurat_object, tk_int_parallel_de_one_vs_all_plots_results)
		setwd(glue("{orig_wd}/{out_dir}"))
			




		
        if (numb_of_samples > 1) {
            print("Multiple samples in analysis. Starting across sample plots.")

        
        
            print("Starting set3 plots: de_pairwise_conserved")
            # Set3: GET PAIRWISE CONSERVED PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_pairwise_conserved")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_pairwise_conserved_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_pairwise_conserved_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_pairwise_cons.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_pairwise_cons_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_pairwise_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_pairwise_cons/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_pairwise_cons/umap_gene_expression_flattened"))		
            # Heatmap 
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg_cons, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg_cons, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_pairwise_cons_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_pairwise_cons/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_pairwise_conserved_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))

                
        
        
        




            print("Starting set4 plots: de_one_vs_all_conserved")
            # Set4: GET ONE_VS_ALL CONSERVED PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_one_vs_all_conserved")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_one_vs_all_conserved_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_one_vs_all_conserved_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_one_vs_all_cons.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_one_vs_all_cons_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_one_vs_all_cons/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_one_vs_all_cons/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_one_vs_all_cons/umap_gene_expression_flattened"))
            # Heatmap
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg_cons, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg_cons, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))             
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_one_vs_all_cons_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_one_vs_all_cons/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_one_vs_all_conserved_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))

                
        
        
        
        
            print("Starting set5 plots: de_across_samples")
            # Set5: DE ACROSS SAMPLES -- within same cluster -- PLOTS AND LISTS
            r_vect <- which(sapply(de_data_list, function(x) x$de_type) == "de_across_samples")
            out <- mclapply(X = r_vect, FUN = tk_int_parallel_de_across_samples_plots, seurat_object = seurat_object, de_data_list = de_data_list, umap_plot = umap_plot_list[[2]],
                mc.cores = as.integer(system("echo $THREADS", intern = TRUE)))
            names(out) <- base::sapply(out, function(x) x[["name"]])
            excel_list_of_df_de <- lapply(out, function(x) x[["excel_list_of_df_de"]])
            excel_list_of_df_de_signif <- lapply(out, function(x) x[["excel_list_of_df_de_signif"]])
            signif_genes_across_all_clusters <- unique(unlist(sapply(out, function(x) x[["signif_genes_across_all_clusters"]])))
            tk_int_parallel_de_across_samples_plots_results <- out
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de))
            excel_list_of_df_de <- excel_list_of_df_de[ordered_names]
            ordered_names <- gtools::mixedsort(names(excel_list_of_df_de_signif))
            excel_list_of_df_de_signif <- excel_list_of_df_de_signif[ordered_names]
            # The openxlsx package is using an old zip method, so you might get warnings (but should still work)
            # Note: zip::zip() is deprecated, please use zip::zipr() instead
            # https://github.com/awalker89/openxlsx/issues/454
            # xlsx tab names cannot be longer than 31 characters, clip them
            excel_list_of_df_de_XLSX <- excel_list_of_df_de
            names(excel_list_of_df_de_XLSX) <- names(excel_list_of_df_de_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            excel_list_of_df_de_signif_XLSX <- excel_list_of_df_de_signif
            names(excel_list_of_df_de_signif_XLSX) <- names(excel_list_of_df_de_signif_XLSX) %>% stringr::str_trunc(width = 30, ellipsis = "")
            openxlsx::write.xlsx(excel_list_of_df_de_XLSX, file = "de_across_samples.xlsx")
            openxlsx::write.xlsx(excel_list_of_df_de_signif_XLSX, file = "de_across_samples_signif.xlsx")
            # Move original UMAP pdfs to the flattened folder, then symlink the PDF back to the original location
            comparison_dirs_up <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_up")))
            comparison_dirs_dn <- base::sapply(X = r_vect, function(r) as.character(paste0("./de_across_samples/umap_gene_expression/", de_data_list[[r]]$name, "_de_results_signif_geneplots_dn")))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_up[r], new_dir = "./de_across_samples/umap_gene_expression_flattened"))
            base::sapply(X = r_vect, function(r) tk_move_and_symlink(original_dir = comparison_dirs_dn[r], new_dir = "./de_across_samples/umap_gene_expression_flattened"))            
            # Heatmap
            # Get only top DE genes from each comparison
            top_up_and_down_reg <- lapply(excel_list_of_df_de_signif, tk_top_up_and_down_reg, n = 5)
            top_up_and_down_reg_unique <- unique(unlist(top_up_and_down_reg)) 
            # Get the only top genes farthest from the mean, regardless if they are pos/neg log2FC
            top_up_and_down_reg_genes <- lapply(excel_list_of_df_de_signif, tk_top_up_or_down_reg, n = 10)
            top_up_and_down_reg_genes_unique <- unique(unlist(top_up_and_down_reg_genes))
            # Heatmap (top up and and down 5 genes)
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_and_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_unique, filename_prefix = "heatmap_top_up_and_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 5 up- and top 5 down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))    
            # Heatmap (top up or down genes, farthest from the mean)       
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap_top_up_or_down_reg"))
            tk_cluster_heatmap(seurat_object, gene_names = top_up_and_down_reg_genes_unique, filename_prefix = "heatmap_top_up_or_down_reg", active.ident = cluster_cell_class_factor_final, heatmap_title = "Top 10 up- or down-regulated genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))             
            if (!dir.exists(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"), recursive = TRUE)}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/heatmap"))
            tk_cluster_heatmap(seurat_object, gene_names = signif_genes_across_all_clusters, filename_prefix = "de_across_samples_signif_heatmap", active.ident = cluster_cell_class_factor_final, heatmap_title = "All significant DE genes from each comparison\nEach row independently scaled")
            setwd(glue("{orig_wd}/{out_dir}"))
            # Dotplots
            seurat_object@active.ident <- cluster_cell_class_factor_final
            if (!file.exists(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))) {dir.create(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))}
            setwd(glue("{orig_wd}/{out_dir}/de_across_samples/dotplot_gene_expression"))
            tk_int_dotplots(seurat_object, tk_int_parallel_de_across_samples_plots_results)
            setwd(glue("{orig_wd}/{out_dir}"))


        } else if (numb_of_samples == 1) {
            print("One sample in analysis. Skipping across sample plots.")
        }



        ## EXPORT DATA
		DefaultAssay(seurat_object) <- orig_assay

        saveRDS(seurat_object@meta.data, "cell_metadata.rds")
        write_tsv(as_tibble(seurat_object@meta.data, rownames = "cell_id"), "cell_metadata.txt")
		# Export results to list of list
		curr_return_list <- list(cluster_cell_class_factor_final = cluster_cell_class_factor_final,
			final_cluster_labels = final_cluster_labels, 
			iteration_data_list = iteration_data_list, 
			merge_data_list = merge_data_list, 
			de_data_list = de_data_list,
			mapping = mapping, 
			mapping_full = mapping_full, 
			lfc_threshold = lfc_threshold, 
			num_genes_diff_between_clusters_threshold = num_genes_diff_between_clusters_threshold,
			umap_plot_list = umap_plot_list,
			seurat_object_meta_data = seurat_object@meta.data)
		return_list[[j]] <- curr_return_list
		names(return_list)[j] <- cluster_res_name[j]
	}
	return(return_list)
}







#' @export
tk_int_parallel_de_across_samples <- function(t, seurat_object = seurat_object, cluster_sample_list = cluster_sample_list, lfc_threshold = lfc_threshold) {
	# DE LISTS ACROSS SAMPLES, within one cluster
	curr_cluster <- cluster_sample_list[[t]]$cluster
	curr_sample_vect <- cluster_sample_list[[t]]$sample
	
	combos <- combn(curr_sample_vect, 2)
	all_results <- list()
	for (i in seq_len(ncol(combos))) {
		curr_sample_a <- combos[1, i]
		curr_sample_b <- combos[2, i]
		ident_1 <- paste0(curr_cluster, ".", curr_sample_a)
		ident_2 <- paste0(curr_cluster, ".", curr_sample_b)
		curr_name <- paste0(curr_cluster, ".", curr_sample_a,"_vs_", curr_cluster, ".", curr_sample_b)
		# Initialize to "not" skip 
		sample_specific_cells <- seurat_object@active.ident == ident_1
		ident_1_cells_num <- dim(GetAssayData(seurat_object, slot = "data")[, sample_specific_cells, drop = FALSE])[2]
		sample_specific_cells <- seurat_object@active.ident == ident_2
		ident_2_cells_num <- dim(GetAssayData(seurat_object, slot = "data")[, sample_specific_cells, drop = FALSE])[2]
		if (ident_1_cells_num < 3 | ident_2_cells_num < 3) {
			# Too few files for DE analysis, skip DE analysis
			print(paste0("Too few files for DE analysis, skipping DE analysis for ident: ", curr_name))
			# Create empty tibble
			de_blank <- matrix(0, 0, ncol = 6)
			colnames(de_blank) <- c("symbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
			de_blank <- as_tibble(de_blank)
			de_blank <- de_blank %>%
				mutate_all(as.character)
			de <- de_blank
			de_signif <- de_blank
			number_of_signif_genes <- 0
			results <- list(name = curr_name, cluster = curr_cluster, sample_a = curr_sample_a, sample_b = curr_sample_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
		} else {
			print(paste0(curr_name, " DE analysis."))
			# Find DE genes
			de <- FindMarkers(object = seurat_object, slot = "data", ident.1 = ident_1, ident.2 = ident_2, min.pct = 0, test.use = "wilcox", logfc.threshold = 0)
			de <- de %>%
				tibble::rownames_to_column(var = "symbol") %>%
				as_tibble()
			# Pass thresholds
			# Sort the table by the first log2FC column group
			de_signif <- de %>%
				dplyr::filter(abs(avg_log2FC) >= lfc_threshold & p_val_adj <= 0.01) %>%
				dplyr::arrange(desc(avg_log2FC))	
			# How many genes are significant?
			number_of_signif_genes <- dim(de_signif)[1]
			results <- list(name = curr_name, cluster = curr_cluster, sample_a = curr_sample_a, sample_b = curr_sample_b, de = de, de_signif = de_signif, number_of_signif_genes = number_of_signif_genes)
		}
		all_results[[i]] <- results
		names(all_results)[i] <- curr_name
	}
	return(all_results)
}










# ---------------------------------------------------------------------
# Todds version of the FeaturePlot()
# ---------------------------------------------------------------------


# REPLACE THIS (with my function):
#ind_plots_list[[i]] <- FeaturePlot(object = seurat_object, slot = "data", features = genes_list_short[i], order = TRUE, cols = colorRampPalette(c("white", "red"))(100), reduction = "umap", combine = TRUE)

#' @export
tk_feature_plot <- function(object = seurat_object, slot = "data", feature = NULL, subtitle = NULL, facet_by = NULL) {
    curr_assay <- DefaultAssay(object)
    # Get clean name for data slot
    if (curr_assay == "SCT" & slot == "counts") {
        slot_name <- "Calc raw counts"
    } else if (curr_assay == "SCT" & slot == "data") {
        slot_name <- "Log calc counts"
    } else if (curr_assay == "SCT" & slot == "scale.data") {    
        slot_name <- "Norm data" # Pearson Residual
    } else {
        slot_name <- "Unknown"
    }

    expr_data <- Seurat::FetchData(object = object, vars = feature, slot = slot)

    umap_data <- as.data.frame(object@reductions$umap@cell.embeddings)
    stopifnot(all(rownames(umap_data) == rownames(expr_data)))
    plot_data <- bind_cols(umap_data, expr_data, object@meta.data)

    plot_data$active.ident <- as.factor(object@active.ident)
    centers <- plot_data %>%
        dplyr::group_by(active.ident) %>%
        dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
        
    if (feature %in% colnames(plot_data)) {
        plot <- plot_data %>%
            # Plot the highest expression dots on top of lower expression dots
            dplyr::arrange(!!rlang::sym(feature)) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(feature))) +
                geom_point(size = 0.5) +
                scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
                # Ensure the plots stay square
                theme_cowplot() +
                theme(aspect.ratio = 1, 
                    plot.title = element_text(hjust = 0.5)) +
                #geom_text(data = centers, mapping = aes(label = active.ident), size = 4, color = "black") +
                #labs(title = glue("Project: {object@project.name}\nAssay: {DefaultAssay(object)}; Data slot: {slot_name}\n{feature}\n{subtitle}"))
                labs(title = glue("{feature}"))
        if (!is.null(facet_by)) {
            plot <- plot +
                facet_grid(rows = facet_by)
        }
    } else {
        # feature is not found in dataset
        warning(paste0(feature, " not found in data, creating blank plot"))
        # Create an empty plot!
        plot <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_cowplot() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{feature}"))
    }
    return(plot)
}
# Example:
# Default_Assay(curr_seurat_object) <- "SCT"
# p <- tk_feature_plot(object = curr_seurat_object_keep, feature = "Cd3e", subtitle = "Cluster Res 0.9")
# pdf()
# print(p)
# dev.off()




# ---------------------------------------------------------------------
# Todds version of the VlnPlot()
# ---------------------------------------------------------------------


#' @export
tk_violin_feature_plot <- function(object = seurat_object, slot = "data", feature = NULL, subtitle = NULL, facet_by = NULL) {
    curr_assay <- DefaultAssay(object)
    # Get clean name for data slot
    if (curr_assay == "SCT" & slot == "counts") {
        slot_name <- "Calc raw counts"
    } else if (curr_assay == "SCT" & slot == "data") {
        slot_name <- "Log calc counts"
    } else if (curr_assay == "SCT" & slot == "scale.data") {    
        slot_name <- "Norm data" # Pearson Residual
    } else {
        slot_name <- "Unknown"
    }

    expr_data <- Seurat::FetchData(object = object, vars = feature, slot = slot)
    expr_data$active.ident <- factor(object@active.ident, levels = gtools::mixedsort(levels(object@active.ident)), labels = gtools::mixedsort(levels(object@active.ident)))

    plot_data <- bind_cols(expr_data, object@meta.data)
    
    if (feature %in% colnames(plot_data)) {
        plot <- plot_data %>%
            ggplot(aes(x = active.ident, y = !!rlang::sym(feature), fill = active.ident)) +
                geom_violin(scale = "width", adjust = 1, trim = TRUE) +
                # Ensure the plots stay square
                theme_cowplot() +
                theme(aspect.ratio = 1, 
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "none") +
                labs(title = glue("{feature}"), x = "Cluster", y = "Expression", fill = "Cluster")
        if (!is.null(facet_by)) {
            plot <- plot +
                facet_grid(rows = facet_by)
        }
    } else {
        # feature is not found in dataset
        warning(paste0(feature, " not found in data, creating blank plot"))
        # Create an empty plot!
        plot <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_cowplot() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{feature}"))
    }
    return(plot)
}
# Example:
# Default_Assay(curr_seurat_object) <- "SCT"
# p <- tk_violin_feature_plot(object = curr_seurat_object_keep, feature = "Cd3e", subtitle = "Cluster Res 0.9")
# pdf()
# print(p)
# dev.off()




# ---------------------------------------------------------------------
# Todds version of the boxplot
# ---------------------------------------------------------------------


#' @export
tk_box_feature_plot <- function(object = seurat_object, slot = "data", feature = NULL, subtitle = NULL, facet_by = NULL) {
    curr_assay <- DefaultAssay(object)
    # Get clean name for data slot
    if (curr_assay == "SCT" & slot == "counts") {
        slot_name <- "Calc raw counts"
    } else if (curr_assay == "SCT" & slot == "data") {
        slot_name <- "Log calc counts"
    } else if (curr_assay == "SCT" & slot == "scale.data") {    
        slot_name <- "Norm data" # Pearson Residual
    } else {
        slot_name <- "Unknown"
    }

    expr_data <- Seurat::FetchData(object = object, vars = feature, slot = slot)
    expr_data$active.ident <- factor(object@active.ident, levels = gtools::mixedsort(levels(object@active.ident)), labels = gtools::mixedsort(levels(object@active.ident)))
    
    plot_data <- bind_cols(expr_data, object@meta.data)
     
    if (feature %in% colnames(plot_data)) {
        plot <- plot_data %>%
            ggplot(aes(x = active.ident, y = !!rlang::sym(feature), fill = active.ident)) +
                geom_boxplot() +
                # Ensure the plots stay square
                theme_cowplot() +
                theme(aspect.ratio = 1, 
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "none") +
                labs(title = glue("{feature}"), x = "Cluster", y = "Expression", fill = "Cluster")
        if (!is.null(facet_by)) {
            plot <- plot +
                facet_grid(rows = facet_by)
        }
    } else {
        # feature is not found in dataset
        warning(paste0(feature, " not found in data, creating blank plot"))
        # Create an empty plot!
        plot <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_cowplot() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{feature}"))
    }
    return(plot)
}
# Example:
# Default_Assay(curr_seurat_object) <- "SCT"
# p <- tk_violin_feature_plot(object = curr_seurat_object_keep, feature = "Cd3e", subtitle = "Cluster Res 0.9")
# pdf()
# print(p)
# dev.off()



# ---------------------------------------------------------------------
# Todds version of the RidgePlot()
# ---------------------------------------------------------------------


#' @export
tk_ridge_feature_plot <- function(object = seurat_object, slot = "data", feature = NULL, subtitle = NULL, facet_by = NULL) {
    curr_assay <- DefaultAssay(object)
    # Get clean name for data slot
    if (curr_assay == "SCT" & slot == "counts") {
        slot_name <- "Calc raw counts"
    } else if (curr_assay == "SCT" & slot == "data") {
        slot_name <- "Log calc counts"
    } else if (curr_assay == "SCT" & slot == "scale.data") {    
        slot_name <- "Norm data" # Pearson Residual
    } else {
        slot_name <- "Unknown"
    }

    expr_data <- Seurat::FetchData(object = object, vars = feature, slot = slot)
    expr_data$active.ident <- factor(object@active.ident, levels = gtools::mixedsort(levels(object@active.ident)), labels = gtools::mixedsort(levels(object@active.ident)))
    
    plot_data <- bind_cols(expr_data, object@meta.data)
    
    if (feature %in% colnames(plot_data)) {
        plot <- plot_data %>%
            ggplot(aes(x = !!rlang::sym(feature), y = active.ident, fill = active.ident)) +
                ggridges::geom_density_ridges(scale = 4) +
                ggridges::theme_ridges() +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0)) +
                coord_cartesian(clip = "off") +
                # Ensure the plots stay square
                theme_cowplot() +
                theme(aspect.ratio = 1, 
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "none") +
                labs(title = glue("{feature}"), x = "Expression", y = "Cluster", fill = "Cluster")
        if (!is.null(facet_by)) {
            plot <- plot +
                facet_grid(rows = facet_by)
        }
    } else {
        # feature is not found in dataset
        warning(paste0(feature, " not found in data, creating blank plot"))
        # Create an empty plot!
        plot <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_cowplot() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{feature}"))
    }
    return(plot)
}
# Example:
# Default_Assay(curr_seurat_object) <- "SCT"
# p <- tk_ridge_feature_plot(object = curr_seurat_object, feature = "Ccl5", subtitle = "Cluster Res 0.9")
# pdf("ridge.pdf")
# print(p)
# dev.off()






# ---------------------------------------------------------------------
# Todds version of the Dotplot()
# ---------------------------------------------------------------------


#' @export
tk_dot_feature_plot <- function(object = seurat_object, slot = "data", feature = NULL, subtitle = NULL, facet_by = NULL) {
    curr_assay <- DefaultAssay(object)
    # Get clean name for data slot
    if (curr_assay == "SCT" & slot == "counts") {
        slot_name <- "Calc raw counts"
    } else if (curr_assay == "SCT" & slot == "data") {
        slot_name <- "Log calc counts"
    } else if (curr_assay == "SCT" & slot == "scale.data") {    
        slot_name <- "Norm data" # Pearson Residual
    } else {
        slot_name <- "Unknown"
    }

    expr_data <- Seurat::FetchData(object = object, vars = feature, slot = slot)
    expr_data$active.ident <- factor(object@active.ident, levels = gtools::mixedsort(levels(object@active.ident)), labels = gtools::mixedsort(levels(object@active.ident)))
    object@active.ident <- factor(object@active.ident, levels = gtools::mixedsort(levels(object@active.ident)), labels = gtools::mixedsort(levels(object@active.ident)))
    
    plot_data <- bind_cols(expr_data, object@meta.data)
    
    if (feature %in% colnames(plot_data)) {
        # NOTE: regardless of data slot, the "average gene expression" is "scale()" across clusters
        # https://github.com/satijalab/seurat/blob/9843b843ed0c3429d86203011fda00badeb29c2e/R/visualization.R#L3464
        # The size of the dots ("percent expressed") will vary depending on which genes are included in the plot
        if (!is.null(facet_by)) {
            plot <- DotPlot(object, features = feature, dot.scale = 8, group.by = facet_by) + RotatedAxis() + 
                scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
                labs(title = glue("{feature}"), x = "Gene", y = "Cluster", fill = "Cluster") +
                    theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        plot.title = element_text(hjust = 0.5))
        } else {
            plot <- DotPlot(object, features = feature, dot.scale = 8) + RotatedAxis() + 
                scale_color_gradientn(name = "Exprs", colors = c("gray95", rev(viridis::plasma(100)))) +
                labs(title = glue("{feature}"), x = "Gene", y = "Cluster", fill = "Cluster") +
                    theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        plot.title = element_text(hjust = 0.5))
        }
    } else {
        # feature is not found in dataset
        warning(paste0(feature, " not found in data, creating blank plot"))
        # Create an empty plot!
        plot <- ggplot() +
            annotate("text", x = 1, y = 1, label = "Nothing to plot for feature") +
            theme_cowplot() +
            theme(aspect.ratio = 1, 
                plot.title = element_text(hjust = 0.5)) +
            labs(title = glue("{feature}"))
    }
    return(plot)
}
# Example:
# Default_Assay(curr_seurat_object) <- "SCT"
# p <- tk_dot_feature_plot(object = curr_seurat_object, feature = "Ccl5", subtitle = "Cluster Res 0.9")
# pdf("ridge.pdf")
# print(p)
# dev.off()












# ---------------------------------------------------------------------
# Misc heatmap related functions
# ---------------------------------------------------------------------




#' @export
tk_top_up_and_down_reg <- function(de_signif_df, n = 5) {
    top_g <- de_signif_df %>%
        #dplyr::filter(pct.1 >= 0.05 & pct.2 >= 0.05) %>%
        dplyr::slice_head(n = n) %>%
        pull(symbol)
    bottom_g <- de_signif_df %>%
        dplyr::filter(pct.1 >= 0.05 & pct.2 >= 0.05) %>%
        dplyr::slice_tail(n = n) %>%
        pull(symbol)   

    top_bot_g <- unique(c(top_g, bottom_g))
    return(top_bot_g)
}
#' @export
tk_top_up_and_down_reg_cons <- function(de_signif_df, n = 5) {
    top_g <- de_signif_df %>%
        dplyr::slice_head(n = n) %>%
        pull(symbol)
    bottom_g <- de_signif_df %>%
        dplyr::slice_tail(n = n) %>%
        pull(symbol)   

    top_bot_g <- unique(c(top_g, bottom_g))
    return(top_bot_g)
}



#' @export
tk_top_up_or_down_reg <- function(de_signif_df, n = 5) {
    df <- de_signif_df #%>%
        #dplyr::filter(pct.1 >= 0.05 & pct.2 >= 0.05)
    if (dim(df)[1] >= 1) {
        farthest_from_mean <- df$symbol[head(order(abs(df$avg_log2FC - mean(df$avg_log2FC)), decreasing = TRUE), n)]
    } else {
        farthest_from_mean <- NULL
    }
    return(farthest_from_mean)
}

#' @export
tk_top_up_or_down_reg_cons <- function(de_signif_df, n = 5) {
    df <- de_signif_df
    if (dim(df)[1] >= 1) {
        farthest_from_mean <- df$symbol[head(order(abs(df$avg_log2FC - mean(df$avg_log2FC)), decreasing = TRUE), n)]
   } else {
        farthest_from_mean <- NULL
    }
   return(farthest_from_mean)
}











# ---------------------------------------------------------------------
# Symlink feature plots
# ---------------------------------------------------------------------

#' @export
tk_move_and_symlink <- function(original_dir, new_dir) {
    # (1) move all files from original_dir to new_dir (if they are not already in new_dir)
    # (2) symlink the real files in new_dir back to original_dir
    original_dir2 <- normalizePath(original_dir)
    
    
    if (dir.exists(original_dir2)) {
        # Make sure new_dir exists
        if (!dir.exists(new_dir)) {dir.create(new_dir, recursive = TRUE)}
        new_dir2 <- normalizePath(new_dir)
        # Get vector of files
        files <- list.files(original_dir2, full.names = TRUE)
        filenames <- basename(files)
        if (length(files) > 0) {
            for (i in seq_along(files)) {
                if (!file.exists(paste0(new_dir2, "/", filenames[i]))) {
                    # If file not in new_dir2, move it there
                    fs::file_move(files[i], paste0(new_dir2, "/", filenames[i]))
                } else {
                    # Delete any files that did not need to be moved
                    file.remove(files[i])
                }
            }
            # Create symlinks instead
            for (i in seq_along(files)) {
                file.symlink(paste0(new_dir2, "/", filenames[i]), original_dir2)
            }
        }
    } else {
        warning(paste(original_dir2, "does not exist, cannot create symlinks"))
    }
}

# Provide relative path names
# tk_move_and_symlink("cluster_1_vs_2_up", "../umap_gene_expression_flattened")






#######################################################################
# Tweak Dotplot from Seurat
#######################################################################

# https://github.com/satijalab/seurat/blob/4e868fcde49dc0a3df47f94f5fb54a421bfdf7bc/R/visualization.R#L3331-L3597



#' Dot plot visualization
#'
#' Intuitive way of visualizing how feature expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level
#' across all cells within a class (blue is high).
#'
#' @param object Seurat object
#' @param assay Name of assay to use, defaults to the active assay
#' @param features Input vector of features, or named list of feature vectors
#' if feature-grouped panels are desired (replicates the functionality of the
#' old SplitDotPlotGG)
#' @param cols Colors to plot: the name of a palette from
#' \code{RColorBrewer::brewer.pal.info}, a pair of colors defining a gradient,
#' or 3+ colors defining multiple gradients (if split.by is set)
#' @param col.min Minimum scaled average expression threshold (everything
#' smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param idents Identity classes to include in plot (default is all)
#' @param group.by Factor to group the cells by
#' @param split.by Factor to split the groups by (replicates the functionality
#' of the old SplitDotPlotGG);
#' see \code{\link{FetchData}} for more details
#' @param cluster.idents Whether to order identities by hierarchical clusters
#' based on given features, default is FALSE
#' @param scale Determine whether the data is scaled, TRUE for default
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#'
#' @return A ggplot object
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius
#' theme element_blank labs scale_color_identity scale_color_distiller
#' scale_color_gradient guides guide_legend guide_colorbar
#' facet_grid unit
#' @importFrom scattermore geom_scattermore
#' @importFrom stats dist hclust
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @export
#' @concept visualization
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features = cd_genes)
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' DotPlot(object = pbmc_small, features = cd_genes, split.by = 'groups')
#'
tk_DotPlot <- function(
  object,
  assay = NULL,
  features,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))


  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = tk_PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  saveRDS(data.plot, "data.plot.rds")
  saveRDS(scale.func, "scale.func.rds")
  return(plot)
}


#' tk_PercentAbove
#' https://github.com/satijalab/seurat/blob/4e868fcde49dc0a3df47f94f5fb54a421bfdf7bc/R/utilities.R#L2127-L2136
#' This was copied exactly from Seurat. The function is not exported from Seurat, so I included it here.
#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#' @export
tk_PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}




