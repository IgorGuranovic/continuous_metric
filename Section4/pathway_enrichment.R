# Clear workspace
rm(list = ls())

# Define file paths
npcs_data_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_cap.csv"
stroma_data_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_stroma.csv"
roi_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/ROI_notes.csv"
output_dir <- "D:/Igor/Figures/Section4"
dir.create(output_dir, showWarnings = FALSE)

# Load libraries
library(Seurat)
library(pheatmap)
library(viridis)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(parallel)
library(doParallel)

# Set up parallel backend
num_cores <- 30
registerDoParallel(cores = num_cores)

# Load ROI notes
roi_notes <- read.csv(roi_path)

# Function to create heatmap
create_heatmap <- function(data_path, output_path, title, dataset_name) {
  data_matrix <- as.matrix(read.csv(data_path, row.names = 1))
  colnames(data_matrix) <- gsub("_count$", "", colnames(data_matrix))
  roi_notes_filtered <- roi_notes[!is.na(roi_notes$continuum), ]
  roi_notes_sorted <- roi_notes_filtered[order(roi_notes_filtered$continuum), ]
  ordered_cells <- roi_notes_sorted$ROI
  data_matrix_log2 <- log2(data_matrix + 1)
  obj <- CreateSeuratObject(counts = data_matrix_log2, project = "SeuratProject")
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = data_matrix_log2)
  meta <- obj@meta.data
  meta$slice <- roi_notes$slice[match(rownames(meta), roi_notes$ROI)]
  obj@meta.data <- meta
  obj <- subset(obj, cells = ordered_cells)
  avg_exp <- rowMeans(GetAssayData(obj, layer = "counts"))
  filtered_genes <- names(avg_exp)[avg_exp > 4]
  obj <- subset(obj, features = filtered_genes)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500)
  var_features_info <- HVFInfo(obj)
  
  # Use all 1500 variable genes for downstream analyses
  obj_1500 <- subset(obj, features = VariableFeatures(obj))
  obj_1500 <- ScaleData(obj_1500, vars.to.regress = "slice")
  scaled_data_1500 <- GetAssayData(obj_1500, layer = "scale.data")
  scaled_data_1500 <- scaled_data_1500[, ordered_cells]
  assign(paste0("scaled_data_", tolower(dataset_name)), scaled_data_1500, envir = .GlobalEnv)
  assign(paste0("ordered_cells_", tolower(dataset_name)), ordered_cells, envir = .GlobalEnv)
  
  # Subset to top 50 genes for heatmap visualization only
  top_var_genes <- rownames(var_features_info)[order(-var_features_info$variance)][1:50]
  obj_50 <- subset(obj, features = top_var_genes)
  obj_50 <- ScaleData(obj_50, vars.to.regress = "slice")
  scaled_data_50 <- GetAssayData(obj_50, layer = "scale.data")
  scaled_data_50 <- scaled_data_50[, ordered_cells]
  heatmap_data <- scaled_data_50
  heatmap_data[heatmap_data > 3] <- 3
  heatmap_data[heatmap_data < -3] <- -3
  heatmap_colors <- viridis(100, option = "plasma", direction = 1)
  
  # Nature-style heatmap
  pdf(output_path, width = 14, height = 12)
  pheatmap(heatmap_data,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = title,
           scale = "none",
           breaks = seq(-3, 3, length.out = 101),
           color = heatmap_colors,
           fontsize = 40,
           fontsize_row = 12,
           fontsize_col = 10,
           fontsize_number = 10,
           border_color = "white",
           cellwidth = 6,
           cellheight = 12,
           treeheight_row = 80,
           legend_breaks = seq(-3, 3, by = 2),
           legend_labels = seq(-3, 3, by = 2),
           fontsize_legend = 12,
           angle_col = 0,
           lwd = 0,
           margin = c(100, 100, 100, 200)
  )
  dev.off()
  message("Heatmap generated: ", output_path)
}

# Function to perform sine fitting and enrichment analysis
perform_sine_fit_and_enrichment <- function(scaled_data, ordered_cells, dataset_name, output_dir, set_ids, descriptions_map, legend_map) {
  continuum_values <- roi_notes$continuum[match(ordered_cells, roi_notes$ROI)]
  TOTAL_CONTINUUM <- 100
  omega <- 2 * pi / TOTAL_CONTINUUM
  
  # Parallel sine fitting
  peak_continuum <- foreach(i = 1:nrow(scaled_data), .combine = c) %dopar% {
    expr <- scaled_data[i, ]
    fit <- lm(expr ~ sin(omega * continuum_values) + cos(omega * continuum_values))
    coef <- coef(fit)
    beta_sin <- coef[2]
    beta_cos <- coef[3]
    phi <- atan2(beta_cos, beta_sin)
    x0 <- (pi/2 - phi) / omega
    x_peak <- (x0 %% TOTAL_CONTINUUM + TOTAL_CONTINUUM) %% TOTAL_CONTINUUM
    x_peak
  }
  
  df_phase <- data.frame(Gene = rownames(scaled_data), PeakContinuum = peak_continuum)
  df_phase <- df_phase[order(df_phase$PeakContinuum), ]
  
  # Preload gene mappings to avoid database access in parallel
  all_genes <- unique(df_phase$Gene)
  mapping <- suppressMessages(bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  message(paste("Mapped", sum(!is.na(mapping$ENTREZID)), "out of", length(all_genes), "genes to Entrez IDs"))
  
  # Sliding window enrichment parameters
  window_width <- 33
  step_size <- 3
  window_starts <- seq(0, TOTAL_CONTINUUM, by = step_size)
  
  # Circular range helper function
  in_circular_range <- function(x, start_val, width, total_len) {
    shift <- (x - start_val) %% total_len
    return(shift < width)
  }
  
  # Enrichment and p.adjust extraction function
  enrich_and_extract_padj <- function(entrez_vec, set_ids) {
    if (length(entrez_vec) < 5) {
      return(rep(NA, length(set_ids)))
    }
    
    # KEGG enrichment
    ekegg <- enrichKEGG(gene = entrez_vec, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH")
    kegg_df <- as.data.frame(ekegg)
    
    # GO enrichment using enrichGO
    egobp <- enrichGO(gene = entrez_vec, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05)
    gobp_df <- as.data.frame(egobp)
    
    # Reactome enrichment
    epath <- enrichPathway(gene = entrez_vec, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
    path_df <- as.data.frame(epath)
    
    # Extract p.adjust values
    out_padj <- numeric(length(set_ids))
    for (i in seq_along(set_ids)) {
      set_id <- set_ids[i]
      if (grepl("^GO:", set_id)) {
        row_idx <- which(gobp_df$ID == set_id)
        out_padj[i] <- if (length(row_idx) == 1) gobp_df$p.adjust[row_idx] else NA
      } else if (grepl("^hsa", set_id)) {
        row_idx <- which(kegg_df$ID == set_id)
        out_padj[i] <- if (length(row_idx) == 1) kegg_df$p.adjust[row_idx] else NA
      } else if (grepl("^R-HSA", set_id)) {
        row_idx <- which(path_df$ID == set_id)
        out_padj[i] <- if (length(row_idx) == 1) path_df$p.adjust[row_idx] else NA
      } else {
        out_padj[i] <- NA
      }
    }
    return(out_padj)
  }
  
  # Parallel sliding window enrichment
  results_list <- foreach(idx = seq_along(window_starts), 
                          .combine = rbind, 
                          .packages = c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"),
                          .export = c("df_phase", "mapping", "set_ids", "TOTAL_CONTINUUM", 
                                      "window_width", "in_circular_range", "enrich_and_extract_padj")) %dopar% {
                                        window_start <- window_starts[idx]
                                        in_window <- in_circular_range(df_phase$PeakContinuum, window_start, window_width, TOTAL_CONTINUUM)
                                        df_win <- df_phase[in_window, ]
                                        gene_symbols <- df_win$Gene
                                        center_val <- (window_start + window_width / 2) %% TOTAL_CONTINUUM
                                        
                                        if (length(gene_symbols) > 0) {
                                          entrez_vec <- mapping$ENTREZID[mapping$SYMBOL %in% gene_symbols]
                                          entrez_vec <- unique(entrez_vec[!is.na(entrez_vec)])
                                          padj_vec <- enrich_and_extract_padj(entrez_vec, set_ids)
                                          neglog_vec <- -log10(padj_vec)
                                          neglog_vec[is.na(neglog_vec)] <- 0
                                        } else {
                                          neglog_vec <- rep(0, length(set_ids))
                                        }
                                        
                                        data.frame(
                                          window_center = rep(center_val, length(set_ids)),
                                          set_id = set_ids,
                                          neglog10 = neglog_vec
                                        )
                                      }
  
  # Combine results
  df_results <- results_list
  df_results$description <- descriptions_map[df_results$set_id]
  df_results$legend_label <- legend_map[df_results$set_id]
  
  # Create Nature-style plot
  p <- ggplot(df_results, aes(x = window_center, y = neglog10, color = legend_label)) +
    geom_line(size = 1.2) +
    theme_classic() +
    labs(x = "Continuum (%)", y = expression(bold(-log[10](P[adj]))), title = dataset_name, color = NULL) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      panel.grid = element_blank()
    ) +
    scale_color_manual(values = c("Wnt" = "#1b9e77", "BMP" = "#7570b3", 
                                  "TGFB" = "#e7298a", "Notch" = "#66a61e", "Insulin" = "#e6ab02", 
                                  "PI3K-Akt" = "#a6761d", "Hippo" = "#ff7f00"))
  
  # Save as TIFF
  tiff(file.path(output_dir, paste0(dataset_name, "_top5_pathways_enrichment.tiff")), 
       width = 3.94, height = 2.76, units = "in", res = 300)
  print(p)
  dev.off()
}

# Define set_ids, descriptions_map, and legend_map for NPCs and stroma
npcs_set_ids <- c("hsa04310", "GO:0030509", "hsa04350", "hsa04330", "GO:0008286")
npcs_descriptions_map <- c(
  "hsa04310" = "Wnt signaling pathway",
  "GO:0030509" = "BMP signaling pathway",
  "hsa04350" = "TGF-beta signaling pathway",
  "hsa04330" = "Notch signaling pathway",
  "GO:0008286" = "insulin receptor signaling pathway"
)
npcs_legend_map <- c(
  "hsa04310" = "Wnt",
  "GO:0030509" = "BMP",
  "hsa04350" = "TGFB",
  "hsa04330" = "Notch",
  "GO:0008286" = "Insulin"
)

stroma_set_ids <- c("hsa04310", "GO:0030509", "hsa04350", "hsa04151", "GO:0035329")
stroma_descriptions_map <- c(
  "hsa04310" = "Wnt signaling pathway",
  "GO:0030509" = "BMP signaling pathway",
  "hsa04350" = "TGF-beta signaling pathway",
  "hsa04151" = "PI3K-Akt signaling pathway",
  "GO:0035329" = "hippo signaling"
)
stroma_legend_map <- c(
  "hsa04310" = "Wnt",
  "GO:0030509" = "BMP",
  "hsa04350" = "TGFB",
  "hsa04151" = "PI3K-Akt",
  "GO:0035329" = "Hippo"
)

# Create heatmaps
create_heatmap(npcs_data_path, file.path(output_dir, "NPCs_top50_variable_genes_heatmap.pdf"), "NPCs", "NPCs")
create_heatmap(stroma_data_path, file.path(output_dir, "Stroma_top50_variable_genes_heatmap.pdf"), "Stroma", "Stroma")

# Perform sine fit and enrichment analysis, saving TIFF files
perform_sine_fit_and_enrichment(scaled_data_npcs, ordered_cells_npcs, "NPCs", output_dir, npcs_set_ids, npcs_descriptions_map, npcs_legend_map)
perform_sine_fit_and_enrichment(scaled_data_stroma, ordered_cells_stroma, "Stroma", output_dir, stroma_set_ids, stroma_descriptions_map, stroma_legend_map)

# Clean up parallel backend
stopImplicitCluster()

message("Processing complete. TIFF files saved in ", output_dir)