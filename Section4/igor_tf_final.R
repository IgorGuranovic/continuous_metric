# Clear workspace
rm(list = ls())

# Define file paths
npcs_data_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_cap.csv"
stroma_data_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_stroma.csv"
roi_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/ROI_notes.csv"
trrust_path <- "D:/Igor/TF/trrust_rawdata.human.tsv"
npcs_heatmap_output <- "D:/Igor/top50_variable_genes_cap_heatmap.pdf"
stroma_heatmap_output <- "D:/Igor/top50_variable_genes_stroma_heatmap.pdf"

# Load libraries
library(Seurat)
library(pheatmap)
library(viridis)
library(tidyr)
library(ggplot2)

# Function to create heatmap
create_heatmap <- function(data_path, output_path, title, dataset_name) {
  data_matrix <- as.matrix(read.csv(data_path, row.names = 1))
  roi_notes <- read.csv(roi_path)
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
  pdf(output_path, width = 14, height = 12)  # Increased size for dendrogram and colorbar
  pheatmap(heatmap_data,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = title,
           scale = "none",
           breaks = seq(-3, 3, length.out = 101),
           color = heatmap_colors,
           fontsize = 40,           # Balanced title size
           fontsize_row = 12,
           fontsize_col = 10,
           fontsize_number = 10,
           border_color = "white",  # White borders for contrast
           cellwidth = 6,
           cellheight = 12,
           treeheight_row = 80,     # Increased for full dendrogram visibility
           legend_breaks = seq(-3, 3, by = 2),  # More breaks for larger colorbar
           legend_labels = seq(-3, 3, by = 2),
           fontsize_legend = 12,    # Larger legend font
           angle_col = 0,
           lwd = 0,                # Re-enabled thin borders
           margin = c(100, 100, 100, 200)  # Increased right margin for colorbar
  )
  dev.off()
  message("Heatmap generated: ", output_path)
}


# Create heatmaps
create_heatmap(npcs_data_path, npcs_heatmap_output, "NPCs", "NPCs")
create_heatmap(stroma_data_path, stroma_heatmap_output, "Stroma", "Stroma")
message("Heatmap generation complete.")

# --- Transcription Factor Enrichment Analysis ---

# Load TRRUST database
trrust_db <- read.table(trrust_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(trrust_db) <- c("TF", "Target", "Direction", "PMID")

# Define TFs
individual_tfs <- c("PAX2", "WT1", "GATA3", "LEF1", "SOX9", "MYC", "TEAD1")
tfap1_family <- c("FOSL1", "FOSL2", "FOSB", "FOS", "JUN", "JUNB", "JUND")
e2f_family <- c("E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8")
selected_tfs <- c(individual_tfs, tfap1_family, e2f_family)

# TFEA function (Enhanced Nature-quality figures)
perform_tfea <- function(scaled_data, ordered_cells, dataset_name, output_dir) {
  tf_ea_output <- file.path(output_dir, paste0(dataset_name, "_tf_ea_scores.csv"))
  singles_plot_output <- file.path(output_dir, paste0(dataset_name, "_tf_ea_plots_singles.pdf"))
  tfap1_plot_output <- file.path(output_dir, paste0(dataset_name, "_tf_ea_plots_tfap1.pdf"))
  e2f_plot_output <- file.path(output_dir, paste0(dataset_name, "_tf_ea_plots_e2f.pdf"))
  ind_cor_output <- file.path(output_dir, paste0(dataset_name, "_individual_tf_cor_heatmap.pdf"))
  tfap1_cor_output <- file.path(output_dir, paste0(dataset_name, "_tfap1_cor_heatmap.pdf"))
  e2f_cor_output <- file.path(output_dir, paste0(dataset_name, "_e2f_cor_heatmap.pdf"))
  
  # TFEA with sliding windows (using 1500 genes)
  window_width <- 25
  step_size <- 3
  total_continuum <- 100
  window_starts <- seq(0, total_continuum, by = step_size)
  roi_notes <- read.csv(roi_path)
  roi_notes_filtered <- roi_notes[!is.na(roi_notes$continuum), ]
  roi_notes_sorted <- roi_notes_filtered[order(roi_notes_filtered$continuum), ]
  actual_continuum_values <- roi_notes_sorted$continuum
  
  results_tf_list <- list()
  for (idx in seq_along(window_starts)) {
    window_start <- window_starts[idx]
    window_end <- window_start + window_width
    cells_in_window <- if (window_end > total_continuum) {
      which((actual_continuum_values >= window_start & actual_continuum_values <= total_continuum) |
              (actual_continuum_values >= 0 & actual_continuum_values < (window_end %% total_continuum)))
    } else {
      which(actual_continuum_values >= window_start & actual_continuum_values < window_end)
    }
    if (length(cells_in_window) == 0) {
      message("Skipping window ", window_start, "-", window_end, ": no cells.")
      next
    }
    
    window_data <- scaled_data[, cells_in_window, drop = FALSE]
    mean_expr <- rowMeans(window_data)
    gene_ranks <- rank(-mean_expr, ties.method = "min")
    genes_in_window <- rownames(scaled_data)
    
    for (tf in selected_tfs) {
      tf_targets <- trrust_db$Target[trrust_db$TF == tf]
      if (length(tf_targets) == 0) {
        message("No targets found for TF ", tf, " in TRRUST.")
        next
      }
      is_target <- genes_in_window %in% tf_targets
      if (sum(is_target) == 0) {
        message("No matching targets for TF ", tf, " in window ", window_start, "-", window_end)
        next
      }
      w_i <- ifelse(is_target, 1, 0)
      w_sorted <- w_i[order(gene_ranks)]
      N <- length(w_sorted)
      e_i <- cumsum(w_sorted) / sum(w_sorted)
      background <- seq(0, 1, length.out = N)
      E <- (2 / N) * sum(e_i - background)
      
      results_tf_list[[length(results_tf_list) + 1]] <- data.frame(
        window_center = (window_start + (window_width / 2)) %% total_continuum,
        TF = tf,
        E_score = E
      )
    }
  }
  if (length(results_tf_list) == 0) {
    stop("No TF E-scores calculated for ", dataset_name, ". Check TRRUST or gene names.")
  }
  df_tf_results <- do.call(rbind, results_tf_list)
  write.csv(df_tf_results, tf_ea_output, row.names = FALSE)
  
  # Plot individual TFs (Enhanced Nature-style)
  ind_data <- df_tf_results[df_tf_results$TF %in% individual_tfs, ]
  if (nrow(ind_data) > 0 && length(unique(ind_data$TF)) > 0) {
    pdf(singles_plot_output, width = 12, height = 8)
    p <- ggplot(ind_data, aes(x = window_center, y = E_score, color = TF)) +
      geom_line(size = 2) +  # Thicker lines
      theme_minimal(base_size = 12) +  # Clean theme
      labs(x = "Continuum (%)", y = "E-score", title = paste(dataset_name, "TFEA: Key TFs")) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +
      scale_color_brewer(palette = "Paired") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 2) +
      theme(
        text = element_text(family = "sans"),
        axis.text = element_text(size = 36, color = "black"),  # Larger tick labels
        axis.title = element_text(size = 48, color = "black", face = "bold"),  # Larger, bold axis titles
        plot.title = element_text(size = 48, hjust = 0.5, face = "bold"),  # Larger, bold plot title
        panel.grid = element_blank(),  # No gridlines
        axis.line = element_line(size = 2, color = "black"),  # Thick axis lines
        legend.position = "top",
        legend.text = element_text(size = 40),  # Larger legend text
        legend.title = element_blank(),
        plot.margin = margin(15, 15, 15, 15)
      ) +
      guides(color = guide_legend(nrow = 2))
    print(p)
    dev.off()
  } else {
    message("Skipping individual TFs plot: no valid data.")
  }
  
  # Plot TFAP1 family (Enhanced Nature-style)
  tfap1_data <- df_tf_results[df_tf_results$TF %in% tfap1_family, ]
  if (nrow(tfap1_data) > 0 && length(unique(tfap1_data$TF)) > 0) {
    pdf(tfap1_plot_output, width = 12, height = 8)
    p <- ggplot(tfap1_data, aes(x = window_center, y = E_score, color = TF)) +
      geom_line(size = 2) +
      theme_minimal(base_size = 12) +
      labs(x = "Continuum (%)", y = "E-score", title = paste(dataset_name, "TFEA: TFAP1 Family")) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +
      scale_color_brewer(palette = "Set1") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 2) +
      theme(
        text = element_text(family = "sans"),
        axis.text = element_text(size = 36, color = "black"),
        axis.title = element_text(size = 48, color = "black", face = "bold"),
        plot.title = element_text(size = 48, hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(size = 2, color = "black"),  # Thick axis lines
        legend.position = "top",
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        plot.margin = margin(15, 15, 15, 15)
      ) +
      guides(color = guide_legend(nrow = 2))
    print(p)
    dev.off()
  } else {
    message("Skipping TFAP1 family plot: no valid data.")
  }
  
  # Plot E2F family (Enhanced Nature-style)
  e2f_data <- df_tf_results[df_tf_results$TF %in% e2f_family, ]
  if (nrow(e2f_data) > 0 && length(unique(e2f_data$TF)) > 0) {
    pdf(e2f_plot_output, width = 12, height = 8)
    p <- ggplot(e2f_data, aes(x = window_center, y = E_score, color = TF)) +
      geom_line(size = 2) +
      theme_minimal(base_size = 12) +
      labs(x = "Continuum (%)", y = "E-score", title = paste(dataset_name, "TFEA: E2F Family")) +
      scale_x_continuous(breaks = seq(0, 100, 20)) +
      scale_color_brewer(palette = "Set2") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 2) +
      theme(
        text = element_text(family = "sans"),
        axis.text = element_text(size = 36, color = "black"),
        axis.title = element_text(size = 48, color = "black", face = "bold"),
        plot.title = element_text(size = 48, hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(size = 2, color = "black"),  # Thick axis lines
        legend.position = "top",
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        plot.margin = margin(15, 15, 15, 15)
      ) +
      guides(color = guide_legend(nrow = 1))
    print(p)
    dev.off()
  } else {
    message("Skipping E2F family plot: no valid data.")
  }
  
  # Correlation heatmaps (Nature-style, unchanged)
  color_breaks <- seq(-1, 1, length.out = 101)
  color_palette <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#FFFFFF", "#F46D43", "#D73027", "#A50026"))(100)
  if (nrow(ind_data) > 0 && length(unique(ind_data$TF)) > 1) {
    ind_wide <- pivot_wider(ind_data, names_from = TF, values_from = E_score, values_fill = 0)
    ind_tfs <- individual_tfs[individual_tfs %in% colnames(ind_wide)]
    ind_cor <- cor(ind_wide[, ind_tfs, drop = FALSE], method = "pearson")
    pdf(ind_cor_output, width = 10, height = 10)
    pheatmap(ind_cor,
             color = color_palette,
             breaks = color_breaks,
             border_color = "black",
             #main = paste(dataset_name, "TF Correlation: Key TFs"),
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize = 12,
             fontsize_row = 12,
             fontsize_col = 12,
             fontsize_number = 10,
             angle_col = 45,
             lwd = 2,
             treeheight_row = 40,
             treeheight_col = 40,
             cellwidth = 35,
             cellheight = 35)
    dev.off()
  } else {
    message("Skipping individual TF correlation heatmap: insufficient data.")
  }
  
  if (nrow(tfap1_data) > 0 && length(unique(tfap1_data$TF)) > 1) {
    tfap1_wide <- pivot_wider(tfap1_data, names_from = TF, values_from = E_score, values_fill = 0)
    tfap1_tfs <- tfap1_family[tfap1_family %in% colnames(tfap1_wide)]
    tfap1_cor <- cor(tfap1_wide[, tfap1_tfs, drop = FALSE], method = "pearson")
    pdf(tfap1_cor_output, width = 8, height = 8)
    pheatmap(tfap1_cor,
             color = color_palette,
             breaks = color_breaks,
             border_color = "black",
             #main = paste(dataset_name, "TF Correlation: TFAP1 Family"),
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize = 12,
             fontsize_row = 12,
             fontsize_col = 12,
             fontsize_number = 10,
             angle_col = 45,
             lwd = 2,
             treeheight_row = 30,
             treeheight_col = 30,
             cellwidth = 45,
             cellheight = 45)
    dev.off()
  } else {
    message("Skipping TFAP1 correlation heatmap: insufficient data.")
  }
  
  if (nrow(e2f_data) > 0 && length(unique(e2f_data$TF)) > 1) {
    e2f_wide <- pivot_wider(e2f_data, names_from = TF, values_from = E_score, values_fill = 0)
    e2f_tfs <- e2f_family[e2f_family %in% colnames(e2f_wide)]
    e2f_cor <- cor(e2f_wide[, e2f_tfs, drop = FALSE], method = "pearson")
    pdf(e2f_cor_output, width = 8, height = 8)
    pheatmap(e2f_cor,
             color = color_palette,
             breaks = color_breaks,
             border_color = "black",
             #main = paste(dataset_name, "TF Correlation: E2F Family"),
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize = 12,
             fontsize_row = 12,
             fontsize_col = 12,
             fontsize_number = 10,
             angle_col = 45,
             lwd = 2,
             treeheight_row = 30,
             treeheight_col = 30,
             cellwidth = 45,
             cellheight = 45)
    dev.off()
  } else {
    message("Skipping E2F correlation heatmap: insufficient data.")
  }
  
  message("TFEA completed for ", dataset_name)
}

# Run TFEA
perform_tfea(scaled_data_npcs, ordered_cells_npcs, "NPCs", "D:/Igor")
perform_tfea(scaled_data_stroma, ordered_cells_stroma, "Stroma", "D:/Igor")
message("TFEA analysis complete for both datasets.")
