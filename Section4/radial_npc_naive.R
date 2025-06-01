# ==========================
# xenium_polar_plots_combined.R
# Generates filled polar plots for average Z-scores of multiple gene sets
# ==========================

# Clear workspace
rm(list = ls())

# Load libraries
library(dplyr)
library(ggplot2)
library(Seurat)
library(grid)  # For manual grob combination
library(magick)
library(pdftools)

# Install required packages if not already installed
if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick")
}
if (!requireNamespace("pdftools", quietly = TRUE)) {
  install.packages("pdftools")
}

# Define file paths
data_path_npc <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_cap.csv"
data_path_stroma <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_stroma.csv"
roi_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/ROI_notes.csv"
output_dir <- "D:/Igor/Figures/Section4"

# Set working directory
setwd(output_dir)

# Load ROI notes (shared across all plots)
roi_notes <- read.csv(roi_path)

# Define configurations for each plot in a data frame
plot_configs <- data.frame(
  name = c("npc_naive", "npc_diff", "npc_prolif", "stroma_shallow", "stroma_deep"),
  title = c("NPC - Naive", "NPC - Differentiation", "NPC - Proliferation", "Stroma - Shallow", "Stroma - Deep"),
  fill_color = c("#1F77B4", "#FFC0CB", "#FFFF00", "#1F77B4", "#FFC0CB"),
  data_path = c(data_path_npc, data_path_npc, data_path_npc, data_path_stroma, data_path_stroma),
  y_limit_lower = c(-0.7, -0.7, -0.7, -0.7, -0.7),
  y_limit_upper = c(0.7, 0.7, 0.7, 0.7, 0.7),
  stringsAsFactors = FALSE
)

# Define gene sets for each plot
gene_sets <- list(
  npc_naive = c("TCF4", "MEOX1", "NBL1", "ABTB2", "ROBO2", "TRIL", "CHRNA1", "PCDH18", "TCF7L2", "PCDH15", "ELAVL4", "CCDC80", "TTC28", "FOXD1", "TMEM100"),
  npc_diff = c("ID3", "ITPR1", "NOTCH2", "GXYLT2", "PAX2", "ID1", "RXRA", "DNM1", "LYPD1"),
  npc_prolif = c("HIST1H1A", "TYMS", "CENPU", "HIST1H4C", "H2AFX", "MCM7", "MCM6", "HELLS", "MCM4", "PCNA", "MCM2", "CHEK1"),
  stroma_shallow = c("SFRP1", "CDCA7L", "PDE5A", "TGFBI", "NDNF", "FOXD1", "CYP1B1", "TRIL"),
  stroma_deep = c("GATA6", "THBS1", "CXCL12", "ID2", "ID3", "TNC")
)

# Function to generate a plot for a given configuration
generate_polar_plot <- function(config, genes_of_interest) {
  # Load data
  data_matrix <- as.matrix(read.csv(config$data_path, row.names = 1))
  
  # Remove "_count" suffix from matrix cell names
  colnames(data_matrix) <- gsub("_count$", "", colnames(data_matrix))
  
  # Sort ROIs by continuum
  roi_notes_filtered <- roi_notes[!is.na(roi_notes$continuum), ]
  roi_notes_sorted <- roi_notes_filtered[order(roi_notes_filtered$continuum), ]
  ordered_cells <- roi_notes_sorted$ROI
  
  # Log2+1 transformation
  data_matrix_log2 <- log2(data_matrix + 1)
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = data_matrix_log2, project = "SeuratProject")
  obj <- SetAssayData(object = obj, assay = "RNA", layer = "data", new.data = data_matrix_log2)
  
  # Add metadata and filter ROIs
  meta <- obj@meta.data
  meta$slice <- roi_notes$slice[match(rownames(meta), roi_notes$ROI)]
  obj@meta.data <- meta
  obj <- subset(obj, cells = ordered_cells)
  
  # Filter genes with average expression > 4
  avg_exp <- rowMeans(GetAssayData(obj, layer = "counts"))
  filtered_genes <- names(avg_exp)[avg_exp > 4]
  obj <- subset(obj, features = filtered_genes)
  
  # Find top 1500 variable features
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500)
  obj <- subset(obj, features = VariableFeatures(obj))
  
  # Scale data
  obj <- ScaleData(obj, vars.to.regress = "slice")
  
  # Order ROIs by continuum in scaled data
  scaled_data <- GetAssayData(obj, layer = "scale.data")
  scaled_data <- scaled_data[, ordered_cells]
  
  # Calculate average Z-score across genes of interest
  heatmap_subset <- scaled_data[genes_of_interest, , drop = FALSE]
  df_avg_z <- data.frame(
    Selection = colnames(heatmap_subset),
    avg_z = colMeans(heatmap_subset, na.rm = TRUE)
  )
  
  # Merge with ROI notes to get continuum
  df_merged <- df_avg_z %>%
    left_join(roi_notes, by = c("Selection" = "ROI"))
  
  # Fit constrained LOESS to ensure same value at continuum 0 and 100
  df_ordered <- df_merged[order(df_merged$continuum), ]
  n10 <- floor(nrow(df_ordered) * 0.1)
  first10 <- df_ordered$avg_z[1:n10]
  last10 <- df_ordered$avg_z[(nrow(df_ordered) - n10 + 1):nrow(df_ordered)]
  boundaryVal <- mean(c(first10, last10), na.rm = TRUE)
  
  df_boundary <- data.frame(continuum = c(0, 100), avg_z = c(boundaryVal, boundaryVal))
  df_for_loess <- rbind(df_merged[, c("continuum", "avg_z")], df_boundary)
  w <- rep(1, nrow(df_for_loess))
  w[(nrow(df_for_loess) - 1):nrow(df_for_loess)] <- 100
  
  loess_model_constrained <- loess(avg_z ~ continuum, data = df_for_loess, span = 0.75, weights = w)
  
  # Generate smooth data for plotting
  df_smooth <- data.frame(
    continuum = seq(min(df_merged$continuum, na.rm = TRUE), max(df_merged$continuum, na.rm = TRUE), length.out = 200)
  )
  df_smooth$avg_z <- predict(loess_model_constrained, newdata = df_smooth)
  
  # Prepare data for filled polar plot
  df_smooth <- df_smooth %>%
    mutate(
      r = avg_z  # Use avg_z directly as radius, including negative values
    )
  
  df_smooth[nrow(df_smooth), "r"] <- df_smooth[1, "r"]
  df_smooth[nrow(df_smooth), "avg_z"] <- df_smooth[1, "avg_z"]
  
  # Create data for filling to minimum radius
  min_radius <- config$y_limit_lower
  df_curve <- df_smooth %>%
    mutate(type = "curve")
  df_baseline <- df_smooth %>%
    mutate(
      r = min_radius,  # Set radius to minimum
      avg_z = min_radius,  # Set fill value for baseline
      type = "baseline"
    )
  # Combine points: curve forward, baseline backward to enclose area
  df_plot <- rbind(df_curve, df_baseline[nrow(df_baseline):1, ])
  
  # Create filled polar plot with solid fill, no legend, and Nature-ified aesthetics (using Script 1 formatting)
  p_base <- ggplot(df_plot, aes(x = continuum, y = r)) +
    geom_polygon(fill = config$fill_color, color = NA) +  # Solid fill, no outline for polygon
    geom_path(data = df_curve, aes(x = continuum, y = r), 
              color = "black", size = 2, inherit.aes = FALSE) +  # Thicker black border for outer points
    coord_polar(theta = "x", clip = "off", direction = 1) +  # direction = 1 ensures clockwise
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100"),
      limits = c(0, 100)  # Ensure full 0-100 range
    ) +
    scale_y_continuous(
      breaks = c(-0.5, 0, 0.5),
      labels = c("-0.5", "0", "0.5"),
      limits = c(config$y_limit_lower, config$y_limit_upper)
    ) +
    labs(
      title = config$title,
      x = "Continuum (%)",
      y = "Gene Set (Z-score)"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_text(size = 36, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 36, face = "bold", margin = margin(r = 15)),
      axis.text = element_text(size = 32, face = "bold", color = "black"),
      plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "gray80", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20),  # Plenty of padding
      legend.position = "none"  # Remove legend
    )
  
  # Create Cartesian elements (dotted lines and curved arrow) using Script 1 formatting
  y_cartesian <- c(-0.1, 0.2, 0.5)  # Positions for dotted lines
  x_start <- -0.6  # Padding from y-axis labels
  x_end <- 0.2  # Center line (continuum 0/100)
  arrow_y <- -1.2  # Below the plot, above x-axis title
  arrow_x_start <- 0.4  # Right side
  arrow_x_end <- 0  # Left side
  
  cartesian_elements <- ggplot() +
    geom_segment(aes(x = x_start, xend = x_end, y = y_cartesian[1], yend = y_cartesian[1]), 
                 linetype = "dotted", color = "black", size = 1) +
    geom_segment(aes(x = x_start, xend = x_end, y = y_cartesian[2], yend = y_cartesian[2]), 
                 linetype = "dotted", color = "black", size = 1) +
    geom_segment(aes(x = x_start, xend = x_end, y = y_cartesian[3], yend = y_cartesian[3]), 
                 linetype = "dotted", color = "black", size = 1) +
    geom_curve(aes(x = arrow_x_start, xend = arrow_x_end, y = arrow_y, yend = arrow_y), 
               curvature = -0.3, color = "black", size = 1.5,
               arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1.5, 1)) +
    theme_void()
  
  # Convert both plots to grobs
  polar_grob <- ggplotGrob(p_base)
  cartesian_grob <- ggplotGrob(cartesian_elements)
  
  # Create a new grid page and draw the polar plot first
  pdf_file <- file.path(output_dir, paste0(config$name, ".pdf"))
  tiff_file <- file.path(output_dir, paste0(config$name, ".tiff"))
  
  pdf(pdf_file, width = 8, height = 8)
  grid.newpage()
  grid.draw(polar_grob)
  grid.draw(cartesian_grob)
  dev.off()
  
  # Convert PDF to TIFF
  img <- image_read_pdf(pdf_file, density = 300)
  image_write(img, path = tiff_file, format = "tiff", compression = "LZW")
  
  # Output messages
  message("Nature-ified filled polar plot generated as PDF: ", pdf_file)
  message("PDF converted to TIFF with 300 DPI and LZW compression: ", tiff_file)
}

# Loop over each configuration to generate the plots
for (i in 1:nrow(plot_configs)) {
  config <- plot_configs[i, ]
  genes <- gene_sets[[config$name]]
  generate_polar_plot(config, genes)
}