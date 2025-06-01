# ==========================
# xenium_polar_plot_centered.R
# Generates a filled polar plot for average Z-scores of proliferation genes
# ==========================

# Clear workspace
rm(list = ls())

# Load libraries
library(dplyr)
library(ggplot2)
library(Seurat)

# Define file paths
data_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/combined_gene_counts_matrix_normalized_cap.csv"
roi_path <- "D:/Igor/Xenium/Data for clustering algorithm/Z-Filtered Datasets/Alex_20/ROI_notes.csv"
output_dir <- "D:/Igor/Figures/Section4"

# Set working directory
setwd(output_dir)

# Load data
data_matrix <- as.matrix(read.csv(data_path, row.names = 1))
roi_notes <- read.csv(roi_path)

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

# Define proliferation genes of interest
genes_of_interest <- c("HIST1H1A", "TYMS", "CENPU", "HIST1H4C", "H2AFX","MCM7","MCM6","HELLS","MCM4","PCNA","MCM2","CHEK1") # "HIST1H1A", "TYMS", "CENPU", "TOP2A" are from enriched in 1, others are from heatmap

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
min_radius <- -1  # Minimum radius based on axis limits
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

# [Previous code unchanged until the plotting section]

# Create filled polar plot with solid fill, no legend, and Nature-ified aesthetics
p5 <- ggplot(df_plot, aes(x = continuum, y = r)) +
  geom_polygon(fill = "#FFFF00", color = NA) +  # Solid blue fill, no outline for polygon
  geom_path(data = df_curve, aes(x = continuum, y = r), 
            color = "black", size = 2, inherit.aes = FALSE) +  # Thicker black border for outer points
  coord_polar(theta = "x", clip = "off") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100"),
    limits = c(0, 100)  # Ensure full 0-100 range
  ) +
  scale_y_continuous(
    breaks = c(-.5, 0, .5),
    labels = c("-0.5", "0", "0.5"),
    limits = c(-1, 1)
  ) +
  labs(
    title = "NPC - Proliferation",
    x = "Continuum (%)",
    y = "Gene Set (Z-score)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans"),
    axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 28, face = "bold", margin = margin(r = 15)),
    axis.text = element_text(size = 24, face = "bold", color = "black"),
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1.5, color = "black"),
    plot.margin = margin(20, 20, 20, 20),  # Plenty of padding
    legend.position = "none"  # Remove legend
  )

# Save plot as high-resolution TIFF for Nature
ggsave("npc_prolif.tiff", plot = p5, width = 8, height = 8, dpi = 300, device = "tiff",
       compression = "lzw")

message("Nature-ified filled polar plot generated: ", file.path(output_dir, "npc_prolif.tiff"))