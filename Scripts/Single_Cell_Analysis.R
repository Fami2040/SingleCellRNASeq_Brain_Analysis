# Load the required libraries for data manipulation and visualization
library(Seurat)    # For single-cell RNA-seq analysis
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation

# Step 1: Set file paths and create output folder
# Define file paths for input datasets and metadata
neurons_file <- "C:/Users/FAMI/Downloads/SingleCellRNASeq_Brain_Analysis/Data/Brain_Neurons-counts/Brain_Neurons-counts.csv"
microglia_file <- "C:/Users/FAMI/Downloads/SingleCellRNASeq_Brain_Analysis/Data/Brain_Microglia-counts/Brain_Microglia-counts.csv"
metadata_file <- "C:/Users/FAMI/Downloads/SingleCellRNASeq_Brain_Analysis/Data/metadata.csv"

# Define folder for saving output results
output_folder <- "C:/Users/FAMI/Downloads/SingleCellRNASeq_Brain_Analysis/Results"

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)  # Ensures all nested folders are created
}

# Step 2: Load datasets
# Load count matrices for neurons and microglia
neurons_counts <- read.csv(neurons_file, row.names = 1)
microglia_counts <- read.csv(microglia_file, row.names = 1)

# Step 3: Create Seurat objects
# Create Seurat objects for neurons and microglia
neurons <- CreateSeuratObject(counts = neurons_counts, project = "Brain_Neurons")
microglia <- CreateSeuratObject(counts = microglia_counts, project = "Brain_Microglia")

# Step 4: Merge Seurat objects
# Combine neurons and microglia data into a single Seurat object
brain_combined <- merge(neurons, y = microglia, add.cell.ids = c("Neurons", "Microglia"))

# Step 5: Quality Control
# Calculate percentage of mitochondrial gene content
brain_combined[["percent.mt"]] <- PercentageFeatureSet(brain_combined, pattern = "^MT-")

# Visualize quality control metrics
qc_plot <- VlnPlot(brain_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(output_folder, "quality_control_violin_plot.png"), plot = qc_plot, width = 8, height = 6)

# Filter cells based on quality thresholds
brain_combined <- subset(brain_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Step 6: Normalize, Identify Variable Features, and Scale Data
# Normalize the data
brain_combined <- NormalizeData(brain_combined)

# Identify the most variable genes
brain_combined <- FindVariableFeatures(brain_combined, selection.method = "vst", nfeatures = 2000)

# Scale the data for PCA
brain_combined <- ScaleData(brain_combined)

# Step 7: Run PCA
# Perform PCA to reduce dimensionality
brain_combined <- RunPCA(brain_combined, features = VariableFeatures(object = brain_combined))

# Elbow plot to determine optimal number of PCs
pca_plot <- ElbowPlot(brain_combined)
ggsave(file.path(output_folder, "pca_elbow_plot.png"), plot = pca_plot, width = 8, height = 6)

# Step 8: Clustering and UMAP
# Perform clustering and visualize using UMAP
brain_combined <- FindNeighbors(brain_combined, dims = 1:10)
brain_combined <- FindClusters(brain_combined, resolution = 0.5)
brain_combined <- RunUMAP(brain_combined, dims = 1:10)

# Save UMAP plot with clusters
umap_plot <- DimPlot(brain_combined, reduction = "umap", label = TRUE)
ggsave(file.path(output_folder, "umap_clusters.png"), plot = umap_plot, width = 8, height = 6)

# Step 9: Join Data Layers
# (Optional Step: Custom function, ensure JoinLayers is defined correctly in your setup)
brain_combined <- JoinLayers(brain_combined)

# Step 10: Add Metadata
# Load external metadata and add it to Seurat object
metadata <- read.csv(metadata_file, row.names = 1)
brain_combined <- AddMetaData(brain_combined, metadata = metadata)

# Step 11: UMAP by Mouse Sex
# Visualize clustering by mouse sex
sex_umap_plot <- DimPlot(brain_combined, reduction = "umap", group.by = "mouse.sex", label = TRUE) +
  ggtitle("UMAP Clustering by Mouse Sex")
ggsave(file.path(output_folder, "umap_by_mouse_sex.png"), plot = sex_umap_plot, width = 8, height = 6)

# Step 12: Bar Plot of Cluster Distribution by Mouse Sex
# Group data by clusters and sex, then visualize cluster proportions
sex_distribution <- brain_combined@meta.data %>%
  group_by(seurat_clusters, mouse.sex) %>%
  summarize(count = n())

sex_bar_plot <- ggplot(sex_distribution, aes(x = seurat_clusters, y = count, fill = mouse.sex)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cluster", y = "Proportion", fill = "Mouse Sex") +
  ggtitle("Cluster Distribution by Mouse Sex")
ggsave(file.path(output_folder, "cluster_distribution_by_sex.png"), plot = sex_bar_plot, width = 8, height = 6)

# Step 13: Assign Identities and Find Differential Markers
# Assign new identities based on 'subtissue' metadata
Idents(brain_combined) <- brain_combined$subtissue

# Export identity table
identity_table <- table(Idents(brain_combined))
write.csv(identity_table, file.path(output_folder, "identity_table.csv"))

# Find markers distinguishing Microglia and Neurons
brain_combined$cell_type <- ifelse(grepl("Microglia", colnames(brain_combined)), "Microglia", "Neurons")
Idents(brain_combined) <- brain_combined$cell_type
markers <- FindMarkers(brain_combined, ident.1 = "Microglia", ident.2 = "Neurons")

# Step 14: Heatmap of Top 10 Markers
# Select top 10 markers for heatmap
top10_markers <- markers %>% top_n(10, avg_log2FC)

# Scale data for selected markers and plot heatmap
brain_combined <- ScaleData(brain_combined, features = rownames(top10_markers))
heatmap_plot <- DoHeatmap(brain_combined, features = rownames(top10_markers)) +
  ggtitle("Top 10 Differentially Expressed Genes")
ggsave(file.path(output_folder, "heatmap_top10_genes.png"), plot = heatmap_plot, width = 8, height = 6)

