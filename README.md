# Single-Cell RNA-seq Analysis of Brain Neurons and Microglia

## 1. Project Description
This project involves the analysis of single-cell RNA sequencing (scRNA-seq) data from brain neurons and microglia. The objective is to identify distinct clusters of cells, explore their biological characteristics, and investigate potential differences based on metadata, such as mouse sex.

## 2. Folder Structure
**Data/**  
Contains input files including count matrices and metadata:
- Brain_Neurons-counts.csv
- Brain_Microglia-counts.csv
- metadata.csv

**Results/**  
Stores output plots and results:
- quality_control_violin_plot.png
- pca_elbow_plot.png
- umap_clusters.png
- umap_by_mouse_sex.png
- cluster_distribution_by_sex.png
- identity_table.csv
- heatmap_top10_genes.png

## 3. Software Requirements
- **R** (version X.X.X or later)
- **RStudio** (optional)

Required R packages:
- Seurat
- ggplot2
- dplyr

## 4. Steps and Workflow

### Step 1: Set File Paths and Create Output Folder
- Define paths for input datasets.
- Create an output folder for saving results.

### Step 2: Load Datasets
- Import count matrices for neurons and microglia.
- Load metadata for mouse-specific information.

### Step 3: Create Seurat Objects
- Construct Seurat objects for neurons and microglia using the count matrices.

### Step 4: Merge Seurat Objects
- Combine the neurons and microglia datasets into a single Seurat object.

### Step 5: Quality Control (QC)
- Calculate mitochondrial gene content.
- Visualize quality metrics using violin plots.
- Filter low-quality cells based on predefined thresholds.

### Step 6: Normalize, Identify Variable Features, and Scale Data
- Normalize the data.
- Identify the top 2000 most variable genes.
- Scale the data for downstream dimensionality reduction.

### Step 7: PCA and Dimensionality Reduction
- Perform PCA.
- Visualize the elbow plot to determine the optimal number of PCs.

### Step 8: Clustering and UMAP
- Perform clustering using the first 10 PCs.
- Visualize clusters using UMAP.

### Step 9: Join Data Layers (Optional)
- Combine various data layers (if applicable).

### Step 10: Add Metadata
- Integrate metadata into the Seurat object.

### Step 11: UMAP by Mouse Sex
- Generate a UMAP plot to explore clustering by mouse sex.

### Step 12: Cluster Distribution by Mouse Sex
- Create a bar plot showing the proportion of cells by cluster and mouse sex.

### Step 13: Assign Identities and Find Differential Markers
- Assign identities based on tissue metadata.
- Export identity table.
- Identify differentially expressed genes between microglia and neurons.

### Step 14: Heatmap of Top 10 Markers
- Generate a heatmap for the top 10 differentially expressed genes.

## 5. Key Outputs

### Plots:
- Quality control violin plot
- PCA elbow plot
- UMAP plots (clusters and sex)
- Cluster distribution bar plot
- Heatmap of top markers

### Tables:
- Identity table
- Top 10 differentially expressed genes

## 6. How to Run the Analysis
1. Clone or download the project repository.
2. Ensure all required R packages are installed.
3. Update file paths in the script to match your system.
4. Run the R script step by step to reproduce the analysis.

## 7. References
- Seurat Documentation: [https://satijalab.org/seurat/](https://satijalab.org/seurat/)
- ggplot2 Documentation: [https://ggplot2.tidyverse.org/](https://ggplot2.tidyverse.org/)
- dplyr Documentation: [https://dplyr.tidyverse.org/](https://dplyr.tidyverse.org/)
