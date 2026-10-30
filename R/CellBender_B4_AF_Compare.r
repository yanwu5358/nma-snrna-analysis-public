# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)

# -----------------------------------------------------------------------------
# 1. Load the H5 Matrices
# -----------------------------------------------------------------------------
raw_path <- "./raw_feature_bc_matrix.h5"
cb_path  <- "./cellbender_out_filtered.h5"

# Read in the count matrices
raw_counts <- Read10X_h5(raw_path)
cb_counts  <- Read10X_h5(cb_path)

# -----------------------------------------------------------------------------
# 2. Create a New Seurat Object & Add Assays
# -----------------------------------------------------------------------------
# CellBender removes empty droplets, so we only want the cells that survived
valid_cells <- colnames(cb_counts)
raw_counts_sub <- raw_counts[, valid_cells]

# Create a brand new Seurat object using the cleaned CellBender counts as the base
sc_obj <- CreateSeuratObject(counts = cb_counts, project = "example_sample", assay = "CellBender")
sc_obj<-CreateSeuratObject(counts = cb_mtx, project = opts$sample, assay = "CellBender")

# Add the uncorrected raw counts as a secondary assay for comparison
sc_obj[["Raw"]] <- CreateAssayObject(counts = raw_counts_sub)

# -----------------------------------------------------------------------------
# 3. Standard Preprocessing (Required to generate UMAP and Clusters)
# -----------------------------------------------------------------------------
DefaultAssay(sc_obj) <- "CellBender"

# Standard workflow on the cleaned data
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))

# Generate UMAP and Clusters (adjust dims and resolution as needed for your specific tissue)
sc_obj <- RunUMAP(sc_obj, dims = 1:20)
sc_obj <- FindNeighbors(sc_obj, dims = 1:20)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)

# Normalize the Raw assay so we can compare expression levels accurately later
sc_obj <- NormalizeData(sc_obj, assay = "Raw")

# -----------------------------------------------------------------------------
# 4. Ambient Contamination Profile Before & After 
# -----------------------------------------------------------------------------
# Identify genes you suspect are ambient contaminants (e.g., highly expressed muscle/neural genes bleeding into other clusters)
ambient_genes <- c("MYOG", "TTN") # Replace with your specific ambient culprits
ambient_genes <-ambient_genes_present
comparison_plots <- list()

for (gene in ambient_genes) {

  # Plot After (CellBender Corrected)
  DefaultAssay(sc_obj) <- "CellBender"
  expr_values <- GetAssayData(sc_obj, layer = "data")[gene, ]
  # Common cutoffs used in FeaturePlot
  q95 <- quantile(expr_values, 0.95)  # 95th percentile (default-like) 
  p_after <- FeaturePlot(sc_obj, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(After CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
 

    # Plot Before (Raw Counts)
  DefaultAssay(sc_obj) <- "Raw"
  p_before <- FeaturePlot(sc_obj, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(Raw / Before CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  # Combine using patchwork
  comparison_plots[[gene]] <- p_before | p_after
}

# Print the comparison plots
pdf_file <- file.path(opts$out, paste0(opts$sample, "_cellbenderB4_AF_plots.pdf"))
pdf(pdf_file, width = 14, height = 50)

print(wrap_plots(comparison_plots, ncol = 1))
dev.off()


comparison_plots <- list()

for (gene in ambient_genes_present) {

  # Plot After (CellBender Corrected)
  DefaultAssay(seu) <- "RNA"
  expr_values <- GetAssayData(seu, layer = "data")[gene, ]
  # Common cutoffs used in FeaturePlot
  q95 <- quantile(expr_values, 0.95)  # 95th percentile (default-like) 
  p_after <- FeaturePlot(seu, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(After CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
 

    # Plot Before (Raw Counts)
  DefaultAssay(seu) <- "Raw"
  p_before <- FeaturePlot(seu, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(Raw / Before CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  # Combine using patchwork
  comparison_plots[[gene]] <- p_before | p_after
}

# Print the comparison plots
pdf_file <- file.path(opts$out, paste0(opts$sample, "_DF_cellbenderB4_AF_plots.pdf"))
pdf(pdf_file, width = 14, height = 50)

print(wrap_plots(comparison_plots, ncol = 1))
dev.off()


### Plot the feature before and after cellbender/and scDblFinder:
valid_cells <- colnames(seu)
raw_counts_sub <- raw_mtx[, valid_cells]

# Add the uncorrected raw counts as a secondary assay for comparison
seu[["Raw"]] <- CreateAssayObject(counts = raw_counts_sub)





marker_genes <- c("MYT1L","MAPT","SYN3","NYAP2","PTPRR","SOX2","EDNRB","HOPX","ALDH1L1","METRN","SOX10","INSC","CDH19","MPZ","CPQ","PRRX1","COL12A1",  "MEOX2", "PDGFRB","MYOG",  "TNNT3", "TRIM55","RYR1", "MYPN")
comparison_plots <- list()

for (gene in marker_genes) {

  # Plot After (CellBender Corrected)
  DefaultAssay(seu) <- "RNA"
  expr_values <- GetAssayData(seu, layer = "data")[gene, ]
  # Common cutoffs used in FeaturePlot
  q95 <- quantile(expr_values, 0.95)  # 95th percentile (default-like) 
  p_after <- FeaturePlot(seu, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(After CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
 

    # Plot Before (Raw Counts)
  DefaultAssay(seu) <- "Raw"
  p_before <- FeaturePlot(seu, features = gene, reduction = "umap", order = TRUE,max.cutoff = q95) +
    ggtitle(paste(gene, "(Raw / Before CellBender)")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  # Combine using patchwork
  comparison_plots[[gene]] <- p_before | p_after
}

# Print the comparison plots
pdf_file <- file.path(opts$out, paste0(opts$sample, "AnnotationMarkers_cellbenderB4_AF_plots.pdf"))
pdf(pdf_file, width = 10, height = 80)

print(wrap_plots(comparison_plots, ncol = 1))
dev.off()

png_file <- file.path(opts$out, paste0(opts$sample, "AnnotationMarkers_cellbenderB4_AF_plots.png"))
png(png_file, width = 10, height = 80,units="in",res=300)

print(wrap_plots(comparison_plots, ncol = 1))
dev.off()




# -----------------------------------------------------------------------------
# 5. Marker Visualizations (DotPlot & UMAP)
# -----------------------------------------------------------------------------
# Define your markers
markers <- c(
  "PAX7", "MYOD1", "MYOG", "TTN", "DES", "MSC", "ACTA1", "ACTN2", "CHRNB1", "CHRND", "MUSK", "LRP4", # Muscle
  "PAX6", "ELAVL4", "ELAVL3", "STMN2", "NSG2", "SYT1",                                             # Neural
  "PDGFRA", "PDGFRB", "GFAP"                                                                       # Other
)

# Ensure we are plotting using the cleaned data
DefaultAssay(sc_obj) <- "CellBender"

# DotPlot
dot_plot <- DotPlot(sc_obj, features = markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(title = "Cluster Specific Markers (CellBender Corrected)", x = "Marker Gene", y = "Cluster") +
  coord_flip() 

print(dot_plot)

DefaultAssay(seu)<-"RNA"
dot_plot <- DotPlot(seu, features = markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_color_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
                        midpoint = 2, limits = c(0, 4)) +
  labs(title = "Cluster Specific Markers (CellBender Corrected)", x = "Marker Gene", y = "Cluster") +
  coord_flip() 
pdf_file <- file.path(opts$out, paste0(opts$sample, "_Markers_DotPlotAFCellbender.pdf"))
pdf(pdf_file, width = 10, height = 10)
print(dot_plot)
dev.off()

DefaultAssay(seu)<-"Raw"
dot_plot_b4 <- DotPlot(seu, features = markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_color_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
                        midpoint = 2, limits = c(0, 4)) +
  labs(title = "Cluster Specific Markers (B4 CellBender)", x = "Marker Gene", y = "Cluster") +
  coord_flip() 
pdf_file <- file.path(opts$out, paste0(opts$sample, "_Markers_DotPlot_B4Cellbender.pdf"))
pdf(pdf_file, width = 10, height = 10)
print(dot_plot_b4)
dev.off()

# UMAP FeaturePlots
# Saving to PDF to handle the large number of genes
pdf("CellBender_Marker_UMAPs_From_H5.pdf", width = 16, height = 12)
FeaturePlot(sc_obj, features = markers, reduction = "umap", ncol = 5, order = TRUE)
dev.off()