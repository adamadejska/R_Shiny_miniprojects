# Seurat tutorial from satijalab.org

library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/Desktop/Seurat_data/Seurat_app/data")
pbmc.data <- Read10X(data.dir = 'hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Preprocessing workflow: QC and selecting cells for future analysis
# Filter cells that have unique feature counts over 2,500 or less than 200
# Filter cells that have >5 % mitochondrial counts
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

# Visualize feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Normalize the data using a global-scaling normalization "LogNormalize"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features (feature selection)
# Which features are highly expressed in some cells and lowly expressed in others?
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data to shift the mean gene expression to 0, and shift the variance across cells to 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear dimensional reduction
# Run PCA, use previously determined variables features as input
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Explore the primary sources of heterogeneity in our dataset
DimHeatmap(pbmc, dims = 1, cells = 500, balanced=TRUE)

# Determine the dimensionality of the dataset to overcome the extensive technical noise
# Top PCs should represent a robust compression of the dataset
# Determine how many PCs we should include

# Compare the distribution of p-values for each PC with a uniform distribution
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

# Elbow plot: ranking of PCs based on the % of variance explained by each one
ElbowPlot(pbmc)

# Cluster the  cells using K-nearest-neighbor graph
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Run a non-linear dimensional reduction (UMAP)
# To learn the underlying manifold of the data in order to place similar cells together in
# a low-dimensional space
pbmc <- RunUMAP(pbmc, dims=1:10)
DimPlot(pbmc, reduction="umap")

# Save the plot
saveRDS(pbmc, file="pbmc_tutorial_umap.rds")

# Find differentially expressed features 
# Find markers that define clusters via differential expression.
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n=5)

# Find markers that distinguish cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, indent.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n=5)

# Find markers for every cluster compared to all remining cells, only positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>% group_by(cluster) %>% slice_max(n=2, order_by=avg_logFC)

# ROC test - 'classification power; for any individual marker (0=random, 1=perfect)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos =  TRUE)
head(cluster0.markers, n=5)

# Visualize marker expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot="counts", log=TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"))

# Generate an expression heat map for given cells and features
pbmc.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assign cell type identity to clusters
# Match clusters with canonical markers for known cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label=TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file="pbmc3k_final.rds")
