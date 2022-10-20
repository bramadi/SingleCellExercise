setwd("~/PBMC")

#Initialize the Seurat object
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Create a Seurat object using the count matrix
# Load the PBMC dataset
pbmc.data <- Read10X_h5("/home/arya/PBMC/data/pbmc_10k_v3_raw_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)
pbmc

#pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells in the control group
head(pbmc@meta.data, 5)

#The QC metrics nFeature_RNA, nCount_RNA, and percent.mt can be visualize:
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

#Demonstrate the link between nCount RNA, nFeature RNA, and percent.mt using dot plots. 
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2

#Eliminate cells with unique feature counts (genes) greater than 5,000 or fewer than 200.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

#Visualize QC metrics again after filtering cells
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

#Use dot plots to show the relationship between nCount_RNA, nFeature_RNA and percent.mt.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2

#nFeature RNA and percent.mo can have different cutoff values.
temp <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
VlnPlot(temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

plot1 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2

rm(temp)

#Normalization of data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#Take 10,000 read counts from the large gene expression matrix in order to visualize the gene expression distribution separately before and after normalization (zeros are not included).
# set seed and put two plots in one figure
set.seed(123)
par(mfrow=c(1,2))
# original expression distribution
raw_geneExp = as.vector(pbmc[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)
# expression distribution after normalization
logNorm_geneExp = as.vector(pbmc[['RNA']]@data) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp)

# same as the above command, skip here
pbmc <- NormalizeData(pbmc)

#Identifying highly variable characteristics (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
 
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
 
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
plot1 + plot2

#The scaling of data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose = FALSE)

#If we wish to scale on the 2000 previously detected variable features, we can omit the features argument from the preceding function call, i.e.
# use the above one, skip here
pbmc <- ScaleData(pbmc

# skip here
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

#Seurat provides several useful ways of visualizing both cells and features that define the PCA.
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Implement non-linear dimension reduction (UMAP/tSNE).
pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)

#The UMAP plot of the single cell clustering results can then be obtained.
DimPlot(pbmc, reduction = "umap")

#Additionally, we can view it using tSNE plot
pbmc <- RunTSNE(pbmc, dims = 1:20, verbose = FALSE)
DimPlot(pbmc, reduction = "tsne")

#We can label individual clusters by setting label = TRUE or by using the LabelClusters method.
DimPlot(pbmc, reduction = "umap", label = TRUE)

plot <- DimPlot(object = pbmc)
LabelClusters(plot = plot, id = 'ident')
