# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         

# setwd("~/gene expression/analysis")

# load libraries
library(Seurat)
library(tidyverse)

# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = 'data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m) # let's take a look at the modalities present

# there are 3 modalities however, we are only interested in the gene expression modality

# to get the counts for the gene expression

cts <-  nsclc.sparse.m$`Gene Expression`
cts[1:10,1:10] # taking a look at the first 10 rows and 10 columns

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)

# this helps us to read our counts into a seurat object. We then name our project as NSCLC
# we then keep all the features that have minimum 3 cells and we want to keep all cells that have at least 200 features

str(nsclc.seurat.obj) # viewing the seurat object
nsclc.seurat.obj
# 29552 features across 42081 samples


# The next step is quality control so that we can filter out the low quality cells
# 1. QC -------

# It is important to look at the number of features or genes in the cell and also 
# the total number of molecules and counts of molecules in a cell
# this gives us an idea of whether the cell is of low quality since low quality of cells will express few genes
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-") # this means we are calculating the percent mitochondrial for all the genes starting wIth MT
# now we have calculated our percentage of mitochondrial and saved in a new column called percent.mt
View(nsclc.seurat.obj@meta.data)

# percent.mt refers to percentage mitochondrial genes which normally contaminate our genes especially in low quality cells
install.packages("SeuratObject")
install.packages("Seurat")
library(Seurat)
library(tidyverse)

pdf('Figure 1 violin_plot.pdf', width = 10, height = 8)
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

dev.off()

# a good quality of cells should follow the linear line drawn in the plot

# the FeatureScatter function allows us to plot two matrix on different axis
# and so we plotted the number of count on x axis and number of genes on y axis
# technically good quality cells should not only have good number of genes detected
# but also good number of molecules detected


# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
# we need filter out doublets and multiple cells which have been sequenced together as single cells so they don't 
# hamper our downstream analysis

# In order to compare the gene expression across multiple cells we have to normalize our data

# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# we want to select or view a few features that exhibit high cell to cell variation
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels

pdf('Figure 2 variablefeatures.pdf', width = 10, height = 8)
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

dev.off()

# there could be unwanted sources of variation in our data and this could be due biological sources or technical sources (such as batch effects)
# eg. difference in cell cycle of these cells and so we want to account for these sources of variations so that during downstream analysis our cells do not cluster
# due to these variations but rather cluster due to biological similarity or biological effects

# 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj) # meaning we want to use all our genes
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
# we perform principal component analysis (PCA) to identify the sources of heterogeneity in our dataset

nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results

print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5) 

# we just want to see the top 5 features
# showing in our principal components


DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
# we are going to choose only the statistically significant principal components
ElbowPlot(nsclc.seurat.obj)


# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)
pdf("Figure 3 Cell cluster.pdf")
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
dev.off()
# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

pdf("Figure 3 Cell cluster_umap.pdf")
DimPlot(nsclc.seurat.obj, reduction = "umap")
dev.off()