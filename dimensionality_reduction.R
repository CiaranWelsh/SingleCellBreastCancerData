# dimensionality reduction using pca and umap methods.
library(Seurat)

# load global vars
source('global_variables.R')

data = readRDS(DATA_NORMED_AND_INTEGRATED)
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:10)

plot_dim_reduction = function(data, label, red='pca', save = F, dims=1:3, 
                              group.by=NULL, split.by=NULL, plot_elbow = F){
  data = ScaleData(data)
  data = RunPCA(data, features = VariableFeatures(data))
  if (plot_elbow){
    elbow_fname = file.path(UMAPS_DIR, paste0('elbow_', label, max(length(dims)), '.jpeg'))
    p = ElbowPlot(data, ndims = 50)
    if (save == T){
      ggsave(elbow_fname)
    } else{
      print(p)
    }
  }
  data = FindNeighbors(data, dims=dims)
  data = FindClusters(data, resolution = 0.5)
  if (red == 'pca'){
    pca_fname = file.path(UMAPS_DIR, paste0('pca_', label, max(dims), '.jpeg'))
    p = DimPlot(data, reduction = 'pca', group.by = group.by, split.by = split.by)
    if (save == T){
      ggsave(pca_fname)
      message(paste0('saved plot to "', pca_fname, '"'))
    } else {
      print(p)
    }
  } else if( red == 'umap') {
    umap_fname = file.path(UMAPS_DIR, paste0('umap_', label, max(dims), '.jpeg'))
    data = RunUMAP(data, dims = dims, n.epochs = 3000)
    p = DimPlot(data, reduction = 'umap', group.by = group.by, split.by = split.by) 
    if (save == T){
      ggsave(umap_fname)
      message(paste0('saved plot to "', umap_fname, '"'))
    } else {
      print(p)
    }
  }
  return (data)
}
?RunUMAP
for (i in 2:20){
  message(paste('plotting', i))
  data = plot_dim_reduction(
    data, 'normed', 'umap', 
    save = T, dims = 1:i #group.by = 'time'
    )
}

data = plot_dim_reduction(data, 'normed', 'pca', save = T)

data.c1_markers <- FindMarkers(data, ident.1 = 1, ident.2=0, min.pct = 0.25)
deg = rownames(data.c1_markers[data.c1_markers$p_val_adj < 0.001,])
deg
data@meta.data
VlnPlot(data, features = c('BRCA1', deg[1:8]),#sample(deg, 8)), 
        ncol = 3, group.by = 'seurat_clusters', 
        log = T)
RidgePlot(data, features = c('BRCA1', deg[1:8]), 
          group.by = 'cell_line')
FeatureScatter(data, feature1 = "UBE2T", feature2 = "CDK1", 
               group.by = 'cell_line'
               )

?RidgePlot
Think about copying the code for one of their plots and producing a time series plotter

FeaturePlot(data, features = c('BRCA1', sample(deg, 3)))
# is there any structure in the data due to time? 
?VlnPlot
?FindMarkers
library(dplyr)
top10 <- data.c1_markers %>% group_by('cell_line') %>% top_n(n = 10, wt = avg_logFC)

top10
DoHeatmap(data, features = top10$gene) + NoLegend()

# 
# data@
# FeatureLocator()
# 
# ?DimPlot
# colnames(pd.mcf7)
# 
# DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
# 
# VlnPlot(data, features = c('IDO1', 'SMAD3', 'PIK3R1'), group.by = 'cell_line')
# ?VlnPlot

# mcf7.anchors = FindIntegrationAnchors(mcf7.seurat, dims=1:15)
# t47d.anchors = FindIntegrationAnchors(t47d.seurat, dims=1:15)
# mcf7.seurat.subset
# 
# '
# machine learning problem
# ------------------------
# - Estimate dublets, dead cells or otherwise poor quality from number of features
#    number of gene counts and percentage mt genes
# - How to generalize to other datasets?
# - Unsupervised would be best
# 
# QC metrics for single cell 
# --------------------------
# - The number of unique genes detected in each cell (nFeature_RNA). 
#   - When this metric is very low, a droplet may have been empty indicating low quality cell
#   - When very high, may be dublet cell
# - Total number of molecules detected (nCount) in a cell indicates similarly to nFeatures_RNA
# - mitochondrial genome
#   - low-quality or dying cells may have extensive mt geneome contamination
#   - 
# '
# x = mcf7.seurat@meta.data[order(mcf7.seurat@meta.data$percent.mt, decreasing = T),]
# y = t47d.seurat@meta.data[order(t47d.seurat@meta.data$percent.mt, decreasing = T),]
# x$nCount_RNA = x$nCount_RNA / 10000
# x['nFeat*nCount'] = order(x$nFeature_RNA * x$nCount_RNA)
# x = x[order(x$`nFeat*nCount`), ]
# x
# y
# ?order
# x
# x
# hist(x$nFeature_RNA)
# ?hist
# plot(x$nFeature_RNA, x$nCount_RNA)
# # Visualize QC metrics as a violin plot
# mcf7.seurat = subset(mcf7.seurat, percent.mt < 10)
# 
# t47d.seurat = subset(t47d.seurat, percent.mt < 25 
#                      & nFeature_RNA > 10000,
# )
# 
# VlnPlot(mcf7.seurat, 
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1'), 
#         ncol = 4
# )
# 
# VlnPlot(t47d.seurat, 
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1'), 
#         ncol = 4
# )
# 
# 
# y = y[order(y$percent.mt, decreasing = T),]
# y
# ?VlnPlot
# 
# ?VlnPlot
# 
# 
# 
# plot1 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# 
# 
# plot1 <- FeatureScatter(t47d.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(t47d.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# 
# 
# # mcf7.seurat <- subset(mcf7.seurat, subset = nFeature_RNA > 3000)
# 
# 
# mcf7.seurat <- NormalizeData(mcf7.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# t47d.seurat <- NormalizeData(t47d.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# mcf7.seurat
# t47d.seurat
# 
# VlnPlot(mcf7.seurat, 
#         features = c(sample(mcf7_non0_genes, 5), 'IDO1'),
#         ncol = 3
# )
# 
# mcf7.seurat
# 
# VlnPlot(t47d.seurat, 
#         features = c(sample(t47d_non0_genes, 5), 'IDO1'),
#         ncol = 3
# )
# 
# mcf7.seurat <- FindVariableFeatures(mcf7.seurat, selection.method = "vst", nfeatures = 2000)
# mcf7.seurat@reductions
# 
# top10 = head(VariableFeatures(mcf7.seurat))
# 
# top10
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(mcf7.seurat)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
# 
# 
# all.genes <- rownames(mcf7.seurat)
# mcf7.seurat <- ScaleData(mcf7.seurat, features = all.genes)
# mcf7.seurat <- RunPCA(mcf7.seurat, features = VariableFeatures(object = mcf7.seurat))
# 
# print(mcf7.seurat[["pca"]], dims = 1:5, nfeatures = 5)
# 
# VizDimLoadings(mcf7.seurat, dims = 1:2, reduction = "pca")
# 
# 
# DimPlot(mcf7.seurat, reduction = "pca")
# 
# DimHeatmap(mcf7.seurat, dims = 1:2, cells = 500, balanced = TRUE)
# 
# mcf7.seurat <- JackStraw(mcf7.seurat, num.replicate = 100)
# mcf7.seurat <- ScoreJackStraw(mcf7.seurat, dims = 1:20)
# JackStrawPlot(mcf7.seurat, dims = 1:20)
# ElbowPlot(mcf7.seurat)
# 
# 
# mcf7.seurat <- FindNeighbors(mcf7.seurat, dims = 1:10)
# mcf7.seurat <- FindClusters(mcf7.seurat, resolution = 0.5)
# head(Idents(mcf7.seurat), 5)
# 
# mcf7.seurat <- RunUMAP(mcf7.seurat, dims = 1:10)
# 
# DimPlot(mcf7.seurat, reduction = "umap")
# 
# 
# cluster1.markers <- FindMarkers(mcf7.seurat, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)

