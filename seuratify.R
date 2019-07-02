library(Seurat)
library(cowplot)

# Take a dataframe as input. Expects the dataframe to have multiple columns labelled the same
# Output list of dataframes that has split the big dataframe by unique colnames
subset_dataframe = function(df){
  unique_colnames = unique(colnames(df))
  l = list()
  for (i in unique_colnames){
    l[[i]] = df[, colnames(df) == i]
  }
  return (l)
}

# Iterate over a list of count matrices. Rename each cell numerically and create a SeuratObject
# 
# args
# ----
# - counts_list: list of matrices to be grouped in a seurat object
# - project_name: string. prepended to labels of output list
# - returns a list of seurat objects
create_seurat_objects = function(counts_list, project_name){
  l = list()
  print(counts_list)
  for (i in names(counts_list)){
    # print(counts_list[[i]])
    c = counts_list[[i]]
    colnames(c) = 1:dim(c)[2]
    l[[i]] = CreateSeuratObject(
      counts=c,
      project = paste(project_name, as.character(i), sep='_'),
      min.cells = 3, min.features = 200
    )
  }
  return (l)
}


########################################################################
source('global_variables.R')

# read data back into R. This data was produced in parse_quants.R
data = readRDS(RDS_FILE)
pd = readRDS(PHENO_DATA)

# get pheno data
pd.mcf7 = pd$pd_mcf7.single_cell
pd.t47d = pd$pd_t47d.single_cell

# 
mcf7_non0_genes = rownames(mcf7$counts[apply(mcf7$counts!=0, 1, all),])
t47d_non0_genes = rownames(t47d$counts[apply(t47d$counts!=0, 1, all),])
t47d_non0_genes

# split data back into cell line data
mcf7 = data$mcf7.single_cell
t47d = data$t47d.single_cell

# extracts counts
mcf7.counts = mcf7$counts
t47d.counts = t47d$counts

# rename count matrices with time column
colnames(mcf7.counts) = as.numeric(pd.mcf7$time)
colnames(t47d.counts) = as.numeric(pd.t47d$time)


mcf7.counts.split = subset_dataframe(mcf7.counts)
t47d.counts.split = subset_dataframe(t47d.counts)



mcf7.seurat = create_seurat_objects(mcf7.counts.split, 'mcf7')
t47d.seurat = create_seurat_objects(t47d.counts.split, 't47d')

# iterate over a list of suerat objects and calculate the percentages of mitochondrial contamination per cell
# 
# args
# ----
# seurat_obj_list: list of SeuratObjects. Output from create_seurat_objects
# returns(list): list of SeuratObjects containing the mt.percentages 
calc_mt_percentage = function(seurat_obj_list){
  l = list()
  for (i in names(seurat_obj_list)){
    seurat_obj_list[[i]][['percent.mt']] = PercentageFeatureSet(seurat_obj_list[[i]], pattern = '^MT-')
  }
  return (l)
}

calc_mt_percentage(mcf7.seurat)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mcf7.seurat[["percent.mt"]] <- PercentageFeatureSet(mcf7.seurat, pattern = "^MT-")
t47d.seurat[["percent.mt"]] <- PercentageFeatureSet(t47d.seurat, pattern = "^MT-")

FindIntegrationAnchors(mcf7.seurat, dims=1:15)

'
machine learning problem
------------------------
- Estimate dublets, dead cells or otherwise poor quality from number of features
   number of gene counts and percentage mt genes
- How to generalize to other datasets?
- Unsupervised would be best

QC metrics for single cell 
--------------------------
- The number of unique genes detected in each cell (nFeature_RNA). 
  - When this metric is very low, a droplet may have been empty indicating low quality cell
  - When very high, may be dublet cell
- Total number of molecules detected (nCount) in a cell indicates similarly to nFeatures_RNA
- mitochondrial genome
  - low-quality or dying cells may have extensive mt geneome contamination
  - 
'
x = mcf7.seurat@meta.data[order(mcf7.seurat@meta.data$percent.mt, decreasing = T),]
y = t47d.seurat@meta.data[order(t47d.seurat@meta.data$percent.mt, decreasing = T),]
x$nCount_RNA = x$nCount_RNA / 10000
x['nFeat*nCount'] = order(x$nFeature_RNA * x$nCount_RNA)
x = x[order(x$`nFeat*nCount`), ]
x
y
?order
x
x
hist(x$nFeature_RNA)
?hist
plot(x$nFeature_RNA, x$nCount_RNA)
# Visualize QC metrics as a violin plot
mcf7.seurat = subset(mcf7.seurat, percent.mt < 10)

t47d.seurat = subset(t47d.seurat, percent.mt < 25 
                     & nFeature_RNA > 10000,
)

VlnPlot(mcf7.seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1'), 
        ncol = 4
)

VlnPlot(t47d.seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1'), 
        ncol = 4
)


y = y[order(y$percent.mt, decreasing = T),]
y
?VlnPlot

?VlnPlot



plot1 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


plot1 <- FeatureScatter(t47d.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(t47d.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# mcf7.seurat <- subset(mcf7.seurat, subset = nFeature_RNA > 3000)


mcf7.seurat <- NormalizeData(mcf7.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
t47d.seurat <- NormalizeData(t47d.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
mcf7.seurat
t47d.seurat

VlnPlot(mcf7.seurat, 
        features = c(sample(mcf7_non0_genes, 5), 'IDO1'),
        ncol = 3
)

mcf7.seurat

VlnPlot(t47d.seurat, 
        features = c(sample(t47d_non0_genes, 5), 'IDO1'),
        ncol = 3
)

mcf7.seurat <- FindVariableFeatures(mcf7.seurat, selection.method = "vst", nfeatures = 2000)
mcf7.seurat@reductions

top10 = head(VariableFeatures(mcf7.seurat))

top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mcf7.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(mcf7.seurat)
mcf7.seurat <- ScaleData(mcf7.seurat, features = all.genes)
mcf7.seurat <- RunPCA(mcf7.seurat, features = VariableFeatures(object = mcf7.seurat))

print(mcf7.seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(mcf7.seurat, dims = 1:2, reduction = "pca")


DimPlot(mcf7.seurat, reduction = "pca")

DimHeatmap(mcf7.seurat, dims = 1:2, cells = 500, balanced = TRUE)

mcf7.seurat <- JackStraw(mcf7.seurat, num.replicate = 100)
mcf7.seurat <- ScoreJackStraw(mcf7.seurat, dims = 1:20)
JackStrawPlot(mcf7.seurat, dims = 1:20)
ElbowPlot(mcf7.seurat)


mcf7.seurat <- FindNeighbors(mcf7.seurat, dims = 1:10)
mcf7.seurat <- FindClusters(mcf7.seurat, resolution = 0.5)
head(Idents(mcf7.seurat), 5)

mcf7.seurat <- RunUMAP(mcf7.seurat, dims = 1:10)

DimPlot(mcf7.seurat, reduction = "umap")


cluster1.markers <- FindMarkers(mcf7.seurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

