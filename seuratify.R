library(Seurat)
library(cowplot)

# split a dataframe into smaller dataframes that have the same colname in the original df
# args
# ----
# df: dataframe to be split. Must have duplicate colnames. In context, the duplicate colnames are time
# returns: a list of dataframes
subset_dataframe = function(df){
  unique_colnames = unique(colnames(df))
  l = list()
  for (i in unique_colnames){
    l[[i]] = df[, colnames(df) == i]
    colnames(l[[i]]) = 1:dim(l[[i]])[2]
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

# iterate over a list of suerat objects and calculate the percentages of mitochondrial contamination per cell
# 
# args
# ----
# seurat_obj_list: list of SeuratObjects. Output from create_seurat_objects
# returns(list): list of SeuratObjects containing the mt.percentages 
calc_mt_percentage = function(seurat_obj_list){
  for (i in names(seurat_obj_list)){
    seurat_obj_list[[i]][['percent.mt']] = PercentageFeatureSet(seurat_obj_list[[i]], pattern = '^MT-')
  }
  return (seurat_obj_list)
}


########################################################################
source('global_variables.R')

# read data back into R. This data was produced in parse_quants.R
data = readRDS(RDS_FILE)
pd = readRDS(PHENO_DATA)

# get pheno data
pd.mcf7 = pd$pd_mcf7.single_cell
pd.t47d = pd$pd_t47d.single_cell

# split data back into cell line data
mcf7 = data$mcf7.single_cell
t47d = data$t47d.single_cell

# extracts counts
mcf7.counts = mcf7$counts
t47d.counts = t47d$counts

# collect genes names that have a signal
mcf7_non0_genes = rownames(mcf7.counts[apply(mcf7.counts!=0, 1, all),])
t47d_non0_genes = rownames(t47d.counts[apply(t47d.counts!=0, 1, all),])
t47d_non0_genes

# rename count matrices with time column
colnames(mcf7.counts) = as.numeric(pd.mcf7$time)
colnames(t47d.counts) = as.numeric(pd.t47d$time)

# split data by time point of estrogen stimulation
mcf7.counts.split = subset_dataframe(mcf7.counts)
t47d.counts.split = subset_dataframe(t47d.counts)

# create some seurat objects. One per time point of stimulation
mcf7.seurat = create_seurat_objects(mcf7.counts.split, 'mcf7')
t47d.seurat = create_seurat_objects(t47d.counts.split, 't47d')
mcf7.seurat
# calculate the percentage of mitochondrial contamination in each time point
mcf7.seurat = calc_mt_percentage(mcf7.seurat)
t47d.seurat = calc_mt_percentage(t47d.seurat)

plot_qc_vln_plots = function(seurat_obj_list, name){
  names = names(seurat_obj_list)

  for (i in 1:length(names)){
    colour_names = names(COLOURS)
    fname = file.path(QC_PLOTS, paste0(name, '_', names[i], '.jpeg'))
    jpeg(filename = fname, height=800, width=800)
    v = VlnPlot(
        seurat_obj_list[[i]],
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1'),
        ncol = 4,
        cols = COLOURS[[colour_names[i]]],
        pt.size = 4
      )  
    # outputs to device
    print(v)
    message(paste('saved to ', fname))
    dev.off()
  }
}
# plot qc metrics
plot_qc_vln_plots(mcf7.seurat, 'mcf7')
plot_qc_vln_plots(t47d.seurat, 't47d')


remove_bad_cells = function(seurat_obj_list, cell='mcf7'){
  if (cell == 'mcf7'){
    for (i in 1:length(seurat_obj_list)){
      seurat_obj_list[[i]] = subset(
        seurat_obj_list[[i]],
        percent.mt < 10 & nCount_RNA > 1e6
        )
    }
  } else if (cell == 't47d'){
    for (i in 1:length(seurat_obj_list)){
      seurat_obj_list[[i]] = subset(
        seurat_obj_list[[i]],
        percent.mt < 20
        )
      }
    }else {
      stop('me no likey')
    }
  return (seurat_obj_list)
}

mcf7.seurat.subset = remove_bad_cells(mcf7.seurat, cell='mcf7')
t47d.seurat.subset = remove_bad_cells(t47d.seurat, cell='t47d')

# replot qc metrics
plot_qc_vln_plots(mcf7.seurat.subset, 'mcf7')
plot_qc_vln_plots(t47d.seurat.subset, 't47d')

# print out object metadata
x = mcf7.seurat.subset$`0`@meta.data
x = x[order(x['nCount_RNA'], decreasing = T), ]
x

y = t47d.seurat.subset$`12`@meta.data
y = y[order(y['nFeature_RNA'], decreasing = T), ]
y

normalise_data = function(seurat_obj_list){
  names = names(seurat_obj_list)
  for (i in 1:length(names)){
    seurat_obj_list[[names[i]]] = NormalizeData(seurat_obj_list[[names[i]]]) 
  }
  return(seurat_obj_list)
}
find_variable_features = function(seurat_obj_list){
  names = names(seurat_obj_list)
  for (i in 1:length(names)){
    seurat_obj_list[[names[i]]] = FindVariableFeatures(
      seurat_obj_list[[names[i]]], selection.method = "vst", nfeatures = 2000
      ) 
  }
  return(seurat_obj_list)
}

mcf7.seurat.norm = normalise_data(mcf7.seurat.subset)
t47d.seurat.norm = normalise_data(t47d.seurat.subset)
t47d.seurat.norm

mcf7.anchors = FindIntegrationAnchors(mcf7.seurat, dims=1:15)
t47d.anchors = FindIntegrationAnchors(t47d.seurat, dims=1:15)
mcf7.seurat.subset

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

