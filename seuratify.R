# create seurat object, do normalisation, scaling and integrate mcf7 and t47d cells into 
#  one dataset
library(Seurat)
library(cowplot)


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

# Create a seurat object
# args
# ----
# counts: a counts matrix. Output from tximport
# pd:     phenodata. Metadata about each of the cells in counts
create_seurat_object = function(counts, pd, project_id, ...){
  seurat_obj = CreateSeuratObject(
    counts, project = project_id, 
    meta.data = pd, 
    min.cells = 2, 
    names.field = 1,
    ...
    )
  seurat_obj[['percent.mt']] = PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return (seurat_obj)
}
mcf7.seurat = create_seurat_object(mcf7.counts, pd.mcf7, 'mcf7')
t47d.seurat = create_seurat_object(t47d.counts, pd.t47d, 't47d')
pd.t47d
# Plot 
qc_violins = function(seurat_obj){
  v = VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'IDO1', 'SMAD3'),
    ncol = 3,
    pt.size = 4,
    group.by = 'time'
  )
  fname = file.path(QC_PLOTS, paste0(seurat_obj@meta.data$orig.ident[1], '.jpeg'))
  ggsave(fname, device = 'jpeg', width = 12, height = 6.67, dpi = 300)
  message(paste('Plot saved to', fname))
}

qc_violins(mcf7.seurat)
qc_violins(t47d.seurat)


remove_bad_cells = function(seurat_obj, cell='mcf7'){
  if (cell == 'mcf7'){
      seurat_obj = subset(
        seurat_obj,
        percent.mt < 10 & nCount_RNA > 1e6 & nFeature_RNA > 10000
        )
    
  } else if (cell == 't47d'){
      seurat_obj = subset(
        seurat_obj,
        percent.mt < 20 & nCount_RNA > 2e6 & nFeature_RNA > 10000
        )
    }else {
      stop('me no likey')
    }
  return (seurat_obj)
}

mcf7.seurat.subset = remove_bad_cells(mcf7.seurat, cell='mcf7')
t47d.seurat.subset = remove_bad_cells(t47d.seurat, cell='t47d')

# replot qc metrics
qc_violins(mcf7.seurat.subset)
qc_violins(t47d.seurat.subset)

normalise_data = function(seurat_obj){
  seurat_obj = NormalizeData(seurat_obj)
  return(seurat_obj)
}

mcf7.seurat.norm = normalise_data(mcf7.seurat.subset)
t47d.seurat.norm = normalise_data(t47d.seurat.subset)

mcf7.seurat.not_norm = mcf7.seurat.subset
t47d.seurat.not_norm = t47d.seurat.subset

plot1 <- FeatureScatter(mcf7.seurat.norm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(t47d.seurat.norm, feature1 = "nCount_RNA", feature2 = "percent.mt")
CombinePlots(list(plot1, plot2))

find_integration_anchors = function(seurat_list){
  anchors = FindIntegrationAnchors(
    seurat_list, 
    anchor.features = 2000, 
    k.filter = 50
    )
  return(anchors)
}
anchors.norm = find_integration_anchors(list(mcf7.seurat.norm, t47d.seurat.norm))
anchors.not_norm = find_integration_anchors(list(mcf7.seurat.not_norm, t47d.seurat.not_norm))
data.norm = IntegrateData(anchors.norm)
data.not_norm = IntegrateData(anchors.not_norm)

saveRDS(data.norm, DATA_NORMED_AND_INTEGRATED)

