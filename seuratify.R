# create seurat object, do normalisation, scaling and integrate mcf7 and t47d cells into 
#  one dataset
library(Seurat)

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
  seurat_obj = FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  return(seurat_obj)
}
?NormalizeData
?ScaleData
mcf7.seurat.norm = normalise_data(mcf7.seurat.subset)
t47d.seurat.norm = normalise_data(t47d.seurat.subset)

# mcf7.seurat.not_norm = mcf7.seurat.subset
# t47d.seurat.not_norm = t47d.seurat.subset


# t47d.seurat.norm <- RunPCA(t47d.seurat.norm, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
t47d.seurat.norm <- RunPCA(t47d.seurat.norm, features = VariableFeatures(t47d.seurat.norm))
t47d.seurat.cc <- CellCycleScoring(t47d.seurat.norm, s.features = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes, set.ident = T)
p1 = DimPlot(t47d.seurat.cc)
t47d.seurat.cc.scaled <- ScaleData(t47d.seurat.cc, vars.to.regress = c("S.Score", "G2M.Score"), 
                            features = rownames(t47d.seurat.cc))
t47d.seurat.cc.scaled <- RunPCA(t47d.seurat.cc.scaled, features = c(cc.genes$s.genes, cc.genes$g2m.genes))

p2 = DimPlot(t47d.seurat.cc.scaled)
CombinePlots(list(p1, p2))
# t47d.seurat.cc
# head(t47d.seurat.cc[[]])

correct_for_cell_cycle= function(seurat_obj, cell, ...){
  seurat_obj = RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  suerat_obj.cc <- CellCycleScoring(seurat_obj, s.features = cc.genes$s.genes,
                                 g2m.features = cc.genes$g2m.genes, set.ident = T)
  p1 = DimPlot(suerat_obj.cc, pt.size = 5)
  suerat_obj.cc.scaled <- ScaleData(suerat_obj.cc, vars.to.regress = c("S.Score", "G2M.Score"), 
                          features = rownames(suerat_obj.cc))
  suerat_obj.cc.scaled = RunPCA(suerat_obj.cc.scaled, features = VariableFeatures(seurat_obj))
  p2 = DimPlot(suerat_obj.cc.scaled, pt.size = 5)
  p = CombinePlots(list(p1, p2))
  print(p)
  fname1 = file.path(QC_PLOTS, paste(cell, 'cell_cycle_pca.jpeg'))
  ggsave(fname1, ...)
  return(suerat_obj.cc.scaled)
}
?VariableFeatures
?DimPlot
mcf7.cc_corrected = correct_for_cell_cycle(mcf7.seurat.norm, 'mcf7', width=15, height = 15)
t47d.cc_corrected = correct_for_cell_cycle(t47d.seurat.norm, 't47d', width=15, height = 15)

# RidgePlot(t47d.cc_corrected, features = cc.genes$g2m.genes[1:6], ncol = 3)
# 
# RidgePlot(t47d.cc_corrected, features = cc.genes$g2m.genes[1:6], ncol = 3, 
#           group.by = 'orig.ident')
# RidgePlot(mcf7.cc_corrected, features = cc.genes$g2m.genes[1:6], ncol = 3, 
#           group.by = 'orig.ident')
# RidgePlot(t47d.seurat, features = cc.genes$g2m.genes[1:6], ncol = 3)
# 
# plot1 <- FeatureScatter(mcf7.seurat.norm, feature1 = "BRCA1", feature2 = 'SMAD3')
# plot2 <- FeatureScatter(mcf7.cc_corrected, feature1 = "BRCA1", feature2 = 'SMAD3')
# plot1 <- FeatureScatter(t47d.seurat.norm, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(t47d.cc_corrected, feature1 = "nCount_RNA", feature2 = "percent.mt")
# CombinePlots(list(plot1, plot2))

# 
find_integration_anchors = function(seurat_list){
  anchors = FindIntegrationAnchors(
    seurat_list, 
    anchor.features = 4000, 
    k.filter = 50
    )
  return(anchors)
}

anchors.norm = find_integration_anchors(list(mcf7.cc_corrected, t47d.cc_corrected))
data.norm = IntegrateData(anchors.norm)
?IntegrateData
# save normalised, scaled, cell cycle corrected and integrated data to rds
saveRDS(data.norm, DATA_NORMED_AND_INTEGRATED)

