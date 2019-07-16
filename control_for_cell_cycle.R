library(Seurat)

source('global_variables.R')


data = readRDS(DATA_NORMED_AND_INTEGRATED)
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

data <- RunPCA(data, features = VariableFeatures(data))
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(data[[]])

variable_features = VariableFeatures(data)
variable_features
features = sample(variable_features, 6)
features
RidgePlot(data, features = features, ncol=3)

VlnPlot(data, features = c('GINS2'))
data <- RunPCA(data, features = c(s.genes, g2m.genes))
DimPlot(data, reduction = 'pca')


















