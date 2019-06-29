library(Seurat)


########################################################################
source('global_variables.R')

data = readRDS(RDS_FILE)
mcf7 = data$mcf7$single_cell
t47d = data$t47d$single_cell


mcf7.seurat <- CreateSeuratObject(counts = mcf7$counts, project = "mcf7", min.cells = 3, min.features = 200)
mcf7.seurat

tot_counts = colSums(mcf7$counts)
barplot(tot_counts)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mcf7.seurat[["percent.mt"]] <- PercentageFeatureSet(mcf7.seurat, pattern = "^MT-")

mcf7.seurat@meta.data

# Visualize QC metrics as a violin plot
VlnPlot(mcf7.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mcf7.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


mcf7.seurat <- subset(mcf7.seurat, subset = nFeature_RNA > 3000)












