# multimodality analysis
library(Seurat)
library(cowplot)
library(KEGGREST)
library(gage)

source('global_variables.R')

data = readRDS(DATA_NORMED_AND_INTEGRATED)
# pi3k_genes = c('PIK3CA', 'AKT1', 'IRS1',
#                'PDK1', 'MTOR', 'TSC2', 
#                'AKT1S1', 'RPS6KB1', 'EIF4EBP1')
# mapk_genes = c(
#   'RAF1', 'MAP2K1', 'MAPK1', 
#   'MAPK14', 'ESR1', 'ESR2'
# )
# 
# pi3k_time = VlnPlot(data, features = pi3k_genes, group.by = 'time')
# mapk_time = VlnPlot(data, features = mapk_genes, group.by = 'time')
# pi3k_cell_line = VlnPlot(data, features = pi3k_genes, group.by = 'cell_line')
# mapk_cell_line = VlnPlot(data, features = mapk_genes, group.by = 'cell_line')
# 
# pi3k_time_file = file.path(VIOLIN_PLOTS_DIR, 'pi3k_genes_time.jpeg')
# mapk_time_file = file.path(VIOLIN_PLOTS_DIR, 'mapk_genes_time.jpeg')
# pi3k_cell_line_file = file.path(VIOLIN_PLOTS_DIR, 'pi3k_genes_cell_line.jpeg')
# mapk_cell_line_file = file.path(VIOLIN_PLOTS_DIR, 'mapk_genes_cell_line.jpeg')
# 
# ggsave(pi3k_time_file, plot = pi3k_time, height = 15, width = 15)
# ggsave(mapk_time_file, plot = mapk_time, height = 15, width = 15)
# ggsave(pi3k_cell_line_file, plot = pi3k_cell_line, height = 15, width = 15)
# ggsave(mapk_cell_line_file, plot = mapk_cell_line, height = 15, width = 15)
# , split.by = 'cell_line')  

plot_all_vln = function(data, features, ncols = 1, nperplot = 3, prefix='', ...){
  features = features[order(features)]
  start = 1
  remainder = length(features) %% nperplot
  nfull_plots = floor(length(features)/nperplot)
  for (i in 1:length(features)){
    fname = file.path(VIOLIN_PLOTS_DIR, paste0(prefix, features[i], '.jpg'))
    tryCatch(
      {
        v = VlnPlot(data, features = features[i], group.by = 'cell_line',
            ncols=ncols)
      },
      error=function(cond) {
          return(NA)
        }
      )
    ?tryCatch
    print(v)
    ggsave(fname, ...)
    message(paste('saved figure to ', fname))
  }
}

# plot_all_vln(data, features)
# VlnPlot(data, features=features[1:3], group.by = 'time', split.by = 'cell_line')


get_kegg_pathway_genes = function(kegg_pathway){
  kegg = keggGet(kegg_pathway)[[1]]
  kegg$GEN
  genes = vector()
  for (i in kegg$GENE){
    if (is.na(as.numeric(i)))
      genes = c(genes, strsplit(i, ';')[[1]][1])
  }
  return(genes)
}
mapk_genes = get_kegg_pathway_genes('hsa04010')
mapk_genes

plot_all_vln(data, mapk_genes, prefix = 'mapk_kegg')
