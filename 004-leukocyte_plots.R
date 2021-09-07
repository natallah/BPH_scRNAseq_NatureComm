input_file <- "data/combined_seurat.RDS" 

library(Seurat)
library(tidyverse)
library(ggplot2)

library(ggforce)
library(ComplexHeatmap)
library(ggrepel)
library(patchwork)
library(viridis)
library(ggdendro)
library(cowplot)
library(ggtree)
library(tidyr)

theme_set(theme_cowplot())
seurat_obj <- readRDS(input_file)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)
seurat_obj@meta.data[["Subtype"]] <- gsub(".*_(.*)", "\\1", seurat_obj@meta.data$Sample)

seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
seurat_obj <- ScaleData(seurat_obj, assay = "ADT")

#dir.create("paperFigs_renee_Dec20")
seurat_obj@meta.data[["Subtype"]] <- gsub(".*_(.*)", "\\1", seurat_obj@meta.data$Sample)
seurat.temp <- seurat_obj
Idents(seurat.temp)<- seurat.temp@meta.data$Subtype
small_cells <- WhichCells(seurat.temp, idents='small')
large_cells <- WhichCells(seurat.temp, idents='large')

seurat.small <- subset(seurat.temp,cells=small_cells)
seurat.large <- subset(seurat.temp,cells=large_cells)

identities.all <- levels(seurat_obj)
my_color_palette.all <- hue_pal()(length(identities.all))

png(filename="sample_small_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat.small, reduction = "umap", group.by = "Sample")+scale_color_manual(values = my_color_palette.all[c(1,2,3,4,5,6,7,9,10,13)])
dev.off()

png(filename="sample_large_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat.large, reduction = "umap", group.by = "Sample")+scale_color_manual(values = my_color_palette.all[c(8,11,12,14)])
dev.off()

png(filename="sample_large_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap", group.by = "Sample")+scale_color_manual(values = my_color_palette.all)
dev.off()

png(filename="type_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap", group.by = "Subtype")
dev.off()

png(filename="umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap")
dev.off()

png(filename="umap_IL7R.png", width=10, height=8,units='in',res=600)
FeaturePlot(seurat_obj, features = "IL7R")
dev.off()

png(filename="CITE_seq_featurePlot_use.png", width=11, height=12,units='in',res=600)
FeaturePlot(seurat_obj, features = c("ADT_CD3", "ADT_CD19", "ADT_CD4", "ADT_CD8a","ADT_CD11b"),min.cutoff = "q05", max.cutoff="q95", ncol = 2)
dev.off()

all.markers <- FindAllMarkers(seurat_obj)

png(filename="DotPlot_seurat.png", width=14, height=20,units='in',res=600)
DotPlot(seurat_obj,features=markers)+scale_color_gradientn(colors=viridis::viridis(20,direction = -1), limits=c(0,4), oob=scales::squish, name='log2 (count + 1)')
dev.off()

gene_cluster$id <- as.factor(gene_cluster$id)
png(filename="DotPlot_noFilter_switchedXY_small.png", width=13, height=6,units='in',res=600)
gene_cluster %>% filter(Gene %in% markers) %>% 
    mutate(`% Expressing`=cell_exp_ct/cell_ct*100) %>%
    mutate(Gene = factor(Gene, levels=markers)) %>%
    #filter(Count > 0, `% Expressing` >1 ) %>%
    ggplot(aes(x=Gene,y=id, color=Count, size = `% Expressing`))+
    geom_point() + scale_color_viridis_c(name= 'log2 (count + 1)')+
    cowplot::theme_cowplot() + theme(axis.line = element_blank())+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
    ylab('')+theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colors=viridis::viridis(20,direction = -1), limits=c(0,4), oob=scales::squish, name='log2 (count + 1)')
dev.off()
library(vctrs)
mat <- gene_cluster %>% filter(Gene %in% markers) %>%
    select(-cell_ct,-cell_exp_ct,-group) %>%
    pivot_wider(names_from= id, values_from = Count) %>%
    data.frame()
row.names(mat) <- mat$Gene
mat <- mat[,-1]
clust <- hclust(dist(mat %>% as.matrix()))

ddgram <- as.dendrogram(clust)
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

macGenes <- c("THBS1","CD36","CD47","IL1B","TGFB1","IL6","CXCL8","EREG","IL4","IL17A","IL18","IFNA1","IFNG",
              "HIF1A","VEGFA","TNF","TNFRSF1A","TNFRSF1B")

g <- FeaturePlot(
    seurat_obj, features = "CD68", 
    order = TRUE,
    cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "genePlots_CD68.png", width = 10, height = 8, units='in' ,dpi=600)


g <- FeaturePlot(
  seurat_obj, features = "TNFRSF1B", 
  order = TRUE,
  cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "genePlots_TNFRSF1B.png", width = 10, height = 8, units='in' ,dpi=600)

g <- VlnPlot(
    seurat_obj, features = "CD4")
ggsave(plot = g, "vlnPlot_CD4.png", width = 10, height = 8, units='in' ,dpi=600)


g <- FeaturePlot(
    seurat_obj, features = "CD4", 
    order = TRUE,
    cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "genePlots_CD4.png", width = 10, height = 8, units='in' ,dpi=600)


png(filename="genePlots_TNFRSF1B_splitVln_nopoints.png",height=3,width=4,res=600,units='in')
VlnPlot(seurat_obj, features="TNFRSF1B",split.by="Subtype",pt.size = 0.000)
dev.off()

png(filename="genePlots_TNFRSF1B_Vln_nopoints.png",height=4,width=4,res=600,units='in')
VlnPlot(seurat_obj, features="TNFRSF1B",pt.size = 0.000)
dev.off()
