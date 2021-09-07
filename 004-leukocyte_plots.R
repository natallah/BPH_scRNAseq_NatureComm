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
seurat_obj <- readRDS("../step1_cellRangerCounts/cellranger_seurat/combinedOct14/combined_seurat.RDS")
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

png(filename="paperFigs_renee_Dec20/sample_large_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap", group.by = "Sample")+scale_color_manual(values = my_color_palette.all)
dev.off()

png(filename="paperFigs_renee_Dec20/type_umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap", group.by = "Subtype")
dev.off()

png(filename="paperFigs_renee_Dec20/umap.png", width=10, height=8,units='in',res=600)
DimPlot(seurat_obj, reduction = "umap")
dev.off()

png(filename="paperFigs_renee_Dec20/umap_IL7R.png", width=10, height=8,units='in',res=600)
FeaturePlot(seurat_obj, features = "IL7R")
dev.off()

png(filename="paperFigs_renee_Dec20/CITE_seq_featurePlot_use.png", width=11, height=12,units='in',res=600)
FeaturePlot(seurat_obj, features = c("ADT_CD3", "ADT_CD19", "ADT_CD4", "ADT_CD8a","ADT_CD11b"),min.cutoff = "q05", max.cutoff="q95", ncol = 2)
dev.off()

#all.markers <- FindAllMarkers(seurat_obj)
all.markers <- read.csv("all_cluster_markers_Oct14.csv",row.names=1)
cluster.0.3 <- c(0:3)
cluster.4 <- c(4)
cluster.5.7 <- c(5:7)
cluster.8 <- c(8)
cluster.9.10 <- c(9,10)
cluster.11 <- c(11)
cluster.12.13 <- c(12:13)
all.markers2 <- all.markers[grepl("^[^RP]",all.markers$gene),]
#write.csv(all.markers, "all_cluster_markers_Oct14.csv")
all.markers.0.3 <- all.markers2[all.markers2$cluster %in% cluster.0.3,]
all.markers.4 <- all.markers2[all.markers2$cluster %in% cluster.4,]
all.markers.5.7 <- all.markers2[all.markers2$cluster %in% cluster.5.7,]
all.markers.8 <- all.markers2[all.markers2$cluster %in% cluster.8,]
all.markers.9.10 <- all.markers2[all.markers2$cluster %in% cluster.9.10,]
all.markers.11 <- all.markers2[all.markers2$cluster %in% cluster.11,]
all.markers.12.13 <- all.markers2[all.markers2$cluster %in% cluster.12.13,]


top0.3 <- all.markers.0.3 %>% group_by(cluster) %>% top_n(n=4, wt=avg_logFC)
top4 <- all.markers.4 %>% group_by(cluster) %>% top_n(n=6, wt=avg_logFC)
top5.7 <- all.markers.5.7 %>% group_by(cluster) %>% top_n(n=4, wt=avg_logFC)
top8 <- all.markers.8 %>% group_by(cluster) %>% top_n(n=6, wt=avg_logFC)
top9.10 <- all.markers.9.10 %>% group_by(cluster) %>% top_n(n=4, wt=avg_logFC)
top11 <- all.markers.11 %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
top12.13 <- all.markers.12.13 %>% group_by(cluster) %>% top_n(n=2, wt=avg_logFC)
top.all <- rbind(top0.3,top4,top5.7,top8,top9.10,top11,top12.13)
top.all[duplicated(top.all$gene),]
dot <- DotPlot(seurat_obj, features=unique(top.all$gene))
gene_cluster <- dot$data
gene_cluster$Gene <- gene_cluster$features.plot


tt <- table(seurat_obj@meta.data$seurat_clusters)
gene_cluster$cell_ct <- c(rep(tt[[1]],50),rep(tt[[2]],50),rep(tt[[3]],50),rep(tt[[4]],50),rep(tt[[5]],50),rep(tt[[6]],50),
                      rep(tt[[7]],50),rep(tt[[8]],50),rep(tt[[9]],50),rep(tt[[10]],50),rep(tt[[11]],50),rep(tt[[12]],50),
                      rep(tt[[13]],50),rep(tt[[14]],50))
gene_cluster$group <- gene_cluster$id
gene_cluster$cell_exp_ct <- (gene_cluster$pct.exp/100)*gene_cluster$cell_ct
temp <- log2(gene_cluster$avg.exp)
temp[is.na(temp)] <- 0
gene_cluster$Count <- temp
#write.csv(gene_cluster,"Gene_Cluster_Oct14_small.csv")
gene_cluster <- read.csv("Gene_Cluster_Oct14.csv")

markers <- gene_cluster$Gene %>% unique()
gene_cluster2 <- gene_cluster

png(filename="paperFigs_renee_Dec20/DotPlot_seurat.png", width=14, height=20,units='in',res=600)
DotPlot(seurat_obj,features=markers)+scale_color_gradientn(colors=viridis::viridis(20,direction = -1), limits=c(0,4), oob=scales::squish, name='log2 (count + 1)')
dev.off()
gene_cluster$id <- as.factor(gene_cluster$id)
png(filename="paperFigs_renee_Dec20/DotPlot_noFilter_switchedXY_small.png", width=13, height=6,units='in',res=600)
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

dir.create("mac_plots_renee_Oct14")
macGenes <- c("THBS1","CD36","CD47","IL1B","TGFB1","IL6","CXCL8","EREG","IL4","IL17A","IL18","IFNA1","IFNG",
              "HIF1A","VEGFA","TNF","TNFRSF1A","TNFRSF1B")

g <- FeaturePlot(
    seurat_obj, features = "CD68", 
    order = TRUE,
    cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "paperFigs_renee_Dec20/genePlots_CD68.png", width = 10, height = 8, units='in' ,dpi=600)


g <- FeaturePlot(
  seurat_obj, features = "TNFRSF1B", 
  order = TRUE,
  cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "paperFigs_renee_Dec20/genePlots_TNFRSF1B.png", width = 10, height = 8, units='in' ,dpi=600)

g <- VlnPlot(
    seurat_obj, features = "CD4")
ggsave(plot = g, "renee_plots_april2021/vlnPlot_CD4.png", width = 10, height = 8, units='in' ,dpi=600)


g <- FeaturePlot(
    seurat_obj, features = "CD4", 
    order = TRUE,
    cols = viridis::viridis(7, direction = -1))
ggsave(plot = g, "renee_plots_april2021/genePlots_CD4.png", width = 10, height = 8, units='in' ,dpi=600)


png(filename="paperFigs_renee_Dec20/genePlots_TNFRSF1B_splitVln_nopoints.png",height=3,width=4,res=600,units='in')
VlnPlot(seurat_obj, features="TNFRSF1B",split.by="Subtype",pt.size = 0.000)
dev.off()

png(filename="paperFigs_renee_Dec20/genePlots_TNFRSF1B_Vln_nopoints.png",height=4,width=4,res=600,units='in')
VlnPlot(seurat_obj, features="TNFRSF1B",pt.size = 0.000)
dev.off()
