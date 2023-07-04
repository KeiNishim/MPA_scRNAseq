#For plotting heatmap with ISG/monocyte signature###
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#For analysis of monocytes
#Road "mono_ind.rds" as s.sub
DefaultAssay(s.sub) <- "RNA"
s.sub$seurat_clusters <- factor(s.sub$seurat_clusters,
                                levels= c("cDC", "CD16 Mono", "CD14 Mono_HLA", "CD14 Mono_ISG","CD14 Mono_VCAN","CD14 Mono_Activated"
                                ))
pp<-DotPlot(s.sub,cols=c("blue","red"),col.max=2.5, col.min=-2.5, scale=T,dot.scale=3,
            features=c("PLBD1", "ALOX5AP", "PADI4", "TSPO", "HP",
                       "S100A8", "S100A12", "S100A9",  "CTSD", "VCAN",
                       "IFI6", "IFI44L", "MX1", "LY6E", "ISG15", "XAF1",
                       "LGALS2", "HLA-DRB1", "HLA-DMA", "CD74", "HLA-DPB1", "FCER1A",
                       "CDKN1C", "LYPD2", "RHOC", "MS4A7", "LST1",
              "CLEC10A", "CD1C", "CST3", "ENHO", "PLD4"),
group.by="seurat_clusters")+theme(axis.text.x = element_text(size = 9, angle=90)) +
    scale_colour_gradient2(low="blue", mid="white", high="red")
pp
ggsave('Fig.2A.png',pp, width=8.5, height=2.3)

DefaultAssay(s.sub)<-"ADT"
s.sub <- ScaleData(s.sub)
s.sub$seurat_clusters <- factor(s.sub$seurat_clusters,
                                levels= c("CD14 Mono_Activated","CD14 Mono_VCAN","CD14 Mono_ISG","CD14 Mono_HLA","CD16 Mono","cDC"
                                ))
p <-VlnPlot(subset(s.sub, downsample=200),group.by="seurat_clusters",ncol=4,
            features=c("ADT-CD14.1","ADT-HLA-DR","ADT-FceRIa", "ADT-CD16"),
            slot="data")
p
ggsave('Fig.S3B.png', p, width=12, height=6.5)
