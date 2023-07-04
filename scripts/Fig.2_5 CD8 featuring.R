#For plotting heatmap with ISG/monocyte signature###
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)


#For analysis of CD8 T cells
#Read CD8_ind.rds as s.sub
DefaultAssay(s.sub) <- "RNA"
s.sub$seurat_clusters <- factor(s.sub$seurat_clusters,
                                levels= c("CD8 T_KIR", "CD8 T_CTL", "CD8 T_EM", "CD8 T_CM", "CD8 T_Naive"
                                ))
pp<-DotPlot(s.sub,cols=c("blue","red"),col.max=2.5, col.min=-2.5, scale=T,dot.scale=3,
            features=c("LEF1","LTB","CCR7","TCF7","SELL","NOSIP",
                       "IL7R","CRIP2","CAPG","MAL","RCAN3","LDHB",
                       "GZMK","DUSP2","CD74","IFNG-AS1","COTL1","GPR183",
                       "GZMH","FGFBP2","NKG7","GZMB","PRF1","EFHD2","GNLY",
                       "TYROBP","KLRC2","KLRC3","KIR2DL3","CTSW"),
            group.by="seurat_clusters")+theme(axis.text.x = element_text(size = 9, angle=90)) +
  scale_colour_gradient2(low="blue", mid="white", high="red")
pp
ggsave('Fig.2E.png',pp, width=8.5, height=2.3)

DefaultAssay(s.sub)<-"ADT"
s.sub <- ScaleData(s.sub)
s.sub$seurat_clusters <- factor(s.sub$seurat_clusters,
                                levels= c("CD8 T_Naive", "CD8 T_CM", "CD8 T_EM", "CD8 T_CTL", "CD8 T_KIR"
                                ))
p <-VlnPlot(subset(s.sub, downsample=200),group.by="seurat_clusters",ncol=5,
              features=c("ADT-CD197","ADT-CD27.1","ADT-CD45RA","ADT-CD45RO", "ADT-CD56"),
                         slot="data")
p
ggsave('Fig.S4B.png', p, width=16, height=6.5)
