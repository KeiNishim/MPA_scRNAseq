#Clustering for monocyte sub-popultion###
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#Road "s.int.rds" as s.int
Idents(s.int) <- "sample"
#####For dimplot using individual samples
s.int <- subset(s.int, ident=c("AI017-1","AI028-1", "AI035-1", "AI039-1","AI054-1", "AI055-1", "AI062-1", "AI075-1", "HC037-1", "HC038-1", "HC039-1", "HC040-1",
                               "HC041-1", "HD020-3", "HD036-1")) 
s.int <- RenameIdents(s.int, "AI017-1"="MPA-2","AI028-1"="MPA-3", "AI035-1"="MPA-8", "AI039-1"="MPA-7",
                      "AI054-1"="MPA-1", "AI055-1"="MPA-6", "AI062-1"="MPA-4", "AI075-1"="MPA-5",
                      "HD020-3"="HD-1", "HD036-1"="HD-2",
                      "HC037-1"="HD-3", "HC038-1"="HD-4", "HC039-1"="HD-5", "HC040-1"="HD-6",
                      "HC041-1"="HD-7","AI054-2"="MPA-5-post","AI075-2"="MPA-8-post")
s.int[["sample"]] <- Idents(s.int)
s.int$sample <- factor(s.int$sample,
                       levels= c("MPA-1","MPA-2","MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8", "HD-1", "HD-2", "HD-3", "HD-4", "HD-5", "HD-6", "HD-7"))

## add clustering step
Idents(s.int) <- s.int[["predicted.celltype.l2"]]
s.int <- subset(s.int, ident=c("CD8 Naive",hCD8 TCMh,"CD8 TEM"))
s.int <- RunUMAP(s.int, reduction = 'ref.spca', dims = 1:30)

s.int <- FindNeighbors(s.int,reduction = "ref.spca",dims = 1:30,k.param = 10,#return.neighbor=TRUE,
                       compute.SNN = TRUE,prune.SNN = 1/20)
s.int <- FindClusters(s.int, resolution=1, reduction="umap")
DimPlot(s.int, repel=T, label=T)

#Manually remove doublet and debris using gene expression data and CITE-seq data. Then re-cluster and annotate population (Omitted)

s.int[["seurat_clusters"]] <- Idents(s.int)
s.int$seurat_clusters <- factor(s.int$seurat_clusters,
                                levels= c("CD8 T_Naive", "CD8 T_CM", "CD8 T_EM", "CD8 T_CTL", "CD8 T_KIR"
                                ))
saveRDS(s.int, "CD8_ind.rds")
xp <- DimPlot(s.int, group.by="seurat_clusters",reduction="umap",shuffle=T, pt.size=0.1)
xp
ggsave('Fig.2F.png',xp, width=7.5, height=7.5)

#Individual mapping
Idents(s.int) <- "sample"
s.sub <- subset(s.int, downsample=1000, idents=c("MPA-1","MPA-2","MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8", 
                                                 "HD-1", "HD-2", "HD-3", "HD-4", "HD-5", "HD-6", "HD-7"))
pp <- DimPlot(s.sub,split.by="sample", ncol=8,group.by="seurat_clusters",reduction="umap",shuffle=T, pt.size=0.2)
pp
ggsave('Fig. S4A.png',pp,width=20,height=7.5)
write.csv(table(s.int@meta.data[["sample"]], s.int@meta.data[["seurat_clusters"]]), "CD8T_clusters.csv")
#################################################

