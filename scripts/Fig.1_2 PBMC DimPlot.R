## Read in libraries --------------------------------
library(Seurat)
library(SeuratDisk)
library(cowplot)
library(dplyr)
library(ggplot2)
library(Polychrome)
library(patchwork)
library(BiocManager)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(scran)

Idents(s.int) <- "predicted.celltype.l2"
s.int <- subset(s.int, ident=c("ASDC", "B intermediate", "B memory", "B naive", "CD14 Mono", "CD16 Mono", "CD4 CTL", "CD4 Naive","CD4 Proliferating", "CD4 TCM",
                               "CD4 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "HSPC", "ILC", "MAIT", "NK", "NK Proliferating",
                               "NK_CD56bright", "Plasmablast", "Treg", "cDC1", "cDC2", "dnT", "gdT", "pDC"))
Idents(s.int) <- "sample"
s.int <- subset(s.int, ident=c("AI017-1","AI028-1", "AI035-1", "AI039-1","AI054-1", "AI055-1", "AI062-1", "AI075-1",
                               "HC037-1", "HC038-1", "HC039-1", "HC040-1",
                               "HC041-1", "HD020-3", "HD036-1")) 
s.int <- RenameIdents(s.int, "AI017-1"="MPA-2","AI028-1"="MPA-3", "AI035-1"="MPA-8", "AI039-1"="MPA-7",
                      "AI054-1"="MPA-1", "AI055-1"="MPA-6", "AI062-1"="MPA-4", "AI075-1"="MPA-5",
                      "HD020-3"="HD-1", "HD036-1"="HD-2",
                      "HC037-1"="HD-3", "HC038-1"="HD-4", "HC039-1"="HD-5", "HC040-1"="HD-6",
                      "HC041-1"="HD-7","AI054-2"="MPA-1-post","AI075-2"="MPA-5-post")
s.int[["sample"]] <- Idents(s.int)
s.int$sample <- factor(s.int$sample,
                       levels= c("MPA-1","MPA-2","MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8",
                                 "HD-1", "HD-2", "HD-3", "HD-4", "HD-5", "HD-6", "HD-7"))

#DEG analysis using count data
counts <- GetAssayData(s.int, slot=c("counts"),assay="RNA")
genes.percent.expression <- rowMeans(counts>0)*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>15])
counts.sub <- counts[genes.filter,]
s.saver <- CreateSeuratObject(counts.sub)
s.saver@meta.data <- s.int@meta.data
AGG <- AggregateExpression(s.saver, group.by="sample", slot="counts")
write.csv(AGG[["RNA"]], "PBMC DEG.csv")

counts <- GetAssayData(s.int, slot=c("counts"),assay="RNA")
genes.percent.expression <- rowMeans(counts>0)*100
genes.filter <- names(genes.percent.expression[genes.percent.expression>0])
counts.sub <- counts[genes.filter,]
s.saver <- CreateSeuratObject(counts.sub)
s.saver@meta.data <- s.int@meta.data
AGG <- AggregateExpression(s.saver, group.by="sample", slot="counts")
write.csv(AGG[["RNA"]], "PBMC All Genes.csv")

#Dimplot for PBMC
s.int <- RunUMAP(s.int, reduction = 'ref.spca', dims = 1:30)
s.int$predicted.celltype.l2 <- factor(s.int$predicted.celltype.l2,
                       levels= c(
                                 "CD4 Naive", "CD4 TCM","CD4 TEM", "CD4 CTL",
                                 "CD4 Proliferating", "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "dnT",
                                 "gdT","Treg","MAIT", "B naive", "B intermediate", "B memory","Plasmablast",
                                 "NK", "NK_CD56bright","NK Proliferating", "CD14 Mono", "CD16 Mono", "cDC1",
                                 "cDC2", "ASDC", "pDC", "ILC", "HSPC"))
xp <- DimPlot(s.int, group.by="predicted.celltype.l2",reduction="umap",shuffle=T, pt.size=0.1)
xp
ggsave('Fig.1B.png',xp, width=7.5, height=5)

xp <- DimPlot(s.int, group.by="predicted.celltype.l2",split.by="disease",reduction="umap",shuffle=T, pt.size=0.1)
xp
ggsave('Fig.1C.png',xp, width=10, height=5)

Idents(s.int) <- "sample"
s.sub <- subset(s.int, downsample=2500, idents=c("MPA-1","MPA-2","MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8", 
                                                 "HD-1", "HD-2", "HD-3", "HD-4", "HD-5", "HD-6", "HD-7"))
pp <- DimPlot(s.sub,split.by="sample", ncol=8,group.by="predicted.celltype.l2",reduction="umap",shuffle=T, pt.size=0.2)
pp
ggsave('Fig. S1.png',pp,width=20,height=7.5)

saveRDS(s.int, "PBMC_ind.rds")

write.csv(table(s.int@meta.data[["sample"]], s.int@meta.data[["predicted.celltype.l2"]]), "PBMC ref population.csv")
