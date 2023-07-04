#For plotting PBMC with ISG/monocyte signature###

#Read PBMC_ind as s.int
Idents(s.int) <- "sample"
s.sub <- subset(s.int, idents=c("MPA-5", "MPA-1", "MPA-2", "MPA-7", "MPA-8","MPA-6", "MPA-4", "MPA-3"),downsample=2870)
p <-DoHeatmap(s.sub, group.by="sample",
              disp.min=-2.5,disp.max=2.5,raster=F,size=5,
              features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                         "SIGLEC1","IFI27", "IFITM1", "ISG20",
                         "CD14","MS4A6A","FOS","DUSP1","IER2","S100A9","VCAN","CSF3R",
                         "LYZ","CEBPD","NAMPT",
                         "S100A8","MNDA","FCN1","S100A12",
                         "CD8A","CD8B", "CCL5","NKG7","GZMB","CST7","GZMH","GNLY","PRF1","KLRD1",
                         "NKG7"
),draw.lines=T, slot="scale.data")
p
ggsave('Fig.3B.png', p, width=12, height=8)

#Enrichment score analysis
features <- list(c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                   "IFITM1", "ISG20"))
features2 <- list(c("CD14","MS4A6A","FOS","DUSP1","IER2","S100A9","VCAN","CSF3R",
                    "LYZ","CEBPD","NAMPT",
                    "S100A8","MNDA","FCN1","S100A12"))
features3 <- list(c("CD8A","CD8B", "CCL5","NKG7","GZMB","CST7","GZMH","GNLY","PRF1","KLRD1",
                    "NKG7"))
s.sub <- AddModuleScore(s.sub, features=features,ctrl=1,name="IFN.score")
s.sub <- AddModuleScore(s.sub, features=features2,ctrl=1,name="mono.score")
s.sub <- AddModuleScore(s.sub, features=features3,ctrl=1,name="CD8.score")
p <- DotPlot(s.sub,scale=T,cols=c("lightgrey","red"), features=c("CD8.score1","mono.score1", "IFN.score1"),
            col.max=2, col.min=-1, dot.scale=7,group.by="sample")+
  coord_flip()+ theme(axis.text.x = element_text(size = 9)) +
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
p
ggsave('Fig.3C.png', p, width=7, height=3.2)
write.csv(p[["data"]], "modulescore.csv")

#For plotting subpopulation with ISG signature###
DefaultAssay(s.int) <- "RNA"
Idents(s.int) <- "predicted.celltype.l2"
s.sub2 <- subset(s.int, ident=c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2"))
s.sub2 <- ScaleData(s.sub2)
Idents(s.sub2) <- "sample"
p1 <-DoHeatmap(subset(s.sub2,idents=c("MPA-1", "MPA-2", "MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8"),downsample=100), 
              group.by="sample",disp.min=-2.5,disp.max=2.5,raster=F,size=5,
              features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                         "IFITM1", "ISG20"
              ),lines.width=T,draw.lines=F, slot="scale.data")
AGG <- AverageExpression(s.sub2, features=features[[1]], group.by="sample")
write.csv(AGG[["RNA"]],"FigS7-1.csv")
s.sub2 <- subset(s.int, ident=c("CD8 Naive", "CD8 Proliferating","CD8 TCM","CD8 TEM"))
s.sub2 <- ScaleData(s.sub2)
Idents(s.sub2) <- "sample"
p2 <-DoHeatmap(subset(s.sub2,idents=c("MPA-1", "MPA-2", "MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8"),downsample=100), 
               group.by="sample",disp.min=-2.5,disp.max=2.5,raster=F,size=5,
               features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                          "IFITM1", "ISG20"
               ),lines.width=T,draw.lines=F, slot="scale.data")
AGG <- AverageExpression(s.sub2, features=features[[1]], group.by="sample")
write.csv(AGG[["RNA"]],"FigS7-2.csv")
s.sub2 <- subset(s.int, ident=c("B intermediate", "B memory","B naive","Plasmablast"))
s.sub2 <- ScaleData(s.sub2)
Idents(s.sub2) <- "sample"
p3 <-DoHeatmap(subset(s.sub2,idents=c("MPA-1", "MPA-2", "MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8"),downsample=100), 
               group.by="sample",disp.min=-2.5,disp.max=2.5,raster=F,size=5,
               features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                          "IFITM1", "ISG20"
               ),lines.width=T,draw.lines=F, slot="scale.data")
AGG <- AverageExpression(s.sub2, features=features[[1]], group.by="sample")
write.csv(AGG[["RNA"]],"FigS7-3.csv")
s.sub2 <- subset(s.int, ident=c("CD4 CTL", "CD4 Naive","CD4 TCM","CD4 TEM"))
s.sub2 <- ScaleData(s.sub2)
Idents(s.sub2) <- "sample"
p4 <-DoHeatmap(subset(s.sub2,idents=c("MPA-1", "MPA-2", "MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8"),downsample=100), 
               group.by="sample",disp.min=-2.5,disp.max=2.5,raster=F,size=5,
               features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                          "IFITM1", "ISG20"
               ),lines.width=T,draw.lines=F, slot="scale.data")
AGG <- AverageExpression(s.sub2, features=features[[1]], group.by="sample")
write.csv(AGG[["RNA"]],"FigS7-4.csv")
s.sub2 <- subset(s.int, ident=c("NK", "NK_CD56bright"))
s.sub2 <- ScaleData(s.sub2)
Idents(s.sub2) <- "sample"
p5 <-DoHeatmap(subset(s.sub2,idents=c("MPA-1", "MPA-2", "MPA-3", "MPA-4", "MPA-5", "MPA-6", "MPA-7", "MPA-8"),downsample=50), 
               group.by="sample",disp.min=-2.5,disp.max=2.5,raster=F,size=5,
               features=c("IFI44L",  "IFI44", "ISG15","XAF1", "IFITM2", "IFI6","MX1","IFIT1", "IFIT2", "IFIT3", 
                          "IFITM1", "ISG20"
               ),lines.width=T,draw.lines=F, slot="scale.data")
AGG <- AverageExpression(s.sub2, features=features[[1]], group.by="sample")
write.csv(AGG[["RNA"]],"FigS7-5.csv")
p1
ggsave('Fig.S7-1.png', p1, width=6, height=4)
p2
ggsave('Fig.S7-2.png', p2, width=6, height=4)
p3
ggsave('Fig.S7-3.png', p3, width=6, height=4)
p4
ggsave('Fig.S7-4.png', p4, width=6, height=4)
p5
ggsave('Fig.S7ap.png', p5, width=6, height=6)
