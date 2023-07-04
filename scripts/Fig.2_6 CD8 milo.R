library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#Read CD8_ind as s.sub
pbmc_small_sce <- as.SingleCellExperiment(s.sub)
traj_milo <- Milo(pbmc_small_sce)
plotUMAP(traj_milo, colour_by="class", point_size=0.1)
traj_milo <- buildGraph(traj_milo, k = 10, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

saveRDS(traj_milo, "CD8_milo.rds")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "class")]
traj_design <- distinct(traj_design)
age <- data.frame(c(79,43,74,52,79,80,69,70,64,62,70,62,33, 27,36))
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
traj_design$class <- c("2_MPA","2_MPA","2_MPA","2_MPA","2_MPA","2_MPA","2_MPA","2_MPA","1_HD","1_HD","1_HD","1_HD","1_HD","1_HD","1_HD")
da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+class, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR)))+ 
  geom_point() +
  geom_hline(yintercept = 1)

## Plot single-cell UMAP
plotReducedDim(traj_milo, dimred = "UMAP", colour_by="class", text_size = 3, point_size=0.1)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,size_range=c(0.8,3), layout="UMAP",alpha=1)+
  scale_fill_gradient2(low='blue', mid='white', high="red")
nh_graph_pl
ggsave('Fig.2G.png', nh_graph_pl +plot_layout(guides="collect"), width=7.5, height=7.5)

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")
da_results$seurat_clusters <- factor(da_results$seurat_clusters,
                                     levels= c("CD8 T_KIR", "CD8 T_CTL","CD8 T_EM", "CD8 T_CM", "CD8 T_Naive"))
pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha=1)+geom_boxplot(outlier.shape = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
pDa
ggsave('Fig.2H.png', pDa, width=7.5, height=7.5)

median(da_results[da_results$seurat_clusters=="CD8 T_Naive", "logFC"])
median(da_results[da_results$seurat_clusters=="CD8 T_CM", "logFC"])
median(da_results[da_results$seurat_clusters=="CD8 T_EM", "logFC"])
median(da_results[da_results$seurat_clusters=="CD8 T_CTL", "logFC"])
median(da_results[da_results$seurat_clusters=="CD8 T_KIR", "logFC"])
