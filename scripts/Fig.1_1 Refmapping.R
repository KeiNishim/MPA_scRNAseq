## 0. Read in arguments, libraries, and custom scripts ################################
## 0.1 Read command line arguments ------------------------------
## expects the following arguments in order:
## - input file path (common path if processing files across multiple directories)
## - output file path
## - tab-delimited file containing table with info on h5 files to be processed. 
##   -> this file needs to contain disease, path, and sample columns
## - reference (default is specified below)
args = commandArgs(trailingOnly=TRUE)
args[1]=[xxx/input file path] 
args[2]=[xxx/output file path]
args[3]=[xxx/metadata.txt] 
print("0.1 Read in arguments: ")

## test if arguments given. if not, return an error
if (length(args) == 0){
    stop("Arguments must be supplied (input file path, output file path, file with files info, reference (default provided))", call.=FALSE)
} else if(length(args) == 3){
    # this is the default reference (Hao et al. PBMC)
    args[4] = [xxx/reference file path]
}
##--------------------------------------------------------------

## 0.2 Read in libraries and custom scripts --------------------------------
print("0.2 Read in libraries and custom scripts")
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

## function to read in data file and normalize
## recommendations from Satija lab is to integrate ADT and RNA assays in two different objects first
load_and_prep_ADT <- function(x){
    tmp.data <- Read10X_h5( x['path'])
    rownames(x=tmp.data[["Antibody Capture"]]) <- paste0("ADT-", rownames(tmp.data[["Antibody Capture"]]))
    tmp<- CreateSeuratObject(counts = tmp.data[["Antibody Capture"]], assay = "ADT")
    ## there migth be a better way to do this ###
    tmp <- AddMetaData(tmp, metadata = x['disease'], col.name="disease")
    tmp <- AddMetaData(tmp, metadata = x['sample'], col.name="sample")
    ###############################################
    tmp <- NormalizeData(tmp, normalization.method="CLR", margin = 2) ## normalization across cells 
    return(tmp)}

load_and_prep_RNA <- function(x){
    tmp.data <- Read10X_h5( x['path'])
    tmp <- CreateSeuratObject(counts = tmp.data[["Gene Expression"]], assay = "RNA")
    ##tmp <- CreateSeuratObject(counts = tmp.data[["Gene Expression"]], assay = "RNA")
    ## there migth be a better way to do this ###
    tmp <- AddMetaData(tmp, metadata = x['disease'], col.name="disease")
    tmp <- AddMetaData(tmp, metadata = x['sample'], col.name="sample")
    ###############################################
    tmp <- PercentageFeatureSet(tmp, pattern="^MT-", col.name="percent.mt")
    tmp <- SCTransform(tmp, method="glmGamPoi", vars.to.regress = c("percent.mt"))
    return(tmp)}
refmapping <- function(query_adt, query_rna, reference){
    ## default projects the predicted cell type label ("celltype.l2")
    ## from the reference to the query
    
    ## assumes that the query_adt and query_rna have been 
    ## processed (e.g. normalized by SCT) 
    
    query_adt[["RNA"]] <- query_rna[["SCT"]]
    DefaultAssay(query_adt) <- "RNA"
    ## find anchors 
    anchors <- FindTransferAnchors(
        reference = reference,
        query = query_adt,
        normalization.method = "SCT", ## modify if necessary
        reference.reduction = "spca", ## Seurat recommendations for CITE data
        ##dims = 1:50 )
        dims = 1:40 ) ## try with dims=40 2022-06-16 to match number of dimensions Nishimura-san used 
    query_adt <- MapQuery(
        anchorset = anchors,
        query = query_adt,
        reference = reference,
        ## transfer cell type labels (=celltype.l2)
        ## ***** if want to transfer other labels/info, modify here *****
        refdata = list(
            celltype.l2 = "celltype.l2"),
            ##umap_cluster = "umap_cluster"),
        reference.reduction = "spca",
        reduction.model = "wnn.umap" )
    return(query_adt)
}

do_refmapping <- function(query_adt_list, query_rna_list, reference){
    for (i in 1:length(query_adt_list)){
        query_adt_list[[i]] <- refmapping(query_adt_list[[i]], query_rna_list[[i]], reference)
    }
    return(query_adt_list)
}

## 0.3 Set data file path and check other paths ------------------------
print("0.3 Set data file path and check paths")
fileroot <- args[1]
print("Specify sample info")
all_files <- read.table(args[3], sep='\t', header=T)
print(all_files)
## check paths
check_filepaths <- function(df){
    filepath_needs_checking = FALSE
    for (i in seq(1,length(df$path))) {
        if(file.exists(df$path[i]) == FALSE){
            print(paste(i, df$sample[i], df$path[i], "path incorrect... CHECK"))
            filepath_needs_checking = TRUE
        }
    }
    return(filepath_needs_checking)
} 
try(if (check_filepaths(all_files)) stop("filepath needs checking"))
print("files checked")
#######################################################################

## 1. Process ADT and RNA data ############################
print("1 Process ADT and RNA data")
print("process ADT data")
all_samples_adt <- apply(all_files, 1, load_and_prep_ADT)
print("process RNA data")
all_samples_rna <- apply(all_files, 1, load_and_prep_RNA)
print("done processing both ADT and RNA data!")
## Save files as RDS
save_filepath = args[2]
date = Sys.Date()

## 2. Reference mapping ##############################################################
## 2.1 Read in PBMC reference ----------------------------------------
## from Hao et al. 2021 Cell #####
## in this case, the reference is already processed and integrated 
print("2.1 Read in reference")
reference <- LoadH5Seurat(args[4])
reference
print("finished reading reference")
## for visualization later, set cell cluster colors based on reference ----------------
set.seed(723451) ## for reproducibility of colors
## cell type colors
celltypeColors <- createPalette(length(unique(reference$celltype.l2)), c("#010101", "#ff0000"), M=10000)
celltypes <- as.vector(unique(reference$celltype.l2))
## rename some cell types
modcelltypes <- replace(celltypes, celltypes=='NK_CD56bright', 'NK-CD56bright')
modcelltypes <- modcelltypes[!modcelltypes %in% c('Doublet')]
## disease colors
disease_levels = unique(all_files$disease)
cohort_colors <- createPalette(length(disease_levels), c("#010101", "#ff0000"), M=10000)
## sample colors
samples <- unique(all_files$sample)
sample_colors <- createPalette(length(samples), c("#010101", "#ff0000"), M=10000)
##-------------------------------------------------------------------------------------

## 2.2 Do reference mapping -----------------------------
print("2.2 Do refmapping")
refmapped <- do_refmapping(all_samples_adt, all_samples_rna, reference) 
print("finished reference mapping")
refmapped
saveRDS(refmapped, file=paste0(save_filepath, 'refmapped_additional_samples.rds'))
print("saved refmapped")
##-------------------------------------------------------

## 2.3 Merge samples for ease of handling ------------------------------------------
print("2.3 Merge samples")
s.int<- merge(refmapped[[1]], refmapped[2:length(refmapped)], merge.dr = c("ref.umap", "ref.spca"))
print("finished merging")
saveRDS(s.int, file=paste0(save_filepath, 's.int.rds'))
write.csv(table(s.int@meta.data[["sample"]], s.int@meta.data[["predicted.celltype.l2"]]), "celltype.l2.csv")

