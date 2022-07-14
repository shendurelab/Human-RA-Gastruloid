#This script integrate scRNA-seq datasets of MOCA E8.5-13.5, RA gastruloids 120hr and mouse gastruloids 120hr (Van de Brink., 2021)

suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(readr)
library(reticulate)
        })

downsample <- function(obj, subset){
    if(subset<1){
    subset=sample(1:ncol(obj),subset*ncol(obj))
    obj=obj[,subset]
        }else{
        subset=sample(1:ncol(obj),subset)
    obj=obj[,subset]}
    return(obj)
}

pijuan <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/EB_analysis/Pijuan.rds')

obj_1=pijuan
options(future.globals.maxSize= 30000 * 1024^2)
df_biomart <- read.table("/net/shendure/vol1/home/weiy666/shendure_lab/sci-fate/data/20201211_sci-plex-fate/mart_export_human_mouse.txt", header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)
project_name="pijuan_integrated"
pijuan_cell <- obj_1@meta.data
pijuan_gene <- obj_1[['RNA']][[]]
pijuan_mtx <- GetAssayData(object = pijuan, assay="RNA", slot = "counts")
df_mouse_human = df_biomart%>%filter(Mouse_ID %in% rownames(pijuan_mtx))
pijuan_mtx_sub = pijuan_mtx[as.vector(df_mouse_human$Mouse_ID),]
rownames(pijuan_mtx_sub) <- df_mouse_human$link
obj_1 <- CreateSeuratObject(counts = pijuan_mtx_sub,
project = project_name,
assay = "RNA",
meta.data = pijuan_cell)  
obj_1$batch="pijuan_mouse" 

processed.dir <- "/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/hGastruloids/"

obj_3<-readRDS(paste0(processed.dir,'hGastruloid_24h_merged_ccr_annotated.RDS'))
obj_3$cluster.24h=Idents(obj_3)
hg_cells <- obj_3@meta.data
hg_gene <- obj_3[['RNA']][[]]
hg_mtx <- GetAssayData(object = obj_3, slot = "counts")
df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
rownames(hg_mtx_sub) = df_human_link$link
obj_3 = CreateSeuratObject(counts = hg_mtx_sub,
project = "hG_24",
assay = "RNA",
meta.data =  hg_cells)
obj_3$batch='b2'

obj_4<-readRDS(paste0(processed.dir,'hGastruloid_48h_merged_ccr_annotated.RDS'))
obj_4$cluster.48h=Idents(obj_4)
hg_cells <- obj_4@meta.data
hg_gene <- obj_4[['RNA']][[]]
hg_mtx <- GetAssayData(object = obj_4, slot = "counts")
df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
rownames(hg_mtx_sub) = df_human_link$link
obj_4 = CreateSeuratObject(counts = hg_mtx_sub,
project = "hG_48",
assay = "RNA",
meta.data =  hg_cells)
obj_4$batch='b2'

obj_5<-readRDS(paste0(processed.dir,'hGastruloid_72h_merged_ccr_annotated.RDS'))
obj_5$cluster.72h=Idents(obj_5)
hg_cells <- obj_5@meta.data
hg_gene <- obj_5[['RNA']][[]]
hg_mtx <- GetAssayData(object = obj_5, slot = "counts")
df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
rownames(hg_mtx_sub) = df_human_link$link
obj_5 = CreateSeuratObject(counts = hg_mtx_sub,
project = "hG_72",
assay = "RNA",
meta.data =  hg_cells)

obj_6<-readRDS(paste0(processed.dir,'hGastruloid_96h_merged_ccr_annotated.RDS'))
obj_6$cluster.96h=Idents(obj_6)
hg_cells <- obj_6@meta.data
hg_gene <- obj_6[['RNA']][[]]
hg_mtx <- GetAssayData(object = obj_6, slot = "counts")
df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
rownames(hg_mtx_sub) = df_human_link$link
obj_6 = CreateSeuratObject(counts = hg_mtx_sub,
project = "hG_96",
assay = "RNA",
meta.data =  hg_cells)
obj_6$batch='b1'


homologs=Reduce(intersect,c(rownames(obj_1),rownames(obj_3),rownames(obj_4),row.names(obj_5),row.names(obj_6)))
obj_1<-obj_1[homologs,]
obj_3<-obj_3[homologs,]
obj_4<-obj_4[homologs,]
obj_5<-obj_5[homologs,]
obj_6<-obj_6[homologs,]


obj = merge(x=obj_1, y=c(obj_3,obj_4,obj_5,obj_6), add.cell.ids= c('MOCA','24h','48h','72h','96h'))
obj.list <- SplitObject(object = obj, split.by = "batch")
for(i in 1:length(obj.list)){
        obj.list[[i]] <- SCTransform(obj.list[[i]],verbose=FALSE)
        }
obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features, 
verbose = FALSE)
obj.list <- lapply(X = obj.list, FUN = RunPCA, features = obj.features)

obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
anchor.features = obj.features, dims = 1:30,reduction = "rpca", k.anchor = 20, verbose = FALSE)

obj.integrated <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", 
verbose = FALSE)
DefaultAssay(object = obj.integrated) <- "integrated"
obj <- obj.integrated
obj <- RunPCA(object = obj, verbose = FALSE)
obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.3)





