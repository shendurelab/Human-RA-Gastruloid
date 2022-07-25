#This script integrate MOCA E8.5-13.5 mouse embryo dataset (Cao et al., 2019, Qiu et al., 2022) with mouse gastruloids at 120hr (Van de Brink et al., 2020) and human RA gastruloids at 120hr
#Link to MOCA and mouse gastruloid datasets are provided in https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/
#Human gastrloid datasets from this paper are provided at National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) under accession numbers GSE208369

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

options(future.globals.maxSize= 50000 * 1024^2)
df_biomart <- read.table(paste0(data_dir,'ref/mart_export_human_mouse.txt'), header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)

downsample <- function(obj, subset){
    if(subset<1){
    subset=sample(1:ncol(obj),subset*ncol(obj))
    obj=obj[,subset]
        }else{
        subset=sample(1:ncol(obj),subset)
    obj=obj[,subset]}
    return(obj)
}

geo_dir = './geo/' #GSE208369
data_dir = './hGastruloid_website/public/' #https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/


MOCA_e8_5b <- readRDS(paste0(data_dir,'MOCA/obj_E8.5b_exp.rds'))
MOCA_e9_5 <- readRDS(paste0(data_dir,'MOCA/obj_E9.5_exp.rds'))
MOCA_e10_5 <- readRDS(paste0(data_dir,'MOCA/obj_E10.5_exp.rds'))
MOCA_e11_5 <- readRDS(paste0(data_dir,'MOCA/obj_E11.5_exp.rds'))
MOCA_e12_5 <- readRDS(paste0(data_dir,'MOCA/obj_E12.5_exp.rds'))
MOCA_e13_5 <- readRDS(paste0(data_dir,'MOCA/obj_E13.5_exp.rds'))


MOCA_e8_5b$timepoint="MOCA_8_5b"
MOCA_e9_5$timepoint="MOCA_9_5"
MOCA_e10_5$timepoint="MOCA_10_5"
MOCA_e11_5$timepoint="MOCA_11_5"
MOCA_e12_5$timepoint="MOCA_12_5"
MOCA_e13_5$timepoint="MOCA_13_5"

MOCA_e8_5b.anno <- readRDS(paste0(data_dir,'MOCA/obj_E8.5b_anno.rds'))
MOCA_e9_5.anno <- readRDS(paste0(data_dir,'MOCA/obj_E9.5_anno.rds'))
MOCA_e10_5.anno <- readRDS(paste0(data_dir,'MOCA/obj_E10.5_anno.rds'))
MOCA_e11_5.anno <- readRDS(paste0(data_dir,'MOCA/obj_E11.5_anno.rds'))
MOCA_e12_5.anno <- readRDS(paste0(data_dir,'MOCA/obj_E12.5_anno.rds'))
MOCA_e13_5.anno <- readRDS(paste0(data_dir,'MOCA/obj_E13.5_anno.rds'))


MOCA_e8_5b@meta.data$e8_5_Anno <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$Anno
MOCA_e9_5@meta.data$e9_5_Anno <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$Anno
MOCA_e10_5@meta.data$e10_5_Anno <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$Anno
MOCA_e11_5@meta.data$e11_5_Anno <- MOCA_e11_5.anno[colnames(MOCA_e11_5),]$Anno
MOCA_e12_5@meta.data$e12_5_Anno <- MOCA_e12_5.anno[colnames(MOCA_e12_5),]$Anno
MOCA_e13_5@meta.data$e13_5_Anno <- MOCA_e13_5.anno[colnames(MOCA_e13_5),]$Anno

MOCA_e8_5b@meta.data$celltype <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$celltype
MOCA_e9_5@meta.data$celltype <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$celltype
MOCA_e10_5@meta.data$celltype <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$celltype
MOCA_e11_5@meta.data$celltype<- MOCA_e11_5.anno[colnames(MOCA_e11_5),]$celltype
MOCA_e12_5@meta.data$celltype <- MOCA_e12_5.anno[colnames(MOCA_e12_5),]$celltype
MOCA_e13_5@meta.data$celltype <- MOCA_e13_5.anno[colnames(MOCA_e13_5),]$celltype

MOCA_e8_5b$batch='sci-b2'
MOCA_e9_5$batch='sci-b1'
MOCA_e10_5$batch='sci-b1'
MOCA_e11_5$batch='sci-b1'
MOCA_e12_5$batch='sci-b1'
MOCA_e13_5$batch='sci-b1'

obj_1 <- merge(x= downsample(MOCA_e8_5b,50000),
                     y=c(downsample(MOCA_e9_5,50000),
                         downsample(MOCA_e10_5,50000),
                         downsample(MOCA_e11_5,50000),
                         downsample(MOCA_e12_5,50000),
                         downsample(MOCA_e13_5,50000)),
          add.cell.ids = c("e8_5","e9_5","e10_5","e11_5","e12_5","e13_5"))

###################### MOCA lift-over 
project_name="MOCA_integrated"
MOCA_cell <- obj_1@meta.data
MOCA_gene <- obj_1[['RNA']][[]]
MOCA_mtx <- GetAssayData(object = obj_1, assay="RNA", slot = "counts")
df_mouse_human = df_biomart%>%filter(Mouse_ID %in% rownames(MOCA_mtx))
MOCA_mtx_sub = MOCA_mtx[as.vector(df_mouse_human$Mouse_ID),]
rownames(MOCA_mtx_sub) <- df_mouse_human$link
obj_1 <- CreateSeuratObject(counts = MOCA_mtx_sub,
project = project_name,
assay = "RNA",
meta.data = MOCA_cell)

##################### Gastruloids lift-over
obj_2 <- readRDS(paste0(geo_dir,'RA_hGas_120h.RDS'))   
hg_cells <- obj_2@meta.data
hg_gene <- obj_2[['RNA']][[]]
hg_mtx <- GetAssayData(object = obj_2, slot = "counts")
df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
rownames(hg_mtx_sub) = df_human_link$link
obj_2 = CreateSeuratObject(counts = hg_mtx_sub,
project = "hGastruloids",
assay = "RNA",
meta.data =  hg_cells)
obj_2$timepoint = "RA_120h"
obj_2$gastruloid1.cluster=obj_2$seurat_clusters
obj_2$batch='10x_hg_1'  #Set batch on each different samples

obj_3=readRDS(paste0(data_dir,'mouse_gastruloid/VDB_mGastruloids_120h_rna.RDS'))
obj_3 <- obj_3[,!is.na(obj_3$VDB_clusters)]
mg_cells <- obj_3@meta.data
mg_gene <- obj_3[['RNA']][[]]
mg_mtx <- GetAssayData(object = obj_3, slot = "counts")
df_mouse_human = df_biomart%>%filter(Mouse_name %in% rownames(mg_mtx))
mg_mtx_sub = mg_mtx[as.vector(df_mouse_human$Mouse_name),]
rownames(mg_mtx_sub) <- df_mouse_human$link
obj_3 = CreateSeuratObject(counts = mg_mtx_sub,
project = "mGastruloids",
assay = "RNA",
meta.data =  mg_cells)
obj_3$timepoint = "mg_120h"
obj_3$batch='10x_mg'

homologs=Reduce(intersect, c(rownames(obj_1),rownames(obj_2),rownames(obj_3)))
obj_1<-obj_1[homologs,]
obj_2<-obj_2[homologs,]
obj_3<-obj_3[homologs,]


obj = merge(x=obj_1, y=c(obj_2,obj_3),add.cell.ids= c('MOCA','"RA_120h','mg_120h'))
obj.list <- SplitObject(object = obj, split.by = "batch")
obj.list <- lapply(X = obj.list, FUN = SCTransform)
obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features, 
verbose = FALSE)
obj.list <- lapply(X = obj.list, FUN = RunPCA, features = obj.features)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
anchor.features = obj.features, dims = 1:30, reduction = "rpca", k.anchor = 20, verbose = FALSE)
obj.integrated <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", 
verbose = FALSE)
DefaultAssay(object = obj.integrated) <- "integrated"
obj <- obj.integrated
obj <- RunPCA(object = obj, verbose = FALSE)
obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.3)

pdf(file='~/MOCA85_135_mG_hG.pdf',width=30,height=20)
p1<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="Anno",pt.size=0.0001,
    label.size = 5, repel = TRUE, shuffle=TRUE, order=TRUE,raster=T)+NoLegend() 
p2<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="cell_type",pt.size=0.0001,
    label.size = 8, repel = TRUE, shuffle=FALSE, order=TRUE,raster=T)
p1+p2
dev.off()




                
