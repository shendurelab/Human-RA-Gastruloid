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

fig.dir = '~/shendure_lab/hGastruloid_figures/' ## change the output dir for figures

downsample <- function(obj, subset){
    if(subset<1){
    subset=sample(1:ncol(obj),subset*ncol(obj))
    obj=obj[,subset]
        }else{
        subset=sample(1:ncol(obj),subset)
    obj=obj[,subset]}
    return(obj)
}


MOCA_e8_5b <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E8.5b_exp.rds')
MOCA_e9_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E9.5_exp.rds')
MOCA_e10_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E10.5_exp.rds')
MOCA_e11_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E11.5_exp.rds')
MOCA_e12_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E12.5_exp.rds')
MOCA_e13_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E13.5_exp.rds')


MOCA_e8_5b$timepoint="MOCA_8_5b"
MOCA_e9_5$timepoint="MOCA_9_5"
MOCA_e10_5$timepoint="MOCA_10_5"
MOCA_e11_5$timepoint="MOCA_11_5"
MOCA_e12_5$timepoint="MOCA_12_5"
MOCA_e13_5$timepoint="MOCA_13_5"

MOCA_e8_5b.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E8.5b_anno.rds")
MOCA_e9_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E9.5_anno.rds")
MOCA_e10_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E10.5_anno.rds")
MOCA_e11_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E11.5_anno.rds")
MOCA_e12_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E12.5_anno.rds")
MOCA_e13_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E13.5_anno.rds")


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

options(future.globals.maxSize= 50000 * 1024^2)
df_biomart <- read.table("/net/shendure/vol1/home/weiy666/shendure_lab/sci-fate/data/20201211_sci-plex-fate/mart_export_human_mouse.txt", header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)

###################### MOCA lift-over 
# You can use my MOCA reference or make your own reference with the following code
# project_name="MOCA_integrated"
# MOCA_cell <- obj_1@meta.data
# MOCA_gene <- obj_1[['RNA']][[]]
# MOCA_mtx <- GetAssayData(object = obj_1, assay="RNA", slot = "counts")
# df_mouse_human = df_biomart%>%filter(Mouse_ID %in% rownames(MOCA_mtx))
# MOCA_mtx_sub = MOCA_mtx[as.vector(df_mouse_human$Mouse_ID),]
# rownames(MOCA_mtx_sub) <- df_mouse_human$link
# obj_1 <- CreateSeuratObject(counts = MOCA_mtx_sub,
# project = project_name,
# assay = "RNA",
# meta.data = MOCA_cell)

# saveRDS(obj_1,'/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/integration/MOCA_E85_E135_ref.RDS')
obj_1 = readRDS('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/integration/MOCA_E85_E135_ref.RDS') 

##################### Gastruloids lift-over
processed.dir <- '/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/hGastruloids/' #Choose your gastruloid directory
obj_2 <- readRDS(paste0(processed.dir,'hGastruloids_RA_120h_rna_processed.RDS'))
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
obj_2$batch='10x_hg_1'  #We are going to correct by batch, so make sure to set batch on each different samples

# obj_4 <- readRDS(paste0(processed.dir,'hGastruloids_RA_120h_rna_processed_b2.RDS'))
# hg_cells <- obj_4@meta.data
# hg_gene <- obj_4[['RNA']][[]]
# hg_mtx <- GetAssayData(object = obj_4, slot = "counts")
# df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
# hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
# rownames(hg_mtx_sub) = df_human_link$link
# obj_4 = CreateSeuratObject(counts = hg_mtx_sub,
# project = "hGastruloids2",
# assay = "RNA",
# meta.data =  hg_cells)
# obj_4$timepoint = "RA_120h_2"
# obj_4$gastruloid2.cluster=obj_4$seurat_clusters
# obj_4$batch='10x_hg_2'

# obj_4 <- readRDS(paste0(processed.dir,'hGastruloids_RA_120h_rna_processed_b2.RDS'))
# hg_cells <- obj_4@meta.data
# hg_gene <- obj_4[['RNA']][[]]
# hg_mtx <- GetAssayData(object = obj_4, slot = "counts")
# df_human_link = df_biomart[df_biomart$Human_name %in% rownames(hg_mtx),]
# hg_mtx_sub = hg_mtx[df_human_link$Human_name,]
# rownames(hg_mtx_sub) = df_human_link$link
# obj_4 = CreateSeuratObject(counts = hg_mtx_sub,
# project = "hGastruloids2",
# assay = "RNA",
# meta.data =  hg_cells)
# obj_4$timepoint = "RA_120h_2"
# obj_4$gastruloid2.cluster=obj_4$seurat_clusters
# obj_4$batch='10x_hg_2'

# mouse.dir <- '/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/mGastruloids/'
# obj_3=readRDS(paste0(mouse.dir,'mGastruloids_120h_rna.RDS'))
# obj_3 <- obj_3[,!is.na(obj_3$VDB_clusters)]
# mg_cells <- obj_3@meta.data
# mg_gene <- obj_3[['RNA']][[]]
# mg_mtx <- GetAssayData(object = obj_3, slot = "counts")
# df_mouse_human = df_biomart%>%filter(Mouse_name %in% rownames(mg_mtx))
# mg_mtx_sub = mg_mtx[as.vector(df_mouse_human$Mouse_name),]
# rownames(mg_mtx_sub) <- df_mouse_human$link
# obj_3 = CreateSeuratObject(counts = mg_mtx_sub,
# project = "mGastruloids",
# assay = "RNA",
# meta.data =  mg_cells)
# obj_3$timepoint = "mg_120h"
# obj_3$batch='10x_mg'


# TLS_obj=readRDS(paste0(mouse.dir,'TLS_timecourse_rna.RDS'))
# TLS_obj$batch <- TLS_obj$orig.ident
# TLS_cell <- TLS_obj[[]]
# TLS_gene <- TLS_obj[['RNA']][[]]
# TLS_mtx <- GetAssayData(object = TLS_obj, slot = "counts")
# df_mouse_human = df_biomart[df_biomart$Mouse_name %in% row.names(TLS_obj),]
# TLS_mtx_sub = TLS_mtx[as.vector(df_mouse_human$Mouse_name),]
# rownames(TLS_mtx_sub) = df_mouse_human$link
# TLS_humanized = CreateSeuratObject(counts = TLS_mtx_sub,
# project = "TLS_humanized",
# assay = "RNA",
# meta.data =  TLS_cell)
# TLS_humanized$orig.ident = "TLS_humanized"
# obj_5=TLS_humanized



#homologs=intersect(rownames(obj_1),rownames(obj_3))
homologs=Reduce(intersect, c(rownames(obj_1),rownames(obj_2))#,rownames(obj_3),rownames(obj_4)))
obj_1<-obj_1[homologs,]
obj_2<-obj_2[homologs,]
# obj_4<-obj_4[homologs,]
# obj_3<-obj_3[homologs,]
#obj_5<-obj_5[homologs,]


obj = merge(x=obj_1, y=c(obj_2,obj_4,obj_3))
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

saveRDS(obj, file='/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/integration/MOCA_8_13_mG_hG2_Integration.RDS')

pdf(file=paste0(fig.dir,'MOCA85_135_mG_hG2_knn.pdf'),width=30,height=20)
p1<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="knn.timepoint",pt.size=0.0001,
    label.size = 5, repel = TRUE, shuffle=TRUE, order=TRUE,raster=T)+NoLegend() 
p2<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="hG.merged.cluster",pt.size=0.0001,
    label.size = 8, repel = TRUE, shuffle=FALSE, order=TRUE,raster=T)#+NoLegend() 
p1+p2
dev.off()
                

pdf(file=paste0(fig.dir,'MOCA85_135_mG_hG2_knn.pdf'),width=30,height=20)
p1<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="knn.timepoint",pt.size=0.0001,
    label.size = 5, repel = TRUE, shuffle=TRUE, order=TRUE,raster=T)+NoLegend() 
p2<-DimPlot(obj, reduction = "umap",  label = TRUE, group.by="hG.merged.cluster",pt.size=0.0001,
    label.size = 8, repel = TRUE, shuffle=FALSE, order=TRUE,raster=T)#+NoLegend() 
p1+p2
dev.off()

# pdf(file=paste0(fig.dir,'MOCA85_105_hG2_MOCA_Anno_hG2_FP_expr.pdf'),width=30,height=10)
# do.call(patchwork::wrap_plots, lapply(df_biomart[df_biomart$Human_name %in% c('FOXA1','FOXA2','CORIN'),]$link, function(x) {FeaturePlot(obj[,obj$timepoint %in% c('MOCA_8_5b','MOCA_9_5','MOCA_10_5')], features=x,pt.size = 0.1, order=TRUE, min.cutoff=10, combine = TRUE) +labs(title=df_biomart[df_biomart$link==x,]$Human_name)})) -> plot
# print(plot)
# dev.off()

# library(plotly)
# library(htmlwidgets)

# df = cbind(Embeddings(object = downsample(obj,100000), reduction = "umap"), downsample(obj,100000)@meta.data)
# fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, size = I(3), color = ~celltype)
# fig = fig %>% layout(
#   scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
#                yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
# saveWidget(fig, paste0('~/nobackups/MOCA_mG_hG_umap2d.html'))

# library(plotly)
# library(htmlwidgets)

# df = cbind(Embeddings(object = obj, reduction = "umap"), obj@meta.data)
# fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, size = I(1), color = ~celltype)
# fig = fig %>% layout(
#   scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
#                yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
# saveWidget(fig, paste0(fig.dir,'MOCA85_105_hG_120hr_cell_type.html'))