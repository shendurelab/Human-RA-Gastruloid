#This script convert scRNA-seq datasets of E6.5-13.5 mouse embryos (Pijuan-sala et al., 2019, Cao et al., 2019, Qiu et al., 2022), mouse gastruloids (van den Brink et al., 2020, Veenvilet et al., 2020) into pseudobulk sample matrices
#Link to MOCA,Pijuan,mouse gastruloids and somitoids datasets are provided in https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/
#Human gastrloid datasets from this paper are provided at National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) under accession numbers GSE208369
#Outputs are saved in the same folder

suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(readr)
        })

geo_dir = './geo/' #GSE208369
data_dir = './hGastruloid_website/public/' #https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/
df_biomart <- read.table(paste0(data_dir,'/ref/mart_export_human_mouse.txt'), header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)

pseudobulk_MOCA_by_embryo <- function(obj){
    obj_mtx <- obj[['RNA']]@counts
    df_mouse_human = df_biomart[df_biomart$Mouse_ID %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[as.vector(df_mouse_human$Mouse_ID),]
    rownames(obj_mtx_sub) = df_mouse_human$link
    
    obj_sum <- NULL
    for(id in obj$embryo_id%>%unique){
        cells <- colnames(obj[,obj$embryo_id==id])
        res <- rowSums(obj_mtx[,cells])
        obj_sum <- cbind(obj_sum,res)
    }
    
    colnames(obj_sum) <- obj$embryo_id%>%unique
    
    return(obj_sum)
}
pseudobulk_mouse <- function(obj){
    obj_mtx <- obj[['RNA']]@counts
    df_mouse_human = df_biomart[df_biomart$Mouse_name %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[as.vector(df_mouse_human$Mouse_name),]
    rownames(obj_mtx_sub) = df_mouse_human$link
    obj_sum <- rowSums(obj_mtx_sub)
    return(obj_sum)
}
pseudobulk_human <- function(obj,liftover=T){
    obj_mtx <- obj[['RNA']]@counts
    if(liftover){
    df_human_mouse = df_biomart[df_biomart$Human_name %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[df_human_mouse$Human_name,]
    rownames(obj_mtx_sub) = df_human_mouse$link
        }else{
        obj_mtx_sub <- obj_mtx
    }
    obj_sum <- rowSums(obj_mtx_sub)
    return(obj_sum)
}

#Pijuan
Pijuan <- readRDS(paste0(data_dir,"pijuan/Pijuan_merged.rds"))
Pijuan$celltype <- sapply(str_split(Pijuan$Anno,pattern=':', n=2),`[`,2)
select_celltype <- c('Epiblast','Primitive Streak','Nascent mesoderm','Anterior Primitive Streak',
                    'Mixed mesoderm','Haematoendothelial progenitors','Notochord','Blood progenitors',
                     'Somatic mesoderm','Splanchnic mesoderm','Caudal lateral epiblast','Caudal neurectoderm',
                     'Paraxial mesoderm A','Surface ectoderm','Gut',
                     'Paraxial mesoderm B','Paraxial mesoderm C','Primitive erythroid','Intermediate mesoderm',
                    'Endothelium','Cardiomyocytes','Spinal cord','Neural crest','Forebrain/Midbrain','Hindbrain',
                     'NMP','Somatic mesoderm A','Somatic mesoderm B','Urothelium','Placodal area','Epidermis')
Pijuan_sub <- Pijuan[,Pijuan$celltype %in% select_celltype]

#MOCA
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
MOCA_e8_5b@meta.data$Anno <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$Anno
MOCA_e9_5@meta.data$Anno <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$Anno
MOCA_e10_5@meta.data$Anno <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$Anno
MOCA_e11_5@meta.data$Anno <- MOCA_e11_5.anno[colnames(MOCA_e11_5),]$Anno
MOCA_e12_5@meta.data$Anno <- MOCA_e12_5.anno[colnames(MOCA_e12_5),]$Anno
MOCA_e13_5@meta.data$Anno <- MOCA_e13_5.anno[colnames(MOCA_e13_5),]$Anno
MOCA_e8_5b@meta.data$celltype <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$celltype
MOCA_e9_5@meta.data$celltype <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$celltype
MOCA_e10_5@meta.data$celltype <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$celltype
MOCA_e11_5@meta.data$celltype <- MOCA_e11_5.anno[colnames(MOCA_e11_5),]$celltype
MOCA_e12_5@meta.data$celltype <- MOCA_e12_5.anno[colnames(MOCA_e12_5),]$celltype
MOCA_e13_5@meta.data$celltype <- MOCA_e13_5.anno[colnames(MOCA_e13_5),]$celltype
MOCA_e8_5b$embryo_id <- gsub('E8.5','E8.5b',MOCA_e8_5b$RT_group)
MOCA_e9_5$embryo_id <- paste0('E9.5_',MOCA_e9_5$embryo_id)
MOCA_e10_5$embryo_id <- paste0('E10.5_',MOCA_e10_5$embryo_id)
MOCA_e11_5$embryo_id <- paste0('E11.5_',MOCA_e11_5$embryo_id)
MOCA_e12_5$embryo_id <- paste0('E12.5_',MOCA_e12_5$embryo_id)
MOCA_e13_5$embryo_id <- paste0('E13.5_',MOCA_e13_5$embryo_id)

all_celltype <- Reduce(union,c(
                               MOCA_e8_5b@meta.data$celltype%>%unique, MOCA_e9_5@meta.data$celltype%>%unique, 
      MOCA_e10_5@meta.data$celltype%>%unique, MOCA_e11_5@meta.data$celltype%>%unique,
      MOCA_e12_5@meta.data$celltype%>%unique,MOCA_e13_5@meta.data$celltype%>%unique))%>%unique
select_cells <- all_celltype[!all_celltype %in% c('Allantois','Extraembryonic visceral endoderm',
                                              'Amniochorionic mesoderm B','Amniochorionic mesoderm A',
                                             'Extraembryonic mesoderm','Mesenchymal stromal cells'
                                             )]
MOCA_e8_5b_sub <- MOCA_e8_5b[,MOCA_e8_5b$celltype %in% select_cells]
MOCA_e9_5_sub <- MOCA_e9_5[,MOCA_e9_5$celltype %in% select_cells]
MOCA_e10_5_sub <- MOCA_e10_5[,MOCA_e10_5$celltype %in% select_cells]
MOCA_e11_5_sub <- MOCA_e11_5[,MOCA_e11_5$celltype %in% select_cells]
MOCA_e12_5_sub <- MOCA_e12_5[,MOCA_e12_5$celltype %in% select_cells]
MOCA_e13_5_sub <- MOCA_e13_5[,MOCA_e13_5$celltype %in% select_cells]

#Gastruloids
#VDB gastruloid
mGastruloids_120h_rna <- readRDS(paste0(data_dir,'mouse_gastruloid/VDB_mGastruloids_120h_rna.RDS'))
Idents(mGastruloids_120h_rna)=mGastruloids_120h_rna$VDB_celltype

#TLS gastruloid
TLS_obj <- readRDS(file=paste0(data_dir,'TLS/TLS_timecourse_rna.RDS'))
TLS_96hr <- TLS_obj[,TLS_obj$timepoint=='96h']
TLS_108hr <- TLS_obj[,TLS_obj$timepoint=='108h']
TLS_120hr<- TLS_obj[,TLS_obj$timepoint=='120h']
Idents(TLS_96hr)<-TLS_96hr$cell_state
Idents(TLS_108hr)<-TLS_108hr$cell_state
Idents(TLS_120hr)<-TLS_120hr$cell_state
TLS_96hr$orig.ident='TLS_96h'
TLS_108hr$orig.ident='TLS_108h'
TLS_120hr$orig.ident='TLS_120h'
# TLS_meta <- read.table(paste0(data_dir,'TLS/TLS_meta_data.tsv'),
#                       sep='\t',header=T)
# TLS_meta$timepoint=str_split_fixed(TLS_meta$TP,'_',2)[,2]
# TLS_meta$cell_id=paste0(TLS_meta$timepoint,'r_',TLS_meta$BC,'-1')
# TLS_obj$cell_id=row.names(TLS_obj)
# df_TLS=TLS_obj@meta.data
# df_TLS=left_join(df_TLS,TLS_meta,by='cell_id')
# TLS_obj@meta.data=df_TLS
# TLS_obj<-TLS_obj[,TLS_obj$cell_state != 'Unknown']
# TLS_96hr <- TLS_obj[,TLS_obj$timepoint=='96h']
# TLS_108hr <- TLS_obj[,TLS_obj$timepoint=='108h']
# TLS_120hr<- TLS_obj[,TLS_obj$timepoint=='120h']
# Idents(TLS_96hr)<-TLS_96hr$cell_state
# Idents(TLS_108hr)<-TLS_108hr$cell_state
# Idents(TLS_120hr)<-TLS_120hr$cell_state

#Human gastruloid
#hG_0h_rna <- readRDS(paste0(geo_dir,'hGas_0h.RDS'))
hG_24h_rna<- readRDS(paste0(geo_dir,'hGas_24h.RDS'))
hG_48h_rna <- readRDS(paste0(geo_dir,'hGas_48h.RDS'))
hG_72h_rna <- readRDS(paste0(geo_dir,'hGas_72h.RDS'))
hG_96h_rna <- readRDS(paste0(geo_dir,'hGas_96h.RDS'))
#hG_RA_96h_rna <- readRDS(paste0(geo_dir,'RA_hGas_96h.RDS'))
hG_RA_120h_rna = readRDS(paste0(geo_dir,'RA_hGas_120h.RDS'))

#Somitoids
Somitoids_168h_CHIR8_10<-readRDS(paste0(data_dir,'human_somitoids/human_somitoids_b1.RDS'))
Somitoids_168h_CHIR5_7<- readRDS(paste0(data_dir,'human_somitoids/human_somitoids_MULTI-seq_b2.RDS'))

#Pseudobulk samples
#MOCA 8.5-13.5
MOCA_list <- list(
                  MOCA_e8_5b=pseudobulk_MOCA_by_embryo(MOCA_e8_5b_sub),
                  MOCA_e9_5=pseudobulk_MOCA_by_embryo(MOCA_e9_5_sub),
                  MOCA_e10_5=pseudobulk_MOCA_by_embryo(MOCA_e10_5_sub),
                  MOCA_e11_5=pseudobulk_MOCA_by_embryo(MOCA_e11_5_sub),
                  MOCA_e12_5=pseudobulk_MOCA_by_embryo(MOCA_e12_5_sub),
                  MOCA_e13_5=pseudobulk_MOCA_by_embryo(MOCA_e13_5_sub)
                  )
gene_list <- NULL
for(i in 1:length(MOCA_list)){
    gene_list[[i]] <- rownames(MOCA_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
data_list_subset<-MOCA_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes,]
}
data_mtx <- do.call(cbind,data_list_subset)  
dim(data_mtx)
saveRDS(Matrix(data_mtx,sparse=TRUE), './outputs/MOCA_mtx_noExE.RDS')

#Human RA gastruloid 96hr, 120hr (mito-depletion vs non mito-depletion)
data_list <- list(
                  hG_RA_120h_rna=pseudobulk_human(hG_RA_120h_rna,liftover=T), 
#                  hG_RA_96h_rna=pseudobulk_human(hG_RA_96h_rna,liftover=T),
                  hG_96h_rna=pseudobulk_human(hG_96h_rna,liftover=T), 
                  hG_72h_rna=pseudobulk_human(hG_72h_rna,liftover=T),
                  hG_48h_rna=pseudobulk_human(hG_48h_rna,liftover=T),
                  hG_24h_rna=pseudobulk_human(hG_24h_rna,liftover=T)
#                  hG_0h_rna=pseudobulk_human(hG_0h_rna,liftover=T)
                    )

#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
data_list_subset<-data_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list_subset)
dim(data_mtx)
saveRDS(Matrix(data_mtx,sparse=TRUE), './outputs/query_mtx_hG_noliftover.RDS')

#TLS mouse gastruloid
# #VDB mouse gastruloid
data_list <- list(
                  mG_120h_rna=pseudobulk_mouse(mGastruloids_120h_rna),
                  TLS_96hr=pseudobulk_mouse(TLS_96hr),TLS_108hr=pseudobulk_mouse(TLS_108hr),
                  TLS_120hr=pseudobulk_mouse(TLS_120hr)
                    )
#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
data_list_subset<-data_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list_subset)
dim(data_mtx)
saveRDS(Matrix(data_mtx,sparse=TRUE), './outputs/query_mtx_mG_noliftover.RDS')

#Somitoids
data_list <- list(
                  Somitoid_168h_CHIR5_7_rna = pseudobulk_human(Somitoids_168h_CHIR5_7,liftover=F),
                  Somitoid_168h_CHIR8_10_rna = pseudobulk_human(Somitoids_168h_CHIR8_10,liftover=F)
                    )
#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
data_list_subset<-data_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list_subset)
dim(data_mtx)
saveRDS(Matrix(data_mtx,sparse=TRUE), './outputs/query_mtx_somitoid_noliftover.RDS')