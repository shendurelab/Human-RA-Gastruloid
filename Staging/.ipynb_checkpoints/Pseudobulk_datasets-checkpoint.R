#This script convert scRNA-seq datasets into pseudobulk sample matrices
suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(readr)
        })

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
Pijuan <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/EB_analysis/Pijuan.rds")
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
rna_dir_path=paste0('/net/shendure/vol10/projects/cxqiu/nobackup/data/gastruloid_mouse_Brink/')
rna_mtx=Read10X(data.dir=rna_dir_path)
mGastruloids_120h_rna_total=CreateSeuratObject(rna_mtx,
                                        min.cells = 0.01*dim(rna_mtx)[2],
                                        min.features = 200)
annotat_dir="/net/shendure/vol1/home/weiy666/shendure_lab/sci-fate/data/VDB_mgastruloid_results"
annotat_cell=readr::read_csv(paste0(annotat_dir,"/VDB_cell_meta.csv"),col_names=TRUE)
annotat_cell=annotat_cell%>%data.frame()
rownames(annotat_cell)=annotat_cell$cell.barcode
select = intersect(colnames(mGastruloids_120h_rna_total),annotat_cell[annotat_cell$batch=='10x',]%>%rownames)
mGastruloids_120h_rna = mGastruloids_120h_rna_total[,select]
mGastruloids_120h_rna=subset(mGastruloids_120h_rna_total,subset= nCount_RNA >1000 & 
           nFeature_RNA > 700 & 
           nFeature_RNA < 8000 & 
           nCount_RNA < 40000)
mGastruloids_120h_rna <- readRDS('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/mGastruloids/mGastruloids_120h_rna.RDS')
Idents(mGastruloids_120h_rna)=mGastruloids_120h_rna$VDB_celltype

#TLS gastruloid
load_TLS <- function(timepoint){
    rna_dir_path=paste0('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/raw_data/TLS/TLS_',timepoint,'hr/filtered_feature_bc_matrix/')
    rna_mtx=Read10X(data.dir=rna_dir_path)
    obj_TLS=CreateSeuratObject(rna_mtx,
                              min.cells = 3,
                              min.features = 200)
    obj_TLS[["percent.mt"]] <- PercentageFeatureSet(obj_TLS, pattern = "^mt-")
    obj_TLS$orig.ident="mouse-TLS"
    obj_TLS=subset(obj_TLS, subset= nCount_RNA >10000 & 
           nFeature_RNA > 3000 & 
           nCount_RNA < 40000 &
           percent.mt < 10)
    
    obj_TLS <- NormalizeData(obj_TLS, verbose = FALSE)
    obj_TLS <- CellCycleScoring(obj_TLS, 
                                 s.features = s.genes.mouse, 
                                 g2m.features = g2m.genes.mouse, 
                                 set.ident = TRUE)
    
    return(obj_TLS)
}
TLS_96hr <- load_TLS('96')
TLS_108hr <- load_TLS('108')
TLS_120hr <- load_TLS('120')
TLS_96hr$orig.ident='TLS_96h'
TLS_108hr$orig.ident='TLS_108h'
TLS_120hr$orig.ident='TLS_120h'
TLS_obj <- merge(x=TLS_96hr,y=c(TLS_108hr,TLS_120hr),add.cell.ids =c('96hr','108hr','120hr'))
TLS_meta <- read.table('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/raw_data/TLS/TLS_meta_data.tsv',
                      sep='\t',header=T)
TLS_meta$timepoint=str_split_fixed(TLS_meta$TP,'_',2)[,2]
TLS_meta$cell_id=paste0(TLS_meta$timepoint,'r_',TLS_meta$BC,'-1')
TLS_obj$cell_id=row.names(TLS_obj)
df_TLS=TLS_obj@meta.data
df_TLS=left_join(df_TLS,TLS_meta,by='cell_id')
TLS_obj@meta.data=df_TLS
TLS_obj<-TLS_obj[,TLS_obj$cell_state != 'Unknown']
TLS_96hr <- TLS_obj[,TLS_obj$timepoint=='96h']
TLS_108hr <- TLS_obj[,TLS_obj$timepoint=='108h']
TLS_120hr<- TLS_obj[,TLS_obj$timepoint=='120h']
Idents(TLS_96hr)<-TLS_96hr$cell_state
Idents(TLS_108hr)<-TLS_108hr$cell_state
Idents(TLS_120hr)<-TLS_120hr$cell_state

#Human gastruloid
hG_0h_rna <- readRDS(file=paste0(processed.dir,'hGastruloids_0h_processed.RDS'))
hG_24h_rna<- readRDS(file=paste0(processed.dir,'hGastruloid_24h_merged_ccr_annotated.RDS'))
hG_48h_rna <- readRDS(paste0(processed.dir,'hGastruloid_48h_merged_ccr_annotated.RDS'))
hG_72h_rna <- readRDS(file=paste0(processed.dir,'hGastruloid_72h_merged_ccr_annotated.RDS'))
hG_96h_rna <- readRDS(file=paste0(processed.dir,'hGastruloid_96h_merged_ccr_annotated.RDS'))
hG_RA_96h_rna <- readRDS(file=paste0(processed.dir,'RA_hGastruloid_96h_merged_ccr_annotated.RDS'))
hG_RA_120h_rna = readRDS(paste0(processed.dir,'RA_hGastruloid_120h_merged_ccr_annotated.RDS'))

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
saveRDS(Matrix(data_mtx,sparse=TRUE), '~/nobackups/MOCA_mtx_noExE.RDS')

#Human RA gastruloid 96hr, 120hr (mito-depletion vs non mito-depletion)
data_list <- list(
                  hG_RA_120h_rna=pseudobulk_human(hG_RA_120h_rna,liftover=T), 
                  hG_RA_96h_rna=pseudobulk_human(hG_RA_96h_rna,liftover=T),
                  hG_96h_rna=pseudobulk_human(hG_96h_rna,liftover=T), 
                  hG_72h_rna=pseudobulk_human(hG_72h_rna,liftover=T),
                  hG_48h_rna=pseudobulk_human(hG_48h_rna,liftover=T),hG_24h_rna=pseudobulk_human(hG_24h_rna,liftover=T),
                  hG_0h_rna=pseudobulk_human(hG_0h_rna,liftover=T)
                    )

#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
#select_genes <-Reduce(union, gene_list)%>%unique
data_list_subset<-data_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes]
    #data_list_subset[[i]][is.na(data_list_subset[[i]])] = 0
}
data_mtx <- do.call(cbind,data_list_subset)
saveRDS(Matrix(data_mtx,sparse=TRUE), '~/nobackups/query_mtx_hG_noliftover.RDS')

#TLS mouse gastruloid
# #VDB mouse gastruloid
data_list <- list(
                  mG_120h_rna=pseudobulk_mouse(mG_120h_rna_meso),
                  TLS_96hr=pseudobulk_mouse(TLS_96_rna_meso),TLS_108hr=pseudobulk_mouse(TLS_108_rna_meso),
                  TLS_120hr=pseudobulk_mouse(TLS_120_rna_meso)
                    )
#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
#select_genes <-Reduce(union, gene_list)%>%unique
data_list_subset<-data_list
for(i in 1:length(data_list_subset)){
    data_list_subset[[i]] <- data_list_subset[[i]][select_genes]
    #data_list_subset[[i]][is.na(data_list_subset[[i]])] = 0
}
data_mtx <- do.call(cbind,data_list_subset)
saveRDS(Matrix(data_mtx,sparse=TRUE), '~/nobackups/query_mtx_mG_noliftover.RDS')