suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(nnls)
library(tidyr)
    library(pheatmap)
        })
options(future.globals.maxSize= 23000 * 1024^2)
set.seed(1234)

processed.dir <- "/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/hGastruloids/"
fig.dir <- "~/shendure_lab/hGastruloid_figures/"

df_biomart <- read.table("/net/shendure/vol1/home/weiy666/shendure_lab/sci-fate/data/20201211_sci-plex-fate/mart_export_human_mouse.txt", header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)

correlation_analysis_nnls_spec_gene <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  cell_names = colnames(MM1)
  gene_median_expr = apply(MM1, 1, median)
  gene_max_expr = apply(MM1, 1, max)
  correlation_matrix = lapply(cell_names, function(x) {
    gene_other_max_expr = apply(MM1[, colnames(MM1) != x],1,  max)
    target_expr = MM1[, colnames(MM1) == x]
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_median_expr + 1), gene_id = row.names(MM1))
    gene_markers = (df_tmp %>% filter(ratio > fold.change) %>%  arrange(desc(ratio)) %>% head(top_gene_num))$gene_id
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_other_max_expr + 1), gene_id = row.names(MM1))
    top_markers = (df_tmp %>% arrange(desc(ratio)) %>% head(spec_gene_num))$gene_id
    selected_row = (row.names(MM1) %in% c(as.character(top_markers), as.character(gene_markers)))
    MM1_filtered = MM1[selected_row, ]
    MM2_filtered = MM2[selected_row, ]
    nnls_result = nnls::nnls(MM2_filtered, MM1_filtered[, x])
    df_coef = data.frame("cell_name_MM2" = colnames(MM2), "beta" = coef(nnls_result), "cell_name_MM1" = x)
    return(df_coef)
  })
  result = do.call(rbind, correlation_matrix)
  result = result %>% spread(cell_name_MM2, beta)
  result_2 = result %>% dplyr::select(-cell_name_MM1)
  rownames(result_2) = result$cell_name_MM1
  return(result_2)
}

correlation_analysis_bidirection <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  com_1 = correlation_analysis_nnls_spec_gene(MM1, MM2, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_1)
  result$source = rownames(result)
  result_1 = result %>% gather(key = target, value = beta_1, 1:ncol(MM2))
  com_2 = correlation_analysis_nnls_spec_gene(MM2, MM1, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_2)
  result$target = rownames(result)
  result_2 = result %>% gather(key = source, value = beta_2, 1:ncol(MM1))
  result = inner_join(result_1, result_2)
  return(result)
}
correlation_analysis_cell_type_selection <- function(tmp_result, marker_num = 1) {
  tmp_result$beta = (2 * (tmp_result$beta_1 + 0.01) * (tmp_result$beta_2 + 0.01))
  tmp = tmp_result %>% dplyr::select(-beta_1) %>% dplyr::select(-beta_2) %>% spread(key = source, value = beta)
  whole_matrix = as.matrix(tmp %>% dplyr::select(-target))
  rownames(whole_matrix) = tmp$target
  return(list(whole_matrix))
}

##Pseudobulk functions
pseudobulk_moca <- function(dataset, project_name){
    obj_1=dataset
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
    
    MOCA_mtx <- obj_1[['RNA']]@counts
    res_moca <- NULL
    for(Anno in obj_1$Anno%>%unique()){
    Anno_cell <- obj_1[,obj_1$Anno == Anno]%>%colnames()
    aggr_expr <- rowSums(MOCA_mtx[,Anno_cell,drop=FALSE])
    res_moca <- cbind(res_moca,aggr_expr)
    }
    colnames(res_moca) <- obj_1$Anno%>%unique()
    res_moca = t(t(res_moca) / colSums(res_moca) * 100000)
    res_moca = log(res_moca + 1)
    return(res_moca)
    }
pseudobulk_hG <- function(dataset,timepoint){
    obj_2=dataset
    obj_2$gastruloid.cluster=Idents(obj_2)
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
    obj_2$timepoint = timepoint
    
    
    hG_mtx <- obj_2[['RNA']]@counts
    res_hG <- NULL
    for(gastruloid.cluster in obj_2$gastruloid.cluster%>%unique){
        cluster_cell <- obj_2[,obj_2$gastruloid.cluster == gastruloid.cluster]%>%colnames()
        aggr_expr <- rowSums(hG_mtx[,cluster_cell])
        res_hG<- cbind(res_hG,aggr_expr)
    }
    colnames(res_hG) <- obj_2$gastruloid.cluster%>%unique
    res_hG = t(t(res_hG) / colSums(res_hG) * 100000)
    res_hG = log(res_hG + 1)
    return(res_hG)
    }
pseudobulk_mG <- function(dataset,timepoint){
    obj_2=dataset
    obj_2$gastruloid.cluster=obj_2$VDB_celltype
    mg_cells <- obj_2@meta.data
    mg_gene <- obj_2[['RNA']][[]]
    mg_mtx <- GetAssayData(object = obj_2, slot = "counts")
    df_mouse_human = df_biomart[df_biomart$Mouse_name %in% rownames(mg_mtx),]
    mg_mtx_sub = mg_mtx[df_mouse_human$Mouse_name,]
    rownames(mg_mtx_sub) = df_mouse_human$link
    obj_2 = CreateSeuratObject(counts = mg_mtx_sub,
    project = "hGastruloids",
    assay = "RNA",
    meta.data =  mg_cells)
    obj_2$timepoint = timepoint
    
    
    mG_mtx <- obj_2[['RNA']]@counts
    res_mG <- NULL
    for(gastruloid.cluster in obj_2$gastruloid.cluster%>%unique){
        print(gastruloid.cluster)
        cluster_cell <- obj_2[,obj_2$gastruloid.cluster == gastruloid.cluster]%>%colnames()
        aggr_expr <- rowSums(mG_mtx[,cluster_cell])
        res_mG<- cbind(res_mG,aggr_expr)
    }
    colnames(res_mG) <- obj_2$gastruloid.cluster%>%unique
    res_mG = t(t(res_mG) / colSums(res_mG) * 100000)
    res_mG = log(res_mG + 1)
    return(res_mG)
    }

#Datasets
#E8.5b (TOME)
MOCA_e8_5b <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E8.5b_exp.rds')
MOCA_e8_5b.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E8.5b_anno.rds")
MOCA_e8_5b@meta.data$Anno <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$Anno
MOCA_e8_5b@meta.data$celltype <- MOCA_e8_5b.anno[colnames(MOCA_e8_5b),]$celltype
#E9.5 (MOCA)
MOCA_e9_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E9.5_exp.rds')
MOCA_e9_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E9.5_anno.rds")
MOCA_e9_5@meta.data$Anno <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$Anno
MOCA_e9_5@meta.data$celltype <- MOCA_e9_5.anno[colnames(MOCA_e9_5),]$celltype
#E10.5 (MOCA)
MOCA_e10_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E10.5_exp.rds')
MOCA_e10_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E10.5_anno.rds")
MOCA_e10_5@meta.data$Anno <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$Anno
MOCA_e10_5@meta.data$celltype <- MOCA_e10_5.anno[colnames(MOCA_e10_5),]$celltype
#Mouse gastruloid (Van de brink et al., 2020)
mGastruloids_120h_rna <- readRDS('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/mGastruloids/mGastruloids_120h_rna.RDS')
#Human gastruloids
hGastruloids_RA_120h_rna_merged <- readRDS(file=paste0(processed.dir,'hGastruloid_120h_merged_ccr_rmd.RDS'))
Idents(hGastruloids_RA_120h_rna_merged) = hGastruloids_RA_120h_rna_merged$cell_type
hGastruloids_RA_96h_rna = readRDS('~/nobackups/RA_hGastruloid_96h_merged_ccr_annotated.RDS')
hGastruloids_96h_rna = readRDS(file=paste0(processed.dir,'hGastruloid_96h_merged_ccr_annotated.RDS'))
hGastruloids_96h_rna_new <- readRDS(paste0('~/nobackups/hGastruloid_96h_merged_ccr_annotated.RDS'))

##Psuedobulk MOCA
res_moca_e85 <- pseudobulk_moca(MOCA_e8_5b,'MOCA_8.5b')
res_moca_e95 <- pseudobulk_moca(MOCA_e9_5,'MOCA_9.5')
res_moca_e105 <- pseudobulk_moca(MOCA_e10_5,'MOCA_10.5')
res_moca_e115 <- pseudobulk_moca(MOCA_e11_5,'MOCA_11.5')

select_cell_type <- c('Anterior floor plate','Definitive endoderm','Endothelium','First heart field',
                      'Intermediate mesoderm','Neural crest','Neuromesodermal progenitors','Notochord',
                      'Paraxial mesoderm A','Paraxial mesoderm B','Posterior floor plate',
                      'Primordial germ cells',
                      'Second heart field','Spinal cord','Splanchnic mesoderm',
                      'Gut','Gut and lung epithelium',
                      'Hepatocytes','Neural crest (PNS glia)','Renal epithelium','Roof plate',
                      'Skeletal muscle progenitors', 'Spinal cord (ventral)','Spinal cord (dorsal)',
                      'Cardiomyocytes','Myocytes','Neuron progenitor cells')
res_moca_e85_filtered <- res_moca_e85[,colnames(res_moca_e85)[gsub('E8.5b:','',colnames(res_moca_e85)) %in% select_cell_type ]]
res_moca_e95_filtered <- res_moca_e95[,colnames(res_moca_e95)[gsub('E9.5:','',colnames(res_moca_e95)) %in% select_cell_type ]]
res_moca_e105_filtered <- res_moca_e105[,colnames(res_moca_e105)[gsub('E10.5:','',colnames(res_moca_e105)) %in% select_cell_type ]]

##Pseudobulk mouse gastruloid
res_mG <- pseudobulk_mG(mGastruloids_120h_rna,'VDB_120h')

##Pseudobulk human gastruloid
res_RA_hG_120h<- pseudobulk_hG(hGastruloids_RA_120h_rna_merged,'RA_120h')
res_RA_hG_96h <- pseudobulk_hG(hGastruloids_RA_96h_rna,'RA_96h')
res_hG_96h<- pseudobulk_hG(hGastruloids_96h_rna,'ori_96h')
res_hG_96h_new <- pseudobulk_hG(hGastruloids_96h_rna_new, 'ori_96h_new')
res_RA_LDN_ra <- pseudobulk_hG(hG_RA_LDN_ra,'RA_LDN')

#RA gastruloids 120hr vs MOCA
homologs=Reduce(intersect, list(rownames(res_RA_hG_120h),rownames(res_moca_e85),rownames(res_moca_e95),rownames(res_moca_e105))) 
res_RA_hG_120h_filtered=res_RA_hG_120h[homologs,]
res_moca_e85_filtered=res_moca_e85_filtered[homologs,]
res_moca_e95_filtered=res_moca_e95_filtered[homologs,]
res_moca_e105_filtered=res_moca_e105_filtered[homologs,]
res_moca <-cbind(res_moca_e85_filtered,cbind(res_moca_e95_filtered,res_moca_e105_filtered))
nnls_result <- correlation_analysis_bidirection(res_RA_hG_120h_filtered,res_moca,fold.change = 2.5, top_gene_num = 200, spec_gene_num = 200)
nnls_output_all <- correlation_analysis_cell_type_selection(nnls_result)
dat <- nnls_output_all[[1]]
df_tmp = data.frame(Anno=rownames(dat), day=str_split_fixed(rownames(dat),':',2)[,1])
df_tmp$day = ordered(df_tmp$day,levels=c('E8.5b','E9.5','E10.5'))
df_tmp = df_tmp%>%arrange(day)
dat = dat[df_tmp$Anno,c('Endothelium','Somite','Differentiating somite','Neural tube',
                       'Intermediate mesoderm','Neuron progenitor cells',
                       'Cardiomyocyte','Endoderm','Myocyte','Neural crest','Renal epithelium',
                       'Cardiac mesoderm')]
options(repr.plot.width=30, repr.plot.height=10)
pheatmap(t(dat), color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 18,
        cluster_rows=F, cluster_cols=F,scale='row')

#mouse gastruloids 120hr vs MOCA
homologs=Reduce(intersect, list(rownames(rownames(res_mG)),rownames(res_moca_e85),rownames(res_moca_e95),rownames(res_moca_e105))) 
res_mG_filtered=res_mG[homologs,colnames(res_mG)!='Allantois']
res_moca_e85_filtered=res_moca_e85_filtered[homologs,]
res_moca_e95_filtered=res_moca_e95_filtered[homologs,]
res_moca_e105_filtered=res_moca_e105_filtered[homologs,]
res_moca <-cbind(res_moca_e85_filtered,cbind(res_moca_e95_filtered,res_moca_e105_filtered))
nnls_result <- correlation_analysis_bidirection(res_mG_filtered,res_moca,fold.change = 2.5, top_gene_num = 200, spec_gene_num = 200)
nnls_output_all <- correlation_analysis_cell_type_selection(nnls_result)
dat <- nnls_output_all[[1]]

df_tmp = data.frame(Anno=rownames(dat), day=str_split_fixed(rownames(dat),':',2)[,1])
df_tmp$day = ordered(df_tmp$day,levels=c('E8.5b','E9.5','E10.5'#,'E11.5'#,
                                         #'E12.5','E13.5'
                                        ))
df_tmp = df_tmp%>%arrange(day)
dat = dat[df_tmp$Anno,c('Endothelium','Endoderm','NMP','Paraxial MD','Differentiation front','Somite','PSM',
                       'PGC-like','Cardiac','Neural tube','Differentiated somite')]

options(repr.plot.width=30, repr.plot.height=10)
pheatmap(t(dat), color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 18,
        cluster_rows=F, cluster_cols=F,scale='row')