#This script performs PCA-based whole organism molecular staging
#Pseudobulk matrices are outputs from Pseudobulk_datasets.R and provided within the same folder
#Link to human embryo CS7 dataset (Tyser et al., 2021) and CS12-16 dataset (Xu et al., 2021) are provided in https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/

suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(pheatmap)
        })
        
set.seed(2021)

data_dir = './hGastruloid_website/public/' #https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/
df_biomart <- read.table(paste0(data_dir,'/ref/mart_export_human_mouse.txt'), header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)

#MOCA 
MOCA_mtx <- readRDS('./MOCA_mtx_noExE.RDS')
MOCA_mtx_meta <- readRDS('./MOCA_mtx_meta.RDS')
embryo_info <- readRDS('./MOCA_E85_embryo_info.rds')
embryo_info$embryo = gsub('E8.5_','E8.5b_',embryo_info$RT_group)
MOCA_cell <- data.frame(embryo=colnames(MOCA_mtx),day=str_split_fixed(colnames(MOCA_mtx),'_',2)[,1])
MOCA_cell$dataset='MOCA'
MOCA_cell <- left_join(MOCA_cell, embryo_info)%>%select(c('embryo','day','dataset','somite_number'))
row.names(MOCA_cell) <- MOCA_cell$embryo
df_mouse_human <- df_biomart[df_biomart$Mouse_ID %in% row.names(MOCA_mtx),]
MOCA_mtx <- MOCA_mtx[df_mouse_human$Mouse_ID,]
row.names(MOCA_mtx) <- df_mouse_human$link

#Pijuan
pijuan_mtx <- readRDS('./pijuan_mtx_noExE.RDS')
pijuan_cell <- data.frame(embryo=colnames(pijuan_mtx),day=str_split_fixed(colnames(pijuan_mtx),'_',2)[,1],
                          dataset='pijuan',somite_number=NA)
row.names(pijuan_cell) <- pijuan_cell$embryo
df_pijuan <- df_biomart[df_biomart$Mouse_ID %in% row.names(pijuan_mtx),]
pijuan_mtx_link <- pijuan_mtx[df_pijuan$Mouse_ID,]
row.names(pijuan_mtx_link) <- df_pijuan$link

#Gastruloids & somitoids
query_mtx_hG  <- readRDS('./query_mtx_hG_noliftover.RDS')
query_mtx_st  <- readRDS('./query_mtx_somitoid_noliftover.RDS')
query_mtx_mG <- readRDS('./query_mtx_mG_noliftover.RDS')
query_cell_hG <- data.frame(embryo=colnames(query_mtx_hG),
                           day=gsub('_rna','',colnames(query_mtx_hG)),
                           dataset='gastruloids',
                          somite_number=NA)
row.names(query_cell_hG)=query_cell_hG$embryo
query_cell_mG <- data.frame(embryo=colnames(query_mtx_mG),
                           day=gsub('_rna','',colnames(query_mtx_mG)),
                           dataset='gastruloids',
                          somite_number=NA)
row.names(query_cell_mG)=query_cell_mG$embryo
query_cell_st <- data.frame(embryo=colnames(query_mtx_st),
                           day=gsub('_rna','',colnames(query_mtx_st)),
                           dataset='somitoids',
                          somite_number=NA)

#HEOA human embryos
#Preprocess HEOA dataset as described in paper
HEOA_dir <- paste0(data_dir,'human_CS12-16_Atlas/')
HEOA_cell <- read_csv(paste0(HEOA_dir,'GSE157329_cell_annotate.csv'))
HEOA_cell <- HEOA_cell%>%data.frame
HEOA_gene <- read_csv(paste0(HEOA_dir,'GSE157329_gene_annotate.csv'))
HEOA_gene <- HEOA_gene%>%data.frame
HEOA_mtx <- readMM(paste0(HEOA_dir,'GSE157329_raw_counts.mtx'))
row.names(HEOA_cell) <- HEOA_cell$cell_id
row.names(HEOA_gene) <- HEOA_gene$gene_id
colnames(HEOA_mtx) <- row.names(HEOA_cell)
row.names(HEOA_mtx) <- row.names(HEOA_gene)
HEOA_cell$embryo_id <- str_split_fixed(HEOA_cell$cell_id,'_',2)[,1]
# prepare genes to be removed: hemoglobin genes, MT genes, sex-specific genes, cell cycle genes, and batch-effect genes.
list_dir = paste0(HEOA_dir,'list/')
cc_gene<- scan(paste0(list_dir,'total_cc.txt'),what=character(0)) # cell cycle genes, merged from: pubmed ID 25378319; pubmed ID 30452682
hb_gene<- scan(paste0(list_dir,'hb_gene.txt'),what=character(0)) # hemoglobin genes
mt_gene<-HEOA_gene$gene_id[grepl(pattern='^MT-', x= HEOA_gene$gene_short_name)] # MT genes
bt_gene<- unique(HEOA_gene[HEOA_gene$gene_short_name %in% scan(paste0(list_dir,'batch.txt'),what=character(0)),]$gene_id)  # batch-effect genes (including sex-specific genes), merged from: pubmed ID 30096314; pubmed ID 31835037
rm_gene<- unique(c(cc_gene, hb_gene, mt_gene, bt_gene))
hvg <- readLines(paste0(list_dir,'hvg.txt'))
pb_mtx <- readRDS(paste0(list_dir,'pseudobulk_heoa.RDS'))
pb_mtx_hvg <- pb_mtx[hvg,]
pb_cell <- HEOA_cell%>%group_by(embryo, stage)%>%summarise(n=n())%>%data.frame
row.names(pb_cell) <- pb_cell$embryo
pb_cell <- pb_cell %>% select(c('embryo','stage')) 
names(pb_cell)[2]='day'
pb_cell$dataset='HEOA'
pb_cell$somite_number=NA

##########marker
#CS7
human_CS7_mtx <- readRDS(paste0(data_dir,'human_CS7/human_CS7_raw_matrix.rds'))
human_CS7_meta <- readRDS(paste0(data_dir,'human_CS7/human_CS7_gastrula_meta.rds'))
human_CS7_mtx <- t(human_CS7_mtx)
human_CS7_rna <- CreateSeuratObject(human_CS7_mtx, project = "human_CS7", 
                   assay = "RNA",
                   min.cells = 10, 
                   min.features = 2000, 
                   meta.data = human_CS7_meta)
human_CS7_rna <- human_CS7_rna[,!human_CS7_rna$sub_cluster %in% c('Hypoblast','YS Endoderm','YS Mesoderm',
                                                                 'Blood Progenitors','Erythro-Myeloid Progenitors',
                                                                 'Erythroblasts','Myeloid Progenitors',
                                                                  'PGC','Non-Neural Ectoderm'
                                                                 ) ]
counts <- GetAssayData(human_CS7_rna, assay = "RNA",slot='count')
rm_gene_2 <- rownames(counts)[which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|
^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|
^IGFBP|^HSP|^FOS|^JUN"))]
select_genes <- setdiff(row.names(counts),union(df_biomart[df_biomart$Human_ID %in% rm_gene,]$Human_name,rm_gene_2))
counts <- counts[select_genes,]
human_CS7_rna_subset <- subset(human_CS7_rna, features = rownames(counts))
human_CS7_rna_subset <- NormalizeData(human_CS7_rna_subset)
human_CS7_rna_subset <- FindVariableFeatures(human_CS7_rna_subset)
human_CS7_rna_subset <- ScaleData(human_CS7_rna_subset)
CS7_hvg <- VariableFeatures(human_CS7_rna_subset)
pb_CS7_n <- rowSums(GetAssayData(human_CS7_rna_subset, assay = "RNA",slot='data'))
CS7_hvg <- VariableFeatures(human_CS7_rna_subset)
cs7_cell <- data.frame(embryo='emb_cs7',day='CS7',dataset='human_embryo_cs7',somite_number=NA)
row.names(cs7_cell) <- cs7_cell$embryo


#########################################
#######PCA Analysis
#########################################

#Subset HEOA dataset to CS12/13
pb_mtx_hvg_sub <- pb_mtx_hvg[,1:2]
row.names(pb_mtx_hvg_sub) <- row.names(pb_mtx_hvg)
#Convert gene IDs of HEOA to match gene names of CS7
df_heoa <- df_biomart[df_biomart$Human_ID %in% row.names(pb_mtx_hvg_sub),]
pb_mtx_hvg_sub <- pb_mtx_hvg_sub[df_heoa$Human_ID,]
row.names(pb_mtx_hvg_sub) <- df_heoa$Human_name

#Combine heoa and CS7 matrix
homologs <- intersect(names(pb_CS7_n),row.names(pb_mtx_hvg_sub))
length(homologs)
heoa_sub_cs7_mtx <- cbind(pb_mtx_hvg_sub[homologs,],pb_CS7_n[homologs])
colnames(heoa_sub_cs7_mtx)[3]='emb_cs7'
heoa_sub_cs7_cell <- rbind(pb_cell[1:2,],cs7_cell)
#Combine query_mtx (hG) and CS7-CS12
homologs <- Reduce(intersect,list(row.names(query_mtx_hG),row.names(heoa_sub_cs7_mtx)))
heoa_cs7_hG_mtx <- cbind(query_mtx_hG[homologs,],heoa_sub_cs7_mtx[homologs,])
heoa_cs7_hG_cell<-rbind(query_cell_hG,heoa_sub_cs7_cell)
row.names(heoa_cs7_hG_cell) <- heoa_cs7_hG_cell$embryo

#Gene sets
gene_set <- readRDS('./PCA_geneset.RDS')
heoa_cs7_hG_mtx_subset <- heoa_cs7_hG_mtx[row.names(heoa_cs7_hG_mtx) %in% gene_set,]

#Co-embed
obj_heoa_cs7_hG <- CreateSeuratObject(counts=heoa_cs7_hG_mtx_subset,meta.data=heoa_cs7_hG_cell)
dim(obj_heoa_cs7_hG)
obj_heoa_cs7_hG  <- NormalizeData(obj_heoa_cs7_hG ,verbose=F)
VariableFeatures(obj_heoa_cs7_hG)  <- row.names(obj_heoa_cs7_hG)
obj_heoa_cs7_hG <- ScaleData(obj_heoa_cs7_hG,verbose=F)
obj_heoa_cs7_hG  <- RunPCA(obj_heoa_cs7_hG,verbose=F,npcs=6)
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(obj_heoa_cs7_hG,group.by='day',pt.size=2.5,label=T,repel=T,label.size=4)+NoLegend()+ggtitle('Human sample all cells')

#Correlation
emb_heoa_cs7_hG <- Embeddings(obj_heoa_cs7_hG,reduction='pca')
cor <- cor(emb_heoa_cs7_hG[,'PC_2'],t(obj_heoa_cs7_hG[['RNA']]@scale.data),method='pearson')
names(cor)<-row.names(obj_heoa_cs7_hG)
cor_top100 <- cor[sort(abs(cor),decreasing=T)%>%head(100)%>%names%>%unique]
cor_top100_pos <- cor_top100[cor_top100>0]%>%names
cor_top100_neg <- cor_top100[cor_top100<0]%>%names
options(repr.plot.width=4, repr.plot.height=20)
mat <- GetAssayData(obj_heoa_cs7_hG,assay='RNA',slot='scale.data')
colnames(mat)=obj_heoa_cs7_hG$day
options(repr.plot.width=8, repr.plot.height=12)
annotation_df <- data.frame(row.names=names(cor_top100),annotation=ifelse(cor_top100>0,'negative','positive'))
ph=pheatmap(mat[names(cor_top100),],cluster_cols = TRUE, annotation_row = annotation_df)

###########################
#Projection of somitoids onto human gastruloid and embryo space
###########################
heoa_cs7_hG_mtx <- GetAssayData(obj_heoa_cs7_hG,assay='RNA',slot='count')
query_mtx_st_link <-query_mtx_st

#Combine query_mtx (MOCA,mG,hG) and CS7-CS12
homologs <- Reduce(intersect,list(row.names(query_mtx_st_link),row.names(heoa_cs7_hG_mtx)))
heoa_cs7_hG_st_mtx <- cbind(query_mtx_st_link[homologs,],heoa_cs7_hG_mtx[homologs,])
heoa_cs7_hG_st_cell<-rbind(query_cell_st,heoa_cs7_hG_cell)
row.names(heoa_cs7_hG_st_cell) <- heoa_cs7_hG_st_cell$embryo

irlba_heoa_cs7_hG <- Loadings(obj_heoa_cs7_hG,reduction='pca')
obj_heoa_cs7_hG_st <- CreateSeuratObject(counts=heoa_cs7_hG_st_mtx,meta.data=heoa_cs7_hG_st_cell)
obj_heoa_cs7_hG_st <- NormalizeData(obj_heoa_cs7_hG_st ,verbose=F)
obj_heoa_cs7_hG_st <- ScaleData(obj_heoa_cs7_hG_st ,verbose=F)
homologs <- intersect(row.names(obj_heoa_cs7_hG_st[['RNA']]@scale.data),row.names(irlba_heoa_cs7_hG))
VariableFeatures(obj_heoa_cs7_hG_st)=row.names(obj_heoa_cs7_hG_st)
heoa_cs7_hG_st_emb <- t(obj_heoa_cs7_hG_st[['RNA']]@scale.data[homologs,]) %*% irlba_heoa_cs7_hG[homologs,]
obj_heoa_cs7_hG_st[['pca']] <- CreateDimReducObject(embeddings = heoa_cs7_hG_st_emb , key = "PCA_", assay = DefaultAssay(obj_heoa_cs7_hG_st))
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(obj_heoa_cs7_hG_st , group.by='day',pt.size=3,label=T,repel=T,label.size=5)+NoLegend()

###########################
#Projection of mouse gastruloids, MOCA onto human gastruloid and embryo space
###########################
irlba_heoa_cs7_hG <- Loadings(obj_heoa_cs7_hG,reduction='pca')
df_irlba <- df_biomart[df_biomart$Human_name %in% row.names(irlba_heoa_cs7_hG),]
irlba_heoa_cs7_hG <- irlba_heoa_cs7_hG[df_irlba$Human_name,]
row.names(irlba_heoa_cs7_hG) <- df_irlba$link

df_human <- df_biomart[df_biomart$Human_name %in% row.names(heoa_cs7_hG_st_mtx),]
heoa_cs7_hG_st_mtx <- heoa_cs7_hG_st_mtx[df_human$Human_name,]
row.names(heoa_cs7_hG_st_mtx) <- df_human$link
#Combine query_mtx (MOCA,mG,hG) and CS7-CS12
MOCA_cell_8_11 <- MOCA_cell[MOCA_cell$day %in% c('E8.5b','E9.5','E10.5'),]
homologs <- Reduce(intersect,list(row.names(MOCA_mtx[,row.names(MOCA_cell_8_11)]),row.names(heoa_cs7_hG_st_mtx)))
moca_heoa_cs7_hG_mtx <- cbind(MOCA_mtx[homologs,row.names(MOCA_cell_8_11)],heoa_cs7_hG_st_mtx[homologs,])
moca_heoa_cs7_hG_cell<-rbind(MOCA_cell_8_11,heoa_cs7_hG_st_cell)
row.names(moca_heoa_cs7_hG_cell) <- moca_heoa_cs7_hG_cell$embryo
#Projection
obj_moca_heoa_cs7_hG <- CreateSeuratObject(counts=moca_heoa_cs7_hG_mtx,meta.data=moca_heoa_cs7_hG_cell)
obj_moca_heoa_cs7_hG <- NormalizeData(obj_moca_heoa_cs7_hG ,verbose=F)
obj_moca_heoa_cs7_hG <- ScaleData(obj_moca_heoa_cs7_hG ,verbose=F)
homologs <- intersect(row.names(obj_moca_heoa_cs7_hG [['RNA']]@scale.data),row.names(irlba_heoa_cs7_hG))
moca_heoa_cs7_hG_emb <- t(obj_moca_heoa_cs7_hG[['RNA']]@scale.data[homologs,]) %*% irlba_heoa_cs7_hG[homologs,]
obj_moca_heoa_cs7_hG[['pca']] <- CreateDimReducObject(embeddings = moca_heoa_cs7_hG_emb , key = "PCA_", assay = DefaultAssay(obj_moca_heoa_cs7_hG))
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(obj_moca_heoa_cs7_hG , group.by='day',pt.size=3,label=T,repel=T,label.size=5)+NoLegend()

###########################
#Projection of mouse gastruloids, MOCA and pijuan onto human gastruloid and embryo space
###########################
#Convert pijuan genes to link
pijuan_cell_7_8 <- pijuan_cell[pijuan_cell$day %in% c(
                                                      'E7','E7.25','E7.5','E7.75','E8','E8.25','E8.5'),]
df_pijuan <- df_biomart[df_biomart$Mouse_ID %in% row.names(pijuan_mtx),]
pijuan_mtx_link <- pijuan_mtx[df_pijuan$Mouse_ID,]
row.names(pijuan_mtx_link) <- df_pijuan$link
#Combine query_mtx (pijuan,MOCA,mG,hG) and CS7-CS12
homologs <- Reduce(intersect,list(row.names(pijuan_mtx_link),row.names(moca_heoa_cs7_hG_mtx)))
pijuan_moca_heoa_cs7_hG_mtx <- cbind(pijuan_mtx_link[homologs,row.names(pijuan_cell_7_8)],moca_heoa_cs7_hG_mtx[homologs,])
pijuan_moca_heoa_cs7_hG_cell<-rbind(pijuan_cell_7_8,moca_heoa_cs7_hG_cell)
row.names(pijuan_moca_heoa_cs7_hG_cell) <- pijuan_moca_heoa_cs7_hG_cell$embryo
#Projection
obj_pijuan_moca_heoa_cs7_hG <- CreateSeuratObject(counts=pijuan_moca_heoa_cs7_hG_mtx,meta.data=pijuan_moca_heoa_cs7_hG_cell)
obj_pijuan_moca_heoa_cs7_hG <- NormalizeData(obj_pijuan_moca_heoa_cs7_hG ,verbose=F)
obj_pijuan_moca_heoa_cs7_hG <- ScaleData(obj_pijuan_moca_heoa_cs7_hG ,verbose=F)
homologs <- intersect(row.names(obj_pijuan_moca_heoa_cs7_hG [['RNA']]@scale.data),row.names(irlba_heoa_cs7_hG))
pijuan_moca_heoa_cs7_hG_emb <- t(obj_pijuan_moca_heoa_cs7_hG[['RNA']]@scale.data[homologs,]) %*% irlba_heoa_cs7_hG[homologs,]
obj_pijuan_moca_heoa_cs7_hG[['pca']] <- CreateDimReducObject(embeddings = pijuan_moca_heoa_cs7_hG_emb , key = "PCA_", assay = DefaultAssay(obj_pijuan_moca_heoa_cs7_hG))
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(obj_pijuan_moca_heoa_cs7_hG , group.by='day',pt.size=2.5,label=T,repel=T,label.size=4)+NoLegend()

###########################
#Projection of mouse gastruloids onto human gastruloid and embryo space
###########################
df_mG <- df_biomart[df_biomart$Mouse_name %in% row.names(query_mtx_mG),]
query_mtx_mG_link <- query_mtx_mG[df_mG$Mouse_name,]
row.names(query_mtx_mG_link) <- df_mG$link
#Combine query_mtx (hG) and CS7-CS12
homologs <- Reduce(intersect,list(row.names(query_mtx_mG_link),row.names(pijuan_moca_heoa_cs7_hG_mtx)))
heoa_cs7_allG_mtx <- cbind(query_mtx_mG_link[homologs,],pijuan_moca_heoa_cs7_hG_mtx[homologs,])
heoa_cs7_allG_cell<-rbind(query_cell_mG,pijuan_moca_heoa_cs7_hG_cell)
row.names(heoa_cs7_allG_cell) <- heoa_cs7_allG_cell$embryo
#Projection
obj_heoa_cs7_allG <- CreateSeuratObject(counts=heoa_cs7_allG_mtx,meta.data=heoa_cs7_allG_cell)
obj_heoa_cs7_allG <- NormalizeData(obj_heoa_cs7_allG ,verbose=F)
obj_heoa_cs7_allG <- ScaleData(obj_heoa_cs7_allG ,verbose=F)
homologs <- intersect(row.names(obj_heoa_cs7_allG [['RNA']]@scale.data),row.names(irlba_heoa_cs7_hG))
heoa_cs7_allG_emb <- t(obj_heoa_cs7_allG[['RNA']]@scale.data[homologs,]) %*% irlba_heoa_cs7_hG[homologs,]
obj_heoa_cs7_allG [['pca']] <- CreateDimReducObject(embeddings = heoa_cs7_allG_emb , key = "PCA_", assay = DefaultAssay(obj_heoa_cs7_hG))
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(obj_heoa_cs7_allG  , group.by='day',pt.size=2.5,label=T,repel=T,label.size=3)+NoLegend()

###########################
#PC2 alignment
###########################

obj_heoa_cs7_allG$Source=ifelse(obj_heoa_cs7_allG$day %in% c('CS7','CS12','CS13',
                                                'E6.5','E6.75', 'E7','E7.25','E7.5','E7.75','E8','E8.25',
                                             'E8.5','E8.5b','E9.5','E10.5','E11.5','E12.5','E13.5'), 'Embryo','Gastruloids')
df_pca2 <- obj_heoa_cs7_allG@meta.data
df_pca2$PCA_2 <- data.frame(Embeddings(obj_heoa_cs7_allG,reduction='pca'))$PCA_2
df_pca2$day <- factor(df_pca2$day, level=c(  'E7','E7.25','E7.5','E7.75','E8','E8.25',
                                             'E8.5',
                                            'E8.5b','E9.5','E10.5',
                                            'CS7','CS12','CS13',
                                            'mG_120h',
                                            'TLS_96hr','TLS_108hr',
                                            'TLS_120hr',
                                            'hG_0h','hG_24h','hG_48h','hG_72h','hG_96h',
                                            'Somitoid_168h_CHIR8_10',
                                            'hG_RA_96h',
                                             'hG_RA_120h'
                                             ))
df_mean = df_pca2%>%group_by(day)%>%summarise(mean_PC2=mean(-PCA_2))
df_mean = df_mean%>%head(10)
df_pca2 <- df_pca2[order(df_pca2$day),]
pdf('./staging.pdf',width=10,height=6)
ggplot(subset(df_pca2,!is.na(day)),aes(x=-PCA_2,y=day,color=day,shape=Source))+geom_point(size=4)+theme_classic()+
theme(text = element_text(size = 15))+ xlab('PC2 (Developmental time)')+
geom_vline(xintercept = df_mean$mean_PC2,linetype="dotted")
dev.off()
