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
df_biomart <- read.table("./mart_export_human_mouse.txt", header=T)
df_biomart$link = paste0(df_biomart$Mouse_ID,"-",df_biomart$Human_ID)
source('./pseudobulk_functions.R')

############Pijuan
Pijuan <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/data/pijuan.rds")
meta = Pijuan$pd
count = Pijuan$count
Pijuan = CreateSeuratObject(count,meta.data = meta)
Pijuan$embryo_id = paste0(Pijuan$day,'_',str_split_fixed(Pijuan$sample,'_',2)[,2])
select_celltype <- c('Epiblast','Primitive Streak','Nascent mesoderm','Anterior Primitive Streak',
                    'Mixed mesoderm','Haematoendothelial progenitors','Notochord','Blood progenitors',
                     'Somatic mesoderm','Splanchnic mesoderm','Caudal lateral epiblast','Caudal neurectoderm',
                     'Paraxial mesoderm A','Surface ectoderm','Gut',
                     'Paraxial mesoderm B','Paraxial mesoderm C','Primitive erythroid','Intermediate mesoderm',
                    'Endothelium','Cardiomyocytes','Spinal cord','Neural crest','Forebrain/Midbrain','Hindbrain',
                     'NMP','Somatic mesoderm A','Somatic mesoderm B','Urothelium','Placodal area','Epidermis')
Pijuan_sub <- Pijuan[,Pijuan$celltype %in% select_celltype]
pijuan_mtx = pseudobulk_mouse_by_embryo(Pijuan_sub,biomart=df_biomart)
saveRDS(pijuan_mtx ,'./pseudobulk_matrices/pijuan_mtx_noExE_liftover.RDS')

############MOCA
MOCA_e8_5b <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E8.5b_exp.rds')
MOCA_e9_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E9.5_exp.rds')
MOCA_e10_5 <- readRDS('/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/exp/obj_E10.5_exp.rds')
MOCA_e8_5b.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E8.5b_anno.rds")
MOCA_e9_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E9.5_anno.rds")
MOCA_e10_5.anno <- readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/tome/revision/anno/obj_E10.5_anno.rds")
#Pseudobulk by embryos
MOCA_e8_5b$celltype <- MOCA_e8_5b.anno$celltype
MOCA_e8_5b$embryo_id <- MOCA_e8_5b$RT_group
MOCA_e9_5$celltype <- MOCA_e9_5.anno$celltype
MOCA_e9_5$embryo_id <- paste0('E9.5_',MOCA_e9_5$embryo_id)
MOCA_e10_5$celltype <- MOCA_e10_5.anno$celltype
MOCA_e10_5$embryo_id <- paste0('E10.5_',MOCA_e10_5$embryo_id)
all_celltype <- Reduce(union,c(MOCA_e8_5b@meta.data$celltype,MOCA_e9_5@meta.data$celltype, 
                               MOCA_e10_5@meta.data$celltype))%>%unique
select_cells <- all_celltype[!all_celltype %in% c('Allantois','Extraembryonic visceral endoderm',
                      'Amniochorionic mesoderm B','Amniochorionic mesoderm A',
                     'Extraembryonic mesoderm','Mesenchymal stromal cells'
                     )]
MOCA_list <- list(
                  MOCA_e8_5b=pseudobulk_mouse_by_embryo(subset(MOCA_e8_5b, subset=celltype %in% select_cells),liftover=TRUE,gene_id='Mouse_ID'),
                  MOCA_e9_5=pseudobulk_mouse_by_embryo(subset(MOCA_e9_5, subset=celltype %in% select_cells),liftover=TRUE,gene_id='Mouse_ID'),
                  MOCA_e10_5=pseudobulk_mouse_by_embryo(subset(MOCA_e10_5, subset=celltype %in% select_cells),liftover=TRUE,gene_id='Mouse_ID')
                  )
gene_list <- NULL
for(i in 1:length(MOCA_list)){
    gene_list[[i]] <- rownames(MOCA_list[[i]]) 
}
select_genes <-Reduce(intersect, gene_list)
for(i in 1:length(MOCA_list)){
    MOCA_list[[i]] <- MOCA_list[[i]][select_genes,]
}
data_mtx <- do.call(cbind,MOCA_list)
saveRDS(data_mtx ,'./pseudobulk_matrices/moca_mtx_noExE_liftover.RDS')

#########HEOA
#Preprocess HEOA dataset as described in paper (1361 highly variable genes)
HEOA_dir <- '/net/shendure/vol1/home/weiy666/nobackups/journal_datasets/Human_CS12-16_Atlas/'
HEOA_cell <- read_csv(paste0(HEOA_dir,'GSE157329_cell_annotate.csv'))
HEOA_cell <- HEOA_cell%>%data.frame
HEOA_gene <- read_csv(paste0(HEOA_dir,'GSE157329_gene_annotate.csv'))
HEOA_gene <- HEOA_gene%>%data.frame
HEOA_mtx <- readMM(paste0(HEOA_dir,'GSE157329_raw_counts.mtx'))
row.names(HEOA_cell) <- HEOA_cell$cell_id
row.names(HEOA_gene) <- HEOA_gene$gene_id
colnames(HEOA_mtx) <- row.names(HEOA_cell)
row.names(HEOA_mtx) <- row.names(HEOA_gene)
HEOA_cell <- HEOA_cell%>%mutate(embryo_id = paste0(stage,'_',embryo))
obj_heoa = CreateSeuratObject(counts=HEOA_mtx,meta.data = HEOA_cell,min.features = 0, min.cells = 0)
pb_mtx = pseudobulk_human_by_embryo(obj_heoa,gene_id='Human_ID',liftover=F)
hvg <- readLines(paste0(HEOA_dir,'output/hvg.txt'))
pb_mtx_hvg <- pb_mtx[hvg,]
saveRDS(pb_mtx_hvg, './pseudobulk_matrices/pseudobulk_heoa_hvg.RDS')

# prepare genes to be removed: hemoglobin genes, MT genes, sex-specific genes, cell cycle genes, and batch-effect genes.
list_dir = paste0(HEOA_dir,'heoa/list/')
cc_gene<- scan(paste0(list_dir,'total_cc.txt'),what=character(0)) # cell cycle genes, merged from: pubmed ID 25378319; pubmed ID 30452682
hb_gene<- scan(paste0(list_dir,'hb_gene.txt'),what=character(0)) # hemoglobin genes
mt_gene<-HEOA_gene$gene_id[grepl(pattern='^MT-', x= HEOA_gene$gene_short_name)] # MT genes
bt_gene<- unique(HEOA_gene[HEOA_gene$gene_short_name %in% scan(paste0(list_dir,'batch.txt'),what=character(0)),]$gene_id)  # batch-effect genes (including sex-specific genes), merged from: pubmed ID 30096314; pubmed ID 31835037
rm_gene<- unique(c(cc_gene, hb_gene, mt_gene, bt_gene))
writeLines(rm_gene,paste0(HEOA_dir,'output/rm_genes.txt'))

#########cs7
#CS7 (excluding genes expressed in less than 10 cells and cells with <2000 UMIs)
human_CS7_mtx <- readRDS('/net/shendure/vol1/home/weiy666/nobackups/journal_datasets/Human_CS7/human_CS7_raw_matrix.rds')
human_CS7_meta <- readRDS('/net/shendure/vol1/home/weiy666/nobackups/journal_datasets/Human_CS7/human_CS7_gastrula_meta.rds')
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
^IGFBP|^HSP|^FOS|^JUN"))] #remove 2209 genes
select_genes <- setdiff(row.names(counts),union(df_biomart[df_biomart$Human_ID %in% rm_gene,]$Human_name,rm_gene_2))
counts <- counts[select_genes,]
human_CS7_rna_subset <- subset(human_CS7_rna, features = rownames(counts))
human_CS7_rna_subset <- NormalizeData(human_CS7_rna_subset)
pb_CS7_n <- rowSums(GetAssayData(human_CS7_rna_subset, assay = "RNA",slot='data'))
saveRDS(pb_CS7_n, './pseudobulk_matrices/pseudobulk_CS7_n.RDS')

#########human gastruloids
processed.dir = '/net/shendure/vol10/www/content/members/hGastruloid_website/public/Profiling/'
hG_0h_rna <- readRDS(file=paste0(processed.dir,'hGas_0h.RDS'))
hG_24h_rna <- readRDS(file=paste0(processed.dir,'hGas_24h.RDS'))
hG_48h_rna <- readRDS(paste0(processed.dir,'hGas_48h.RDS'))
hG_72h_rna <- readRDS(file=paste0(processed.dir,'hGas_72h.RDS'))
hG_96h_rna <- readRDS(file=paste0(processed.dir,'hGas_96h.RDS'))
hG_RA_96h_rna <- readRDS(file=paste0(processed.dir,'RA_hGas_96h.RDS'))
hG_RA_120h_rna <- readRDS(file=paste0(processed.dir,'RA_hGas_120h.RDS'))

data_list <- list(
                  hG_RA_120h_rna=pseudobulk_human(hG_RA_120h_rna,liftover=F,gene_id='Human_name'), 
                  hG_RA_96h_rna=pseudobulk_human(hG_RA_96h_rna,liftover=F,gene_id='Human_name'),
                  hG_96h_rna=pseudobulk_human(hG_96h_rna,liftover=F,gene_id='Human_name'), 
                  hG_72h_rna=pseudobulk_human(hG_72h_rna,liftover=F,gene_id='Human_name'),
                  hG_48h_rna=pseudobulk_human(hG_48h_rna,liftover=F,gene_id='Human_name'),
                  hG_24h_rna=pseudobulk_human(hG_24h_rna,liftover=F,gene_id='Human_name'),
                  hG_0h_rna=pseudobulk_human(hG_0h_rna,liftover=F,gene_id='Human_name'))

#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}
select_genes <-Reduce(intersect, gene_list)
for(i in 1:length(data_list)){
    data_list[[i]] <- data_list[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list)
saveRDS(data_mtx, './pseudobulk_matrices/query_human_matrix.RDS')

########Primate
## Adding primate samples to the PCA space
primate.dir <- '/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/raw_data/Zhai_Primate_CS8-11/MFE56636MTX'
expression_matrix <- Read10X(data.dir = primate.dir,gene.column=1)
primate_meta <- read.csv('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/raw_data/Zhai_Primate_CS8-11/MFE56636-meta.csv',
                       row.names=1) 

df_biomart <- read.table("./mart_export_human_primate.txt", header=T,sep = '\t')
df_biomart <- df_biomart%>%select(c('Human.gene.name','Gene.name'))%>%
            filter(!duplicated(Gene.name),!duplicated(Human.gene.name),Gene.name!='')
row.names(df_biomart) = df_biomart$Gene.name
df_primate_human = df_biomart[df_biomart$Gene.name %in% row.names(expression_matrix),]
expression_matrix_lf <- expression_matrix[df_primate_human$Gene.name,]
row.names(expression_matrix_lf) <- df_primate_human$Human.gene.name
primate = CreateSeuratObject(counts = expression_matrix_lf,meta.data=primate_meta,min.cells=10) #23283 genes

# #Select embryonic tissues
primate_subset <- primate[,primate$cell_type %in% c('APS','BP','Cardi.','Cardi.Meso','Caud.Meso','DE','EC',
                                                    'ECT','SE1','SE2',
                                                    'EPI','Ery1','Ery2','FB/MB/HB','Gut',
                                                   'Inter.Meso','LP.Meso','Mac','Mes','Nas.Meso','NC',
                                                   'NMP','Node','Para.Meso','PGC','Pharyg.Meso','PS',
                                                    'PSM','Rostr.Meso','SC')]
#Correct for batch effect related genes (mtRNA/lncRNA/cell cycle/proliferation)
counts <- GetAssayData(primate_subset, assay = "RNA",slot='count')
rm_gene_2 <- rownames(counts)[which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|
^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|
^IGFBP|^HSP|^FOS|^JUN"))]
select_genes <- setdiff(row.names(counts),union(df_biomart[df_biomart$Human_ID %in% rm_gene,]$Human_name,rm_gene_2))
counts <- counts[select_genes,] #21,873 genes
primate_subset <- subset(primate_subset, features = rownames(counts))

# primate_subset$embryo_id = primate_subset$sample
pb_primate <- pseudobulk_human_by_embryo(primate_subset,liftover=T,gene_id='Human_name')
saveRDS(pb_primate,'./pseudobulk_matrices/pb_primate_matrix_liftover.RDS')

######mouse Gastruloids
#VDB gastruloid (subset cells < 10% of the dataset and with >200 genes, 
# further filter out cells with 1000~40000 UMIs, and 700~8000 genes profiled by 10x)
mGastruloids_120h_rna <- readRDS('/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/mGastruloids/mGastruloids_120h_rna.RDS')
Idents(mGastruloids_120h_rna)=mGastruloids_120h_rna$VDB_celltype
#TLS gastruloid(subset cells < 10% of the dataset and with >200 genes,
#further filter out cells with 10000~40000 UMIs, >3000 genes, and mt% < 10)
TLS_obj <- readRDS(file='/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/mGastruloids/TLS_timecourse_rna.RDS')
TLS_96hr <- TLS_obj[,TLS_obj$timepoint=='96h']
TLS_108hr <- TLS_obj[,TLS_obj$timepoint=='108h']
TLS_120hr<- TLS_obj[,TLS_obj$timepoint=='120h']
Idents(TLS_96hr)<-TLS_96hr$cell_state
Idents(TLS_108hr)<-TLS_108hr$cell_state
Idents(TLS_120hr)<-TLS_120hr$cell_state
data_list <- list(
                  mG_120h_rna=pseudobulk_mouse(mGastruloids_120h_rna,liftover=T,gene_id='Mouse_name'),
                  TLS_96hr=pseudobulk_mouse(TLS_96hr,liftover=T,gene_id='Mouse_name'),
                  TLS_108hr=pseudobulk_mouse(TLS_108hr,liftover=T,gene_id='Mouse_name'),
                  TLS_120hr=pseudobulk_mouse(TLS_120hr,liftover=T,gene_id='Mouse_name'))

#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
for(i in 1:length(data_list)){
    data_list[[i]] <- data_list[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list)
saveRDS(data_mtx,'./pseudobulk_matrices/mouse_gastruloids_matrix_liftover.RDS')

#########human stembryos
data_dir <- '/net/shendure/vol8/projects/Wei_hGastruloid/nobackup/processed_data/other_models/'
#Somitoids
Somitoids_168h_CHIR8_10<-readRDS(paste0(data_dir,'human_somitoids/human_somitoids_b1.RDS'))
Somitoids_168h_CHIR5_7<- readRDS(paste0(data_dir,'human_somitoids/human_somitoids_MULTI-seq_b2.RDS'))
#Adding axioloid samples
axioloid.dir <- paste0(data_dir,'human_axioloids/')
axl_120h_MG_mtx <- Read10X(data.dir = paste0(axioloid.dir,'GSM5975296_120h_MG'),gene.column=1)
axl_120h_MG = CreateSeuratObject(counts = axl_120h_MG_mtx)
#Adding segmentoid samples
segmentoid.dir <- paste0(data_dir,'human_segmentoid/segmentoid_98h')
seg_98h_list <- Read10X(data.dir = segmentoid.dir,gene.column=1)
seg_98h = CreateSeuratObject(counts = seg_98h_list$`Gene Expression`)
#Adding EO (Libby et al)
neural_org.dir <- paste0(data_dir,'Libby_neural_organoid/')
neural_org_mtx <- Read10X(data.dir = neural_org.dir,gene.column=1)
neural_meta <- read.table(paste0(neural_org.dir,'/GSE155383_RPI4_final_calls.csv.gz'),header=T,row.names = 1,
                          sep=',')
colnames(neural_meta) <- 'MULTI-barcode'
neural_org  = CreateSeuratObject(counts = neural_org_mtx,
                                meta.data = neural_meta)
#Adding axial organoid (yaman et al)
yaman <- readRDS(paste0(data_dir,'yaman_axial_organoid/GSM6806916_seurat_sct_umap.RDS'))

#Adding EMLO
emlo <- readRDS(paste0(data_dir,'EMLO/EMLO_seurat_obj_subset_rm.RDS'))

#Griboud et al, Nature Biotech
pb_griboud = readRDS('~/nobackups/query_mtx_griboud.RDS')
pb_griboud_cell = data.frame(row.names=colnames(pb_griboud),embryo=colnames(pb_griboud),
                              day=colnames(pb_griboud),dataset='griboud',
                             somite_number=NA)

#pseudobulk
data_list <- list(
                  Somitoid_168h_CHIR5_7 = pseudobulk_human(Somitoids_168h_CHIR5_7,liftover=T,gene_id='Human_name'),
                  Somitoid_168h_CHIR8_10 = pseudobulk_human(Somitoids_168h_CHIR8_10,liftover=T,gene_id='Human_name'),
                  axl_120h_MG = pseudobulk_human(axl_120h_MG,liftover=T,gene_id='Human_ID'),
                  seg_98h = pseudobulk_human(seg_98h,liftover=T,gene_id='Human_ID'),
                  libby_EO = pseudobulk_human(neural_org,liftover=T,gene_id='Human_name'),
                  yaman_axial = pseudobulk_human(yaman,liftover=T,gene_id='Human_name'),
                  emlo = pseudobulk_human(emlo,liftover=T,gene_id='Human_name')
                    )



#Merge into single matrix
gene_list <- NULL
for(i in 1:length(data_list)){
    gene_list[[i]] <- names(data_list[[i]]) 
}

select_genes <-Reduce(intersect, gene_list)
for(i in 1:length(data_list)){
    data_list[[i]] <- data_list[[i]][select_genes]
}
data_mtx <- do.call(cbind,data_list)

saveRDS(data_mtx,'./pseudobulk_matrices/human_stembryos_matrix_liftover.RDS')