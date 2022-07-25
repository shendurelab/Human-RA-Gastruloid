#This script integrate human embryo CS7 dataset (Tyser et al., 2021) and human gastruloid at 24hr
#Link to human embryo CS7 dataset is provided in https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/
#Human gastrloid datasets from this paper are provided at National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) under accession numbers GSE208369

suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(reticulate)
library(harmony)
library(data.table)
        })

geo_dir = './geo/' #GSE208369
data_dir = './hGastruloid_website/public/' #https://shendure-web.gs.washington.edu/content/members/hGastruloid_website/public/

hGastruloids_24h_rna <- readRDS(paste0(geo_dir,'hGas_24h.RDS'))

human_CS7_mtx <- readRDS(paste0(data_dir,'human_CS7/human_CS7_raw_matrix.rds'))
human_CS7_meta <- readRDS(paste0(data_dir,'human_CS7/human_CS7_gastrula_meta.rds'))
human_CS7_rna <- CreateSeuratObject(t(human_CS7_mtx), project = "human_CS7", 
                   assay = "RNA",
                   min.cells = 10, 
                   min.features = 0, 
                   meta.data = human_CS7_meta)

counts <- GetAssayData(human_CS7_rna, assay = "RNA")
counts <- counts[-(which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|^IGFBP|^HSP|^FOS|^JUN"))), ]
human_CS7_rna <- subset(human_CS7_rna, features = rownames(counts))

human_CS7_rna <- NormalizeData(human_CS7_rna, verbose=F)
human_CS7_rna <- human_CS7_rna%>%FindVariableFeatures(selection.method = "vst",nfeatures=4000)%>%ScaleData(assay="RNA") 
human_CS7_rna <- RunPCA(human_CS7_rna,verbose=F)
human_CS7_rna <- RunUMAP(human_CS7_rna,,verbose=F,dims=1:30)

#Downsample 24hr gastruloid cluster to 58 cells/cluster to similar size of Human CS7 dataset
hGastruloids_24h_rna$cell.names = colnames(hGastruloids_24h_rna)
hGastruloids_24h_rna$Anno=Idents(hGastruloids_24h_rna)
df_tmp <- hGastruloids_24h_rna@meta.data%>%data.frame%>%group_by(seurat_clusters)%>%sample_n(size=58, replace=F)   
select_cells <- df_tmp$cell.names
hGastruloids_24h_rna_subset <- hGastruloids_24h_rna[,select_cells]

obj = merge(human_CS7_rna,hGastruloids_24h_rna_subset, add.cell.ids=c('CS7','24h'))
obj = NormalizeData(obj, verbose=F)
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
obj <- CellCycleScoring(obj, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = TRUE)

obj <- obj %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
obj <- obj %>% ScaleData(vars.to.regress=c('S.Score','G2M.Score'),assay="RNA") 
DefaultAssay(obj) <- "RNA"
obj <- obj%>% RunPCA(verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "orig.ident", assay.use="RNA",theta=3)
obj <- RunUMAP(obj, reduction = "harmony", min.dist = 0.3, dims = 1:30)

#Visualization
p1 <- DimPlot(obj, shuffle=TRUE, label=F, repel=T,label.size=5,pt.size=2,cols=,group.by='sub_cluster') + 
theme(text=element_text(size=30))+ labs(title='Human CS7 embryo')
colors <- ggplot_build(p1)$data[[1]]$colour
my_color = ggplot_build(p1)$data[[1]]$colour%>%unique
my_color = my_color[my_color!='grey50']
options(repr.plot.width=15, repr.plot.height=8)
DimPlot(obj, shuffle=TRUE, label=F, repel=T,label.size=5,pt.size=1.5,
        cols=my_color,
        group.by='sub_cluster',
       ) + 
theme(text=element_text(size=30))+ labs(title='Human CS7 embryo')+
guides(colour = guide_legend(override.aes = list(size=6)))
options(repr.plot.width=15, repr.plot.height=8)
DimPlot(obj, shuffle=TRUE, label=F, repel=T,label.size=5,pt.size=1.5,
        cols=c(my_color[c(1,4,6,5,16,2,13,14)]),
        group.by='Anno',
       ) + 
theme(text=element_text(size=30))+ labs(title='Human gastruloids (24 hrs)')+
guides(colour = guide_legend(override.aes = list(size=6)))