#This code import and process scRNA-seq from 24,48,72 and 96 hours of orginal human gastruloids (Moris et al., 2021)

suppressPackageStartupMessages({
options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(reticulate)
library(harmony)
library(data.table)
        })
set.seed(1234)

#Read 10x scRNA-seq output into Seurat object and generate metrics (UMI counts, feature counts, percentage of mitochndrial reads and doublet scores) for each cell. 
RNA.qc <- function(timepoint,dir_path=""){
    if(dir_path==""){
    rna_dir_path=paste0("/net/shendure/vol10/projects/hGastruloid/count_10X_RNA_",
                        timepoint,
                        "_gastruloid/outs/filtered_feature_bc_matrix")
        }else{
        rna_dir_path=dir_path
    }
    matrix=paste0(rna_dir_path,"/matrix.mtx.gz")
    #Make seurat object
    rna_mtx=Read10X(data.dir=rna_dir_path)
    hGastruloids_rna=CreateSeuratObject(rna_mtx)
    
    #Attach doublet info (scrublet does not seem to work well in our datasets)
    source_python("/net/shendure/vol1/home/weiy666/shendure_lab/sci-fate/bin/scrub_doublets.py")
    ds = predict_doublets(matrix)
    names(ds)=rownames(hGastruloids_rna[[]])
    hGastruloids_rna[["doublet.score"]]=ds
    
    #Calculate mitochondrial reads percentage
    hGastruloids_rna[["percent.mt"]] <- PercentageFeatureSet(hGastruloids_rna, pattern = "^MT-")
    #Add timepoint information
    hGastruloids_rna$orig.ident=timepoint
    #QC plot
    options(repr.plot.width=30, repr.plot.height=10)
    print(VlnPlot(hGastruloids_rna,features=c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3))
    return(hGastruloids_rna)
    }

#Filter cells based on the metrics from RNA.qc. Perform normalization, scaling with regression of cell cycle score and percentage of mitochondria reads, dimensional reduction, clustering and visualization. 
RNA.analysis <- function(obj, timepoint, min_count=1000, max_count=20000, min_feature=300, max_feature=3500, mt=10, 
                         doublet=0.6, mt_correct=FALSE){
    obj <- NormalizeData(obj, verbose = FALSE)
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    obj <- CellCycleScoring(obj, 
                                             s.features = s.genes, 
                                             g2m.features = g2m.genes, 
                                             set.ident = TRUE)
    obj=subset(obj, subset= nCount_RNA >min_count & nFeature_RNA > min_feature & 
               nFeature_RNA < max_feature & 
               nCount_RNA < max_count &
               percent.mt < mt &
              doublet.score < doublet)
    print(dim(obj))
    if(mt_correct==TRUE){
        obj=SCTransform(obj, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), verbose=FALSE, assay = 'RNA', new.assay.name = 'SCT')
    }else{
    obj=SCTransform(obj, vars.to.regress = c("S.Score", "G2M.Score"), verbose=FALSE, assay = 'RNA', new.assay.name = 'SCT')
        }
    DefaultAssay(obj) <- "SCT"
    obj=RunPCA(obj,npcs=30)
    obj=RunUMAP(obj, reduction = "pca", dims = 1:30,seed.use = 42)
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.8)
    
##########visualization#########
    options(repr.plot.width=30, repr.plot.height=10)
    p1=Seurat::DimPlot(obj, reduction = "umap",  label = TRUE, label.size = 8, repel = TRUE) + NoLegend()+ ggtitle("Cluster")
    p2=Seurat::FeaturePlot(obj, reduction = "umap", feature = "nCount_RNA") + ggtitle("nCount_RNA")
    p3=Seurat::DimPlot(obj, reduction = "umap", group.by = "Phase",label = FALSE, 
        label.size = 3, repel = TRUE) + ggtitle("cell.cycle")
    print(p1+p2+p3)
    
    return(obj)
}

#Gastruloid samples
#24hr human gastruloids
hGastruloids_24h_rna <- RNA.qc("24h", dir_path="/net/shendure/vol10/projects/hGastruloid/count_10X_RNA_24h_gastruloid_2ndBatch/outs/filtered_feature_bc_matrix/")
hGastruloids_24h_rna <- RNA.analysis(hGastruloids_24h_rna,min_count=4000, max_count=30000, 
                      min_feature=1500, max_feature=6000, mt=10, doublet=0.7, mt_correct=T)
#Remove clusters derived from technical noise
hGastruloids_24h_rna <- hGastruloids_24h_rna[,!hGastruloids_24h_rna$seurat_clusters %in% c(7,9)]
#Remove cell cycle, mitochondrial, non-coding RNA and ribosomal RNA genes that bias clustering results
counts <- GetAssayData(hGastruloids_24h_rna, assay = "RNA")
counts <- counts[-(which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|^IGFBP|^HSP|^FOS|^JUN"))), ]
hGastruloids_24h_rna <- subset(hGastruloids_24h_rna, features = rownames(counts))
#Normalization, scaling, dimensional reduction and clustering
hGastruloids_24h_rna=SCTransform(hGastruloids_24h_rna, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), verbose=FALSE, assay = 'RNA', new.assay.name = 'SCT')
DefaultAssay(hGastruloids_24h_rna) <- "SCT"
hGastruloids_24h_rna=RunPCA(hGastruloids_24h_rna,npcs=30)
hGastruloids_24h_rna=RunUMAP(hGastruloids_24h_rna, reduction = "pca", dims = 1:30,seed.use = 42)
hGastruloids_24h_rna <- FindNeighbors(hGastruloids_24h_rna, dims = 1:30)
hGastruloids_24h_rna <- FindClusters(hGastruloids_24h_rna, resolution = 0.8)

#48hr human gastruloids
hGastruloids_48h_rna <- RNA.qc("48h",dir_path="/net/shendure/vol10/projects/hGastruloid/count_10X_RNA_48h_gastruloid_2ndBatch/outs/filtered_feature_bc_matrix/")
hGastruloids_48h_rna <- RNA.analysis(hGastruloids_48h_rna,min_count=2000, max_count=30000, 
                       min_feature=1000, max_feature=6000, , doublet=0.69,mt=10)
#Remove cell cycle, mitochondrial, non-coding RNA and ribosomal RNA genes that bias clustering results
counts <- GetAssayData(hGastruloids_48h_rna, assay = "RNA")
counts <- counts[-(which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|^IGFBP|^HSP|^FOS|^JUN"))), ]
hGastruloids_48h_rna <- subset(hGastruloids_48h_rna, features = rownames(counts))
#Normalization, scaling, dimensional reduction and clustering
hGastruloids_48h_rna=SCTransform(hGastruloids_48h_rna, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), verbose=FALSE, assay = 'RNA', new.assay.name = 'SCT')
DefaultAssay(hGastruloids_48h_rna) <- "SCT"
hGastruloids_48h_rna=RunPCA(hGastruloids_48h_rna,npcs=30)
hGastruloids_48h_rna=RunUMAP(hGastruloids_48h_rna, reduction = "pca", dims = 1:30,seed.use = 42)
hGastruloids_48h_rna <- FindNeighbors(hGastruloids_48h_rna, dims = 1:30)
hGastruloids_48h_rna <- FindClusters(hGastruloids_48h_rna, resolution = 0.8)

#72hr human gastruloids
# Batch1
hGastruloids_72h_rna_b1 <- RNA.qc("72h")
hGastruloids_72h_rna_b1 <- RNA.analysis(hGastruloids_72h_rna_b1, doublet=0.36)
# Batch2
hGastruloids_72h_rna_b2 <- RNA.qc("72h",dir_path="/net/shendure/vol10/projects/hGastruloid/count_10X_RNA_72h_gastruloid_2ndBatch/outs/filtered_feature_bc_matrix/")
hGastruloids_72h_rna_b2 <- RNA.analysis(hGastruloids_72h_rna_b2,,min_count=1500, max_count=30000, 
                       min_feature=700, max_feature=4000, mt=10, doublet=0.7)
#Merge batches
hGastruloids_72h_rna_b1$batch='b1'
hGastruloids_72h_rna_b2$batch='b2'
hGastruloids_72h_rna <- merge(hGastruloids_72h_rna_b1,hGastruloids_72h_rna_b2, add.cell.ids=c('72h-1','72h-2'))
DefaultAssay(hGastruloids_72h_rna) <- "RNA"
#Remove cell cycle, mitochondrial, non-coding RNA and ribosomal RNA genes that bias clustering results
counts <- GetAssayData(hGastruloids_72h_rna, assay = "RNA")
counts <- counts[-(which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|^IGFBP|^HSP|^FOS|^JUN"))), ]
hGastruloids_72h_rna <- subset(hGastruloids_72h_rna, features = rownames(counts))
hGastruloids_72h_rna <- NormalizeData(hGastruloids_72h_rna, verbose=FALSE)
hGastruloids_72h_rna <- FindVariableFeatures(hGastruloids_72h_rna,assay="RNA") 
hGastruloids_72h_rna <- hGastruloids_72h_rna %>% ScaleData(vars.to.regress=c('S.Score','G2M.Score',"percent.mt"),assay="RNA") 
hGastruloids_72h_rna <- hGastruloids_72h_rna%>% RunPCA(verbose = FALSE)
#Batch correction with harmony
hGastruloids_72h_rna <- RunHarmony(hGastruloids_72h_rna, group.by.vars = "batch", assay.use="RNA")
hGastruloids_72h_rna <- RunUMAP(hGastruloids_72h_rna, reduction = "harmony", min.dist = 0.3, dims = 1:30)
hGastruloids_72h_rna <- FindNeighbors(hGastruloids_72h_rna, reduction = "harmony", dims = 1:30)
hGastruloids_72h_rna <- FindClusters(hGastruloids_72h_rna)
options(repr.plot.width=10, repr.plot.height=10)
DimPlot(hGastruloids_72h_rna, shuffle=TRUE, group.by = "seurat_clusters", label=T) + 
theme(text=element_text(size=30))

#96hr human gastruloids
#options(repr.plot.width=10, repr.plot.height=10)
hGastruloids_96h_rna <- RNA.qc("96h")
hGastruloids_96h_rna <- RNA.analysis(hGastruloids_96h_rna,min_count=2000,max_count=30000,doublet=0.55,mt_correct=TRUE)
#Remove cell cycle, mitochondrial, non-coding RNA and ribosomal RNA genes that bias clustering results
counts <- GetAssayData(hGastruloids_96h_rna[,hGastruloids_96h_rna$seurat_clusters!=12], assay = "RNA")
counts <- counts[-(which(str_detect(rownames(counts), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP|^IGFBP|^HSP|^FOS|^JUN"))), ]
hGastruloids_96h_rna <- subset(hGastruloids_96h_rna, features = rownames(counts))
#Normalization, scaling, dimensional reduction and clustering
hGastruloids_96h_rna=SCTransform(hGastruloids_96h_rna, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), verbose=FALSE, assay = 'RNA', new.assay.name = 'SCT')
DefaultAssay(hGastruloids_96h_rna) <- "SCT"
hGastruloids_96h_rna=RunPCA(hGastruloids_96h_rna,npcs=30)
hGastruloids_96h_rna=RunUMAP(hGastruloids_96h_rna, reduction = "pca", dims = 1:30,seed.use = 42)
hGastruloids_96h_rna <- FindNeighbors(hGastruloids_96h_rna, dims = 1:30)
hGastruloids_96h_rna <- FindClusters(hGastruloids_96h_rna, resolution = 0.8)

#Time-course
hGastruloids_24h_rna$orig.ident="24h"
hGastruloids_48h_rna$orig.ident="48h"
hGastruloids_72h_rna_b1$orig.ident="72h"
hGastruloids_72h_rna_b2$orig.ident="72h"
hGastruloids_96h_rna$orig.ident="96h"

hGastruloids_24h_rna$batch="b2"
hGastruloids_48h_rna$batch="b2"
hGastruloids_72h_rna_b1$batch="b1"
hGastruloids_72h_rna_b2$batch="b2"
hGastruloids_96h_rna$batch="b1"                                

human.gastr.tp=merge(x=hGastruloids_24h_rna,
                     y=c(
                      hGastruloids_48h_rna,
                      hGastruloids_72h_rna_b1,
                      hGastruloids_72h_rna_b2,
                      hGastruloids_96h_rna
                     ),
          add.cell.ids = c("24h",
                           "48h",
                           "72h-1",
                           "72h-2",
                           "96h"
                          ))
human.gastr.tp <- NormalizeData(human.gastr.tp, verbose = FALSE)
#Harmony
DefaultAssay(human.gastr.tp) <- "RNA"
human.gastr.tp <- human.gastr.tp %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)%>% 
            ScaleData(vars.to.regress=c('S.Score','G2M.Score',"percent.mt"))%>% 
            RunPCA(verbose = FALSE) 
options(repr.plot.width=10, repr.plot.height=10)
human.gastr.tp <- RunHarmony(human.gastr.tp, group.by.vars = "batch",theta=2)%>%
        RunUMAP(reduction = "harmony", min.dist = 0.3, dims = 1:30)%>%
        FindNeighbors(reduction = "harmony",dims=1:30)%>%
           FindClusters(resolution=0.3)
DimPlot(human.gastr.tp, shuffle=TRUE, group.by = "orig.ident") + 
theme(text=element_text(size=30))
