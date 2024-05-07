biomart <- read.table("./mart_export_human_mouse.txt", header=T)
biomart$link = paste0(biomart$Mouse_ID,"-",biomart$Human_ID)

pseudobulk_mouse_by_embryo <- function(obj,gene_id,liftover=T){
    obj_mtx <- obj[['RNA']]@counts
    if(liftover){
    df_mouse_human = biomart[biomart[[gene_id]] %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[df_mouse_human[[gene_id]],]
    rownames(obj_mtx_sub) = df_mouse_human$link
    }else{
    obj_mtx_sub <- obj_mtx
    }
    embryo_mtx = model.matrix(~embryo_id-1, data = obj@meta.data%>%select(embryo_id))
    colnames(embryo_mtx) = gsub('embryo_id','',colnames(embryo_mtx))
    obj_sum = obj_mtx_sub %*% embryo_mtx
    
    return(obj_sum)
}

pseudobulk_human_by_embryo <- function(obj,gene_id,liftover=T){
    obj_mtx <- obj[['RNA']]@counts
    if(liftover){
    df_human_mouse = biomart[biomart[[gene_id]] %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[df_human_mouse[[gene_id]],]
    rownames(obj_mtx_sub) = df_human_mouse$link
    }else{
    obj_mtx_sub <- obj_mtx
    }
    
    embryo_mtx = model.matrix(~embryo_id-1, data = obj@meta.data%>%select(embryo_id))
    colnames(embryo_mtx) = gsub('embryo_id','',colnames(embryo_mtx))
    obj_sum = obj_mtx_sub %*% embryo_mtx
    
    return(obj_sum)
}

pseudobulk_mouse <- function(obj,liftover=T,gene_id){
    obj_mtx <- obj[['RNA']]@counts
    if(liftover){
    df_mouse_human = biomart[biomart[[gene_id]] %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[as.vector(df_mouse_human[[gene_id]]),]
    rownames(obj_mtx_sub) = df_mouse_human$link
    }else{
        obj_mtx_sub <- obj_mtx
    }
    
    obj_sum <- rowSums(obj_mtx_sub)
    
    return(obj_sum)
}

pseudobulk_human <- function(obj,liftover=T,gene_id){
    obj_mtx <- obj[['RNA']]@counts
    if(liftover){
    df_human_mouse = biomart[biomart[[gene_id]] %in% rownames(obj_mtx),]
    obj_mtx_sub = obj_mtx[df_human_mouse[[gene_id]],]
    rownames(obj_mtx_sub) = df_human_mouse$link
    }else{
        obj_mtx_sub <- obj_mtx
    }
    obj_sum <- rowSums(obj_mtx_sub)
    
    
    return(obj_sum)
}
