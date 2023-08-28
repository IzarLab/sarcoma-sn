library(numbat)
library(dplyr)
library(Seurat)
library(ggplot2)
library(glue)
library(data.table)
library(ggtree)
library(stringr)
library(tidygraph)
library(patchwork)

library(viridis)

library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
# Load Enrichr
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()


samps <- c("167", "322", "559", "708", "S408", "S410", "S914", "S956")

clone_db <- list()
samps <- c("167", "322", "559", "708", "S408", "S410", "S914", "S956")

for (pat in samps){
    
    mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3", 
              '5'="#ff9768", '6'='#ae1717', '7'='#0f04b5', '8'="#f87382",
              '9'='green','10'='yellow','11'='deeppink')
    nb = Numbat$new(out_dir = paste0('',pat, ''))
    # Sync-in sample Seurat object from AWS
    seu<- readRDS(paste0("../data_Sarcoma", pat, "GEX_genes_300_UMI_600_annotated_for_infercnv.rds"))
    seu$barcode_orig <- rownames(seu@meta.data)
    # Single-cell CNV calls
    cnv_calls<- nb$joint_post %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y)
    cnv_calls %>% group_by(cnv_state) %>% arrange(p_cnv,desc=F)
    # Clone info
    clones<-dim(table(nb$clone_post$clone_opt))
    clone_info<-nb$clone_post
    seu$cell<-seu$barcode_orig
    seu@meta.data<-left_join(seu@meta.data,clone_info,by='cell')
    rownames(seu@meta.data)<-seu$barcode_orig

    print(DimPlot(seu, group.by = 'clone_opt',
                shuffle = T, raster=T,cols = mypal[1:clones]))


    Idents(seu) <- 'clone_opt'
    markers.seu <- FindAllMarkers(seu)  
    for (clone in 1:clones){

        genes = markers.seu[markers.seu$cluster == clone & markers.seu$p_val_adj < 0.05 & 
                  markers.seu$avg_log2FC >0, ]$gene

        setEnrichrSite("Enrichr") # Human genes
        websiteLive <- TRUE
        dbs <- listEnrichrDbs()



        dbs <- c("MSigDB_Hallmark_2020")
        Sys.sleep(3)
        enriched <- enrichr(genes, dbs)
        if (!is.null(enriched[[1]])) {

            mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
            mut_enr_S410 = mut_enr_ch
            mut_enr <- mut_enr_ch
            h_mut_enr <- mut_enr_S410[1:25,]#[1:500,]
            mut_enr_S410$clone <- clone
            mut_enr_S410$samp <- pat

            clone_db <- append(clone_db, list(mut_enr_S410))
        }

    }

}

saveRDS(clone_db, "all_clone_pathways_to_graph.RDS")