# Script to integrate scRNA-Seq datasets to correct for batch effects

### load libraries----
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratData)
library(patchwork)

### load data----
dirs <- list.dirs(path = 'data/', recursive = F, full.names = F)

for(x in dirs){
  name <-  gsub('_filtered_feature_bc_matrix','', x)
  
  cts <-  ReadMtx(mtx = paste0('data/',x, '/matrix.mtx.gz'), features = paste0('data/',x,'/features.tsv.gz'), cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
  }


### merging dataset---- 


merge_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor,
                             HB53_background, HB53_tumor),
      add.cell.ids = ls()[3:9],
      Project = 'HB')


### QC & Filtering----

view(merge_seurat@meta.data)


#Create a sample columns
merge_seurat$sample <- rownames(merge_seurat@meta.data)

# split sample comlumn (where there is "_")
merge_seurat@meta.data <- separate(merge_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcodes'),
         sep = '_')
# calculate MT%
merge_seurat$mtper <- PercentageFeatureSet(merge_seurat, pattern = '^MT-')

#explore QC

VlnPlot(merge_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(merge_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') # to compare and analyse the low quality data


# filtering

merge_seurat_fil <- subset(merge_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mtper < 10)


### perform standard workflow steps to figure out if we see any batch effects ----
merge_seurat_fil <- NormalizeData(object = merge_seurat_fil)
merge_seurat_fil <- FindVariableFeatures(object = merge_seurat_fil)
merge_seurat_fil <- ScaleData(object = merge_seurat_fil)
merge_seurat_fil <- RunPCA(object = merge_seurat_fil)
ElbowPlot(merge_seurat_fil)
merge_seurat_fil <- FindNeighbors(object = merge_seurat_fil, dims = 1:20)
merge_seurat_fil <- FindClusters(object = merge_seurat_fil)
merge_seurat_fil <- RunUMAP(object = merge_seurat_fil, dims = 1:20)

### Plots----
p1 <- DimPlot(merge_seurat_fil, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merge_seurat_fil, reduction = 'umap', group.by = 'Type', cols = c('red','green','blue'))


grid.arrange(p1, p2, ncol = 2, nrow = 2)



### perform integration to correct for batch effects ----
object_list <- SplitObject(merge_seurat_fil, split.by = 'Patient')
for(i in 1:length(object_list)){
  object_list[[i]] <- NormalizeData(object = object_list[[i]])
  object_list[[i]] <- FindVariableFeatures(object = object_list[[i]]) 
}

### select integration features----
features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 1000)

# find integration anchors (CCA)

anchors <- FindIntegrationAnchors(object.list = object_list,
                       anchor.features = features)
#integrate data
seurat_integrated <- IntegrateData(anchorset = anchors)