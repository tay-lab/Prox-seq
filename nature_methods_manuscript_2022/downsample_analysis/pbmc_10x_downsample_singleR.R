# Libraries
library(Matrix)
library(dplyr)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
library(cowplot)
library(Seurat)
# library(leiden) # leiden algorithm for clustering
library(SingleR) # annotate PBMC
# library(celldex) # reference gene set


my.theme <- theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
                  axis.text = element_text(size=12, color='black'), plot.title = element_text(size=12, face="plain"),
                  legend.position = "none")

# Load celldex reference
ref <- celldex::NovershternHematopoieticData()

for (j in c("0_1","0_2","0_4","0_6","0_8"))
{
  # Import data from 10X output
  my.folder <- paste0("count_downsample_",j)
  counts <- Read10X(paste0("downsample/cDNA/",my.folder,"/filtered_feature_bc_matrix"))
  
  # # Check number of UMIs per cell
  # hist(colSums(counts), main="Histogram of UMIs/cell")
  
  # Convert to seurat object
  pbmc <- CreateSeuratObject(counts = counts)
  
  #*********************************************#
  #------- use seurat --------
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # VlnPlot(pbmc0.1, features=c('nCount_RNA','nFeature_RNA','percent.mt'))
  
  # Need to separate cells based on percent mitochondria
  pbmc <- subset(pbmc, percent.mt<=15)
  
  # Normalize and scale
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- ScaleData(pbmc, features=rownames(pbmc))
  
  # Automated annotation
  singleR.results <- SingleR(as.SingleCellExperiment(pbmc), assay.type.test = 1,
                             ref = ref, labels = ref$label.main)
  table(singleR.results$labels)
  # plotScoreHeatmap(singleR.results)
  pbmc[["singleR.labels"]] <- singleR.results$labels
  
  # For NovershternHematopoieticData()
  for (i in c("Basophils","Dendritic cells","Granulocytes","HSCs","Megakaryocytes","NK T cells")) {
    pbmc[["singleR.labels"]][pbmc[["singleR.labels"]] == i] <- "Others"
  }
  table(pbmc[["singleR.labels"]])

  # Export singleR for each downsample
  write.csv(pbmc$singleR.labels, paste0("singleR_labels_",j,".csv"))
}
