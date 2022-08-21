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

# # Import data from 10X output
# mat <- readMM("Z:/Van/20210831 - PBMC 10X Prox-seq/alignment_result/filtered_feature_bc_matrix/matrix.mtx.gz")
# feature.names <- read.delim("Z:/Van/20210831 - PBMC 10X Prox-seq/alignment_result/filtered_feature_bc_matrix/features.tsv.gz",
#                             header = FALSE, stringsAsFactors = FALSE)
# barcode.names <- read.delim("Z:/Van/20210831 - PBMC 10X Prox-seq/alignment_result/filtered_feature_bc_matrix/barcodes.tsv.gz",
#                             header = FALSE, stringsAsFactors = FALSE)
# 
# # Rename matrix
# colnames(mat) = barcode.names$V1
# rownames(mat) = feature.names$V2

# Import 10x data
mat <- Read10X("filtered_feature_bc_matrix")

# Histogram
hist(colSums(mat), main="Histogram of UMIs/cell")
hist(colSums(mat>0), main="Histogram of genes/cell")


#*********************************************#
# Discard cells with more than 10,000 UMIs
mat <- mat[, colSums(mat)<=10e3]

# # Discard cells with fewer than 1,000 UMIs
# mat <- mat[, colSums(mat)>=1e3]

# Discard cells with fewer than 500 genes
mat <- mat[, colSums(mat>0)>=500]

# Keep genes detected in at least 5 cells
mat <- mat[rowSums(mat>0)>=5,]

#*********************************************#
#------- use seurat --------
pbmc <- CreateSeuratObject(mat)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features=c('nCount_RNA','nFeature_RNA','percent.mt'))

# Need to separate cells based on percent mitochondria
pbmc <- subset(pbmc, percent.mt<=15)

# Export RNA data
write.table(data.frame(nGene=pbmc$nFeature_RNA, nUMI=pbmc$nCount_RNA), gzfile("RNA_QC.txt.gz"), sep='\t')
# Export cell barcodes
write.csv(colnames(pbmc), "RNA_cell_barcodes.csv")

# Normalize and scale
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features=rownames(pbmc))

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc, ndims=50)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.2) # resolution = 0.2

# UMAP
pbmc <- RunTSNE(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "tsne", label=TRUE, label.size=8, pt.size=0.2) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black')) + NoLegend()
ggsave("figures/rna_tsne.png", dpi=600, width=4.2, height = 3.3, units = "in")
# write.csv(pbmc@reductions$tsne@cell.embeddings, "temp.csv")
# write.csv(Idents(pbmc), "temp.csv")

# Export cluster identity
write.csv(Idents(pbmc), "RNA_cluster_id.csv")

# Cluster markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
# write.csv(pbmc.markers, "rna_markers_resolution0.2.csv")

# celldex PBMC data: https://bioconductor.org/packages/3.14/data/experiment/vignettes/celldex/inst/doc/userguide.html
ref <- celldex::NovershternHematopoieticData()

# Automated annotation
singleR.results <- SingleR(as.SingleCellExperiment(pbmc), assay.type.test = 1,
                           ref = ref, labels = ref$label.main)
table(singleR.results$labels)
plotScoreHeatmap(singleR.results)
pbmc[["singleR.labels"]] <- singleR.results$labels

# For NovershternHematopoieticData()
for (i in c("Basophils","Dendritic cells","Granulocytes","HSCs","Megakaryocytes","NK T cells")) {
  pbmc[["singleR.labels"]][pbmc[["singleR.labels"]] == i] <- "Others"
}
table(pbmc[["singleR.labels"]])
# write.csv(pbmc$singleR.labels, "singleR_labels.csv")

# scRNAseq annotation
# sce <- scRNAseq::KotliarovPBMCData(mode=c("rna"))

# Plot t-SNE
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR.labels", pt.size=0.2, cols=brewer.pal(6,"Dark2")) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black'))
ggsave("figures/singleR_tsne.png", dpi=600,
       width=4.5, height = 3.3, units = "in") # Figure 4a


# Plot genes
FeaturePlot(pbmc, reduction = "tsne", features = c("rna_CD3D","rna_CD4","rna_CD8A","rna_CD79A"))

#------ Import PLA data ------
pla <- read.table("10x_PLA_count_matrix.txt.gz",
                  header=T, row.names=1, sep="\t", check.names = FALSE)

pla2 <- pla[,colnames(pbmc)]

# Total protein abundance
AB1 <- unlist(lapply(rownames(pla2), FUN=function(x) strsplit(x,":")[[1]][1]))
AB2 <- unlist(lapply(rownames(pla2), FUN=function(x) strsplit(x,":")[[1]][2]))
AB <- unique(c(AB1,AB2))
pla2.totalprotein <- data.frame(matrix(0,nrow=length(AB),ncol=ncol(pla2)))
colnames(pla2.totalprotein) <- colnames(pla2)
rownames(pla2.totalprotein) <- AB
for (i in AB)
{
  pla2.totalprotein[i,] <- colSums(pla2[AB1==i,]) + colSums(pla2[AB2==i,])
}

# Add protein data to RNA
pbmc[["protein"]] <- CreateAssayObject(counts = pla2.totalprotein)
DefaultAssay(pbmc) <- "protein"
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin = 2)
# pbmc <- FindVariableFeatures(pbmc)
# pbmc <- ScaleData(pbmc)

# Plot
margin <- theme(plot.margin = unit(c(1,1,1,4), "mm"))
p1 <- FeaturePlot(pbmc, reduction="tsne", features=c("rna_CD3E"), pt.size=0.2) +
  theme_void() + margin
p2 <- FeaturePlot(pbmc, reduction="tsne", features=c("rna_CD4"), pt.size=0.2) +
  theme_void() + margin
p3 <- FeaturePlot(pbmc, reduction="tsne", features=c("rna_CD8A"), pt.size=0.2) +
  theme_void() + margin
p4 <- FeaturePlot(pbmc, reduction="tsne", features=c("rna_CD9"), pt.size=0.2) +
  theme_void() + margin
p5 <- FeaturePlot(pbmc, reduction="tsne", features=c("rna_KLRB1"), pt.size=0.2) +
  theme_void() + margin
fig <- grid.arrange(p1,p2,p3,p4,p5, ncol=5, nrow=1)
ggsave("figures/tsne_5rna.png", plot=fig, dpi=600,
       width=9, height = 1.8, units = "in") # Figure 4c
# write.csv(GetAssayData(pbmc, slot = "data",
#                        assay = "RNA")[c("CD3E","CD4","CD8A","CD9","KLRB1"),],
#           "temp.csv")

p1 <- FeaturePlot(pbmc, reduction="tsne", features=c("protein_CD3"), pt.size=0.2, cols=c("lightgrey","green4")) +
  theme_void() + margin
p2 <- FeaturePlot(pbmc, reduction="tsne", features=c("protein_CD4"), pt.size=0.2, cols=c("lightgrey","green4")) +
  theme_void() + margin
p3 <- FeaturePlot(pbmc, reduction="tsne", features=c("protein_CD8"), pt.size=0.2, cols=c("lightgrey","green4")) +
  theme_void() + margin
p4 <- FeaturePlot(pbmc, reduction="tsne", features=c("protein_CD9"), pt.size=0.2, cols=c("lightgrey","green4")) +
  theme_void() + margin
p5 <- FeaturePlot(pbmc, reduction="tsne", features=c("protein_CD161"), pt.size=0.2, cols=c("lightgrey","green4")) +
  theme_void() + margin
fig <- grid.arrange(p1,p2,p3,p4,p5, ncol=5, nrow=1)
ggsave("figures/tsne_5proteins.png", plot=fig, dpi=600,
       width=9, height = 1.8, units = "in") # Figure 4c
# write.csv(GetAssayData(pbmc, slot = "data",
#                        assay = "protein")[c("CD3","CD4","CD8","CD9","CD161"),],
#           "temp.csv")

#----- Pearson correlation of the remaining proteins -----
margin0 <- theme(plot.margin = unit(c(1,2,1,4), "mm"))
pa <- FeatureScatter(pbmc, feature1 = "rna_CD3E", feature2 = "protein_CD3",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD3E") + ylab("CD3 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
pb <- FeatureScatter(pbmc, feature1 = "rna_CD4", feature2 = "protein_CD4",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD4") + ylab("CD4 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
pc <- FeatureScatter(pbmc, feature1 = "rna_CD8A", feature2 = "protein_CD8",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD8A") + ylab("CD8 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
pd <- FeatureScatter(pbmc, feature1 = "rna_CD9", feature2 = "protein_CD9",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD9") + ylab("CD9 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
pe <- FeatureScatter(pbmc, feature1 = "rna_KLRB1", feature2 = "protein_CD161",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("KLRB1") + ylab("CD161 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p1 <- FeatureScatter(pbmc, feature1 = "rna_TNFRSF9", feature2 = "protein_41-BB",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("TNFRSF9") + ylab("41-BB protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p2 <- FeatureScatter(pbmc, feature1 = "rna_CD80", feature2 = "protein_B7",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD80") + ylab("B7 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p3 <- FeatureScatter(pbmc, feature1 = "rna_CCR5", feature2 = "protein_CCR5",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CCR5") + ylab("CCR5 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p4 <- FeatureScatter(pbmc, feature1 = "rna_CCR6", feature2 = "protein_CCR6",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CCR6") + ylab("CCR6 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p5 <- FeatureScatter(pbmc, feature1 = "rna_CCR7", feature2 = "protein_CCR7",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CCR7") + ylab("CCR7 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p6 <- FeatureScatter(pbmc, feature1 = "rna_IL2RB", feature2 = "protein_CD122",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("IL2RB") + ylab("CD122 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p7 <- FeatureScatter(pbmc, feature1 = "rna_IL7R", feature2 = "protein_CD127",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("IL7R") + ylab("CD127 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p8 <- FeatureScatter(pbmc, feature1 = "rna_IL2RG", feature2 = "protein_CD132",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("IL2RG") + ylab("CD132 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p9 <- FeatureScatter(pbmc, feature1 = "rna_CD14", feature2 = "protein_CD14",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD14") + ylab("CD14 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p10 <- FeatureScatter(pbmc, feature1 = "rna_FCGR3A", feature2 = "protein_CD16",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("FCGR3A") + ylab("CD16 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p11 <- FeatureScatter(pbmc, feature1 = "rna_CD19", feature2 = "protein_CD19",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD19") + ylab("CD19 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p12 <- FeatureScatter(pbmc, feature1 = "rna_IL2RA", feature2 = "protein_CD25",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("IL2RA") + ylab("CD25 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p13 <- FeatureScatter(pbmc, feature1 = "rna_CD27", feature2 = "protein_CD27",
                     group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD27") + ylab("CD27 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p14 <- FeatureScatter(pbmc, feature1 = "rna_CD28", feature2 = "protein_CD28",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD28") + ylab("CD28 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p15 <- FeatureScatter(pbmc, feature1 = "rna_CD38", feature2 = "protein_CD38",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD38") + ylab("CD38 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p16 <- FeatureScatter(pbmc, feature1 = "rna_PTPRC", feature2 = "protein_CD45RA",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("PTPRC") + ylab("CD45RA protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p17 <- FeatureScatter(pbmc, feature1 = "rna_PTPRC", feature2 = "protein_CD45RO",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("PTPRC") + ylab("CD45RO protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p18 <- FeatureScatter(pbmc, feature1 = "rna_CD46", feature2 = "protein_CD46",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CD46") + ylab("CD46 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p19 <- FeatureScatter(pbmc, feature1 = "rna_NCAM1", feature2 = "protein_CD56",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("NCAM1") + ylab("CD56 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p20 <- FeatureScatter(pbmc, feature1 = "rna_SELL", feature2 = "protein_CD62L",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("SELL") + ylab("CD62L protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p21 <- FeatureScatter(pbmc, feature1 = "rna_CTLA4", feature2 = "protein_CTLA-4",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CTLA4") + ylab("CTLA-4 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p22 <- FeatureScatter(pbmc, feature1 = "rna_CXCR3", feature2 = "protein_CXCR3",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CXCR3") + ylab("CXCR3 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p23 <- FeatureScatter(pbmc, feature1 = "rna_CXCR4", feature2 = "protein_CXCR4",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("CXCR4") + ylab("CXCR4 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p24 <- FeatureScatter(pbmc, feature1 = "rna_HLA-DRA", feature2 = "protein_HLADR",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("HLA-DRA") + ylab("HLADR protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p25 <- FeatureScatter(pbmc, feature1 = "rna_ICOS", feature2 = "protein_ICOS",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("ICOS") + ylab("ICOS protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p26 <- FeatureScatter(pbmc, feature1 = "rna_LAG3", feature2 = "protein_LAG3",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("LAG3") + ylab("LAG3 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p27 <- FeatureScatter(pbmc, feature1 = "rna_TNFRSF4", feature2 = "protein_OX40",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("TNFRSF4") + ylab("OX40 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p28 <- FeatureScatter(pbmc, feature1 = "rna_TNFSF4", feature2 = "protein_OX40L",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("TNFSF4") + ylab("OX40L protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p29 <- FeatureScatter(pbmc, feature1 = "rna_PDCD1", feature2 = "protein_PD1",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("PDCD1") + ylab("PD1 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p30 <- FeatureScatter(pbmc, feature1 = "rna_SLAMF6", feature2 = "protein_Slamf6",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("SLAMF6") + ylab("Slamf6 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p31 <- FeatureScatter(pbmc, feature1 = "rna_TIGIT", feature2 = "protein_TIGIT",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("TIGIT") + ylab("TIGIT protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
p32 <- FeatureScatter(pbmc, feature1 = "rna_HAVCR2", feature2 = "protein_TIM3",
                      group.by = "singleR.labels", pt.size=1, cols=brewer.pal(6,"Dark2")) +
  xlab("HAVCR2") + ylab("TIM3 protein") +
  margin0 + theme(legend.position="none", axis.title.x=element_text(face="italic")) + my.theme
# No gamma-delta gene
# p33 <- FeatureScatter(pbmc, feature1 = "rna_TRG", feature2 = "protein_gamma-delta",
#                       group.by = "singleR.labels", pt.size=1, cols=brewer.pal(5,"Set1")) +
#   xlab("TRG mRNA") + ylab("gamma-delta protein") +
#   margin0 + theme(legend.position="none") + my.theme
fig <- grid.arrange(pa,pb,pc,pd,pe,
                    p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
                    p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
                    p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
                    p31,p32, ncol=5)
ggsave("figures/all_proteins_scatter.png", plot=fig, dpi=600,
       width=11.5, height = 17, units = "in") # Supplementary figure 7
# temp <- c("CD3E","CD4","CD8A","CD9","KLRB1","TNFRSF9","CD80","CCR5","CCR6","CCR7",
#           "IL2RB","IL7R","IL2RG","CD14","FCGR3A","CD19","IL2RA","CD27","CD28",
#           "CD38","PTPRC","CD46","NCAM1","SELL","CTLA4","CXCR3","CXCR4","HLA-DRA",
#           "ICOS","LAG3","TNFRSF4","PDCD1","SLAMF6","TIGIT","HAVCR2")
# write.csv(t(GetAssayData(pbmc, slot="data", assay="RNA")[temp,]), "temp.csv")
# temp <- c("CD3","CD4","CD8","CD9","CD161","41-BB","B7","CCR5","CCR6","CCR7","CD122",
#           "CD127","CD132","CD14","CD16","CD19","CD25","CD27","CD28","CD38","CD45RA",
#           "CD45RO","CD46","CD56","CD62L","CTLA-4","CXCR3","CXCR4","HLADR","ICOS",
#           "LAG3","OX40","OX40L","PD1","Slamf6","TIGIT","TIM3")
# write.csv(t(GetAssayData(pbmc, slot="data", assay="protein")[temp,]), "temp.csv")

#------ Add PLA data ------ 
pbmc[["pla"]] <- CreateAssayObject(counts = pla2)
DefaultAssay(pbmc) <- "pla"

# Cluster by PLA count
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin = 2, assay = "pla")
pbmc <- ScaleData(pbmc, assay = "pla")
pbmc <- RunPCA(pbmc, features = rownames(pbmc), reduction.name = "plaPCA", reduction.key = "plaPCA_", verbose = FALSE)
ElbowPlot(pbmc, ndims=50, reduction = "plaPCA")
pbmc<- FindNeighbors(pbmc, dims=1:20, reduction = "plaPCA")
pbmc <- FindClusters(pbmc, resolution = 0.06)
pbmc <- RunTSNE(pbmc, assay = "pla", reduction = "plaPCA", reduction.name = "plaTSNE", reduction.key = "plaTSNE_")
DimPlot(pbmc, reduction = "plaTSNE", group.by = "singleR.labels", pt.size=0.2, cols=brewer.pal(6,"Dark2")) +
  xlab("PLA product t-SNE 1") + ylab("PLA product t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black'))
ggsave("figures/singleR_pla_tsne.png", dpi=600,
       width=4.5, height = 3.3, units = "in") # Figure 4b
# write.csv(pbmc@reductions$plaTSNE@cell.embeddings, "temp.csv")

# Show PLA clustering
margin1 <- theme(plot.margin = unit(c(1,4,1,0), "mm"))
p1 <- FeaturePlot(pbmc, reduction="plaTSNE", features=c("pla_CD3:CD3"), pt.size=0.2, cols=c("lightgrey","coral2")) +
  ggtitle("CD3:CD3") + theme_void() + margin1 + theme(plot.title=element_text(size=10,hjust=0.5))
p2 <- FeaturePlot(pbmc, reduction="plaTSNE", features=c("pla_CD3:CD4"), pt.size=0.2, cols=c("lightgrey","coral2")) +
  ggtitle("CD3:CD4") + theme_void() + margin1 + theme(plot.title=element_text(size=10,hjust=0.5))
p3 <- FeaturePlot(pbmc, reduction="plaTSNE", features=c("pla_CD3:CD8"), pt.size=0.2, cols=c("lightgrey","coral2")) +
  ggtitle("CD3:CD8") + theme_void() + margin1 + theme(plot.title=element_text(size=10,hjust=0.5))
p4 <- FeaturePlot(pbmc, reduction="plaTSNE", features=c("pla_CD9:CD9"), pt.size=0.2, cols=c("lightgrey","coral2")) +
  ggtitle("CD9:CD9") + theme_void() + margin1 + theme(plot.title=element_text(size=10,hjust=0.5))
p5 <- FeaturePlot(pbmc, reduction="plaTSNE", features=c("pla_CD9:CD8"), pt.size=0.2, cols=c("lightgrey","coral2")) +
  ggtitle("CD9:CD8") + theme_void() + margin1 + theme(plot.title=element_text(size=10,hjust=0.5))
fig <- grid.arrange(p1,p2,p3,p4,p5, ncol=5, nrow=1)
ggsave("figures/pla_tsne_5plas.png", plot=fig,
       dpi=600, width=9.1, height = 1.7, units = "in") # Figure 4d
# write.csv(GetAssayData(pbmc, slot = "data",
#                        assay = "pla")[c("CD3:CD3","CD3:CD4","CD3:CD8","CD9:CD9","CD9:CD8"),],
#           "temp.csv")

#------ CD9 interactions ------
# Import annotation by CD9:CD9 PLA product level
cd9.homodimer.annot <- read.csv("CD9_homodimer_cell_types.csv", row.names=1)
cd9.homodimer.annot <- cd9.homodimer.annot[colnames(pbmc),,drop=FALSE] # Rearrange
pbmc$snn.clusters <- Idents(pbmc)
Idents(pbmc) <- cd9.homodimer.annot
pbmc$cd9.homodimer.clusters <- cd9.homodimer.annot
pbmc$cd9.homodimer.clusters[(pbmc$cd9.homodimer.clusters!="CD8_CD9") & (pbmc$cd9.homodimer.clusters!="CD8")] <- "others"
DimPlot(pbmc, reduction = "tsne", group.by="cd9.homodimer.clusters", pt.size=0.4,
        cols=c("gray75","#d95f02","#1b9e77"), order=c("CD8_CD9","CD8","others")) +
  theme_void() + theme(plot.margin = unit(c(0,1,0,1), "mm"))
ggsave("figures/CD9_homo_hetero_tsne.png",
       dpi=600, width=2.8, height = 2.25, units = "in") # Figure 4g
# write.csv(pbmc$cd9.homodimer.clusters, "temp.csv")

# CD8 vs CD8_CD9 markers
cd9.homodimer.protein.markers <- FindMarkers(pbmc, ident.1 = "CD8_CD9", ident.2 = "CD8",
                                     min.pct = 0.25, logfc.threshold = 0.5, assay = "protein")
cd9.homodimer.rna.markers <- FindMarkers(pbmc, ident.1 = "CD8_CD9", ident.2 = "CD8",
                                         min.pct = 0.25, logfc.threshold = 0.5, assay = "RNA")


#----- Import CD9:CD8 complex count -----
cd9.cd8.complex.count <- read.csv("CD9.CD8_complex_count.csv", check.names = FALSE, row.names=1)
# cd9.cd8.complex.count <- cd9.cd8.complex.count[colnames(pbmc),,drop=FALSE]

# Find RNA markers for CD9:CD8-expressing vs non-expressing cells
pbmc.complex <- pbmc[,rownames(cd9.cd8.complex.count)]
temp <- rep("Yes", times = ncol(pbmc.complex))
temp[cd9.cd8.complex.count[,"CD9:CD8"] == 0] <- "No"
Idents(pbmc.complex) <- factor(temp, levels=c("Yes","No"))
cd9.cd8.complex.rna.markers <- FindMarkers(pbmc.complex, ident.1="Yes", ident.2="No",
                                           min.pct = 0.25, logfc.threshold = 0.5, assay="RNA")

# Export RNA markers
write.csv(cd9.cd8.complex.rna.markers, "RNA_markers_CD9.CD8_complex.csv")

# RNA markers
margin2 <- theme(plot.margin = unit(c(0,1,0,1), "mm"))
p2 <- VlnPlot(pbmc.complex, features = c("rna_GZMB"), pt.size=0.5) +
  scale_x_discrete(labels=c("Positive","Negative")) + xlab("") + ylab("Relative expression") +
  ggtitle("GZMB") + 
  my.theme + theme(plot.title=element_text(face="italic")) + margin2
p3 <- VlnPlot(pbmc.complex, features = c("rna_NKG7"), pt.size=0.5) +
  scale_x_discrete(labels=c("Positive","Negative")) + xlab("") + ylab("Relative expression") +
  ggtitle("NKG7") + 
  my.theme + theme(plot.title=element_text(face="italic")) + margin2
p4 <- VlnPlot(pbmc.complex, features = c("rna_CCR7"), pt.size=0.5) +
  scale_x_discrete(labels=c("Positive","Negative")) + xlab("") + ylab("Relative expression") +
  ggtitle("CCR7") +
  my.theme + theme(plot.title=element_text(face="italic")) + margin2
p5 <- VlnPlot(pbmc.complex, features = c("rna_SELL"), pt.size=0.5) +
  scale_x_discrete(labels=c("Positive","Negative")) + xlab("") + ylab("Relative expression") +
  ggtitle("SELL") +
  my.theme + theme(plot.title=element_text(face="italic")) + margin2
fig <- grid.arrange(p2,p3,p4,p5, ncol = 1)
ggsave("figures/CD9.CD8_complex_rna_markers.png", plot=fig,
       dpi=600, width=2.4, height = 8.8, units = "in") # Figure 4h
# write.csv(Idents(pbmc.complex), "temp.csv")
# write.csv(GetAssayData(pbmc.complex, slot = "data", 
#                        assay = "RNA")[c("GZMB","NKG7","CCR7","SELL"),],
#           "temp.csv")

# Protein markers
p1 <- VlnPlot(pbmc.complex, features = c("protein_CCR7"), pt.size=0.5) +
  scale_x_discrete(labels=c("Positive","Negative")) + xlab("") + ylab("Relative expression") + 
  ggtitle("CCR7 protein") +
  my.theme + margin2
fig <- grid.arrange(p1)
ggsave("figures/CD9.CD8_complex_protein_markers.png", plot=fig,
       dpi=600, width=2.4, height = 2.3, units = "in") # Figure 4i
# write.csv(GetAssayData(pbmc.complex, slot = "data", assay = "protein")["CCR7",],
#           "temp.csv")
