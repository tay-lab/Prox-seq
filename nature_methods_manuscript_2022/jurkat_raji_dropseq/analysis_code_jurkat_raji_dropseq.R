# Libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(cowplot)
require(gridExtra)
# library(matrix)
library(dplyr)
library(pheatmap)
# library(effsize)
# library(DESeq2)
library(Seurat)
library(RColorBrewer)
library(class)

my.theme <- theme(text=element_text(size=12),
                  axis.text=element_text(size=12,color='black'),
                  plot.title=element_text(size=12,hjust=0.5,face="plain"))

# Saving figure: height of a figure panel should be 300 pixels (3 inches for pdf), then scale it down

# Import PLA data -----
df.neg.pla <- read.table('PLA_dge.txt.gz', sep="\t", header=T, row.names=1)
df.neg.rna <- read.table('cDNA_dge.txt.gz', sep="\t", header=T, row.names=1)

# Remove zero columns
df.neg.rna <- df.neg.rna[,colSums(df.neg.rna)>0]

# Change some protein names
rownames(df.neg.pla) <- gsub("B7 [(]CD80[)]","B7",rownames(df.neg.pla))
rownames(df.neg.pla) <- gsub("ICAM1 [(]CD54[)]","ICAM1",rownames(df.neg.pla))
rownames(df.neg.pla) <- gsub("HLA DR","HLADR",rownames(df.neg.pla))

# Calculate total abundance of proteins -----
AB1 <- unlist(lapply(rownames(df.neg.pla), FUN=function(x) strsplit(x,":")[[1]][1]))
AB2 <- unlist(lapply(rownames(df.neg.pla), FUN=function(x) strsplit(x,":")[[1]][2]))
AB <- unique(c(AB1,AB2))
df.neg.totalprotein <- data.frame(matrix(0,nrow=length(AB),ncol=ncol(df.neg.pla)))
colnames(df.neg.totalprotein) <- colnames(df.neg.pla)
rownames(df.neg.totalprotein) <- AB
for (i in AB)
{
  df.neg.totalprotein[i,] <- colSums(df.neg.pla[AB1==i,]) + colSums(df.neg.pla[AB2==i,])
}

# Analyze with seurat -----
# Remove cells with fewer than 200 genes
df.neg.rna <- df.neg.rna[,colSums(df.neg.rna>0)>=200]
# Remove genes detected in fewer than 5 cells
df.neg.rna <- df.neg.rna[rowSums(df.neg.rna>0)>=5,]
# Remove cells with more than 20% mitochondria genes (without this filter, there will be subclusters with high % mitochondria)
df.neg.rna <- df.neg.rna[,colSums(df.neg.rna[grep("^MT-",rownames(df.neg.rna)),])/colSums(df.neg.rna) <= 0.2]

# Get a list of cell barcodes that exist in both RNA and PLA data set
cell.barcodes <- intersect(colnames(df.neg.rna),colnames(df.neg.totalprotein))

# Export PLA count of matched cell barcodes
# write.table(df.neg.pla[,cell.barcodes], gzfile("RNA_matched_PLA_count_matrix.txt.gz"), sep='\t')
# write.table(df.neg.rna[,cell.barcodes], gzfile("RNA_matched_RNA_count_matrix.txt.gz"), sep='\t')

# Create seurat object
seurat.neg <- CreateSeuratObject(counts=as.sparse(df.neg.rna[,cell.barcodes]))
seurat.neg[["percent.mt"]] <- PercentageFeatureSet(seurat.neg, pattern = "^MT-")
seurat.neg <- NormalizeData(seurat.neg)
seurat.neg <- FindVariableFeatures(seurat.neg)
seurat.neg <- ScaleData(seurat.neg)
seurat.neg <- RunPCA(seurat.neg, verbose = FALSE)
ElbowPlot(seurat.neg, ndims=50)
seurat.neg <- FindNeighbors(seurat.neg, dims=1:10)
seurat.neg <- FindClusters(seurat.neg, resolution=0.1) # 2 clusters
seurat.neg[["rna.idents"]] <- Idents(seurat.neg) # save cluster id from RNA
seurat.neg <- RunTSNE(seurat.neg, dims=1:10)
DimPlot(seurat.neg, reduction ="tsne", cols=brewer.pal(3,"Dark2"), pt.size = 1) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +  my.theme
ggsave("RNA_tsne.pdf", width=3.25, height = 2.8, units = "in") # Fig 1e
# write.csv(seurat.neg@reductions$tsne@cell.embeddings, "temp.csv")

# Annotate cell type
seurat.neg[["rna.cell.type"]] <- "Jurkat"
seurat.neg$rna.cell.type[Idents(seurat.neg)==1] <- "Raji"

# Find gene markers
neg.rna.markers <- FindAllMarkers(seurat.neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
neg.rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Add protein abundance data -----
seurat.neg[["protein"]] <- CreateAssayObject(counts=df.neg.totalprotein[,cell.barcodes])
DefaultAssay(seurat.neg) <- "protein"
seurat.neg <- NormalizeData(seurat.neg, assay="protein", normalization.method = "CLR", margin = 2)
seurat.neg <- ScaleData(seurat.neg, assay = "protein")

margin <- theme(plot.margin = unit(c(1,1,1,4), "mm"))
p1 <- FeaturePlot(seurat.neg, reduction="tsne", features=c("rna_CD3E"), pt.size=1) +
  theme_void() + margin + ggtitle("CD3E") + theme(plot.title = element_text(size=12,hjust=0.5,face="italic"))
p2 <- FeaturePlot(seurat.neg, reduction="tsne", features=c("protein_CD3"), pt.size=1, cols=c("lightgrey","green4")) +
  theme_void() + margin + ggtitle("CD3 protein") + theme(plot.title = element_text(size=12,hjust=0.5))
p3 <- FeaturePlot(seurat.neg, reduction="tsne", features=c("rna_HLA-DRA"), pt.size=1) +
  theme_void() + margin + ggtitle("HLA-DRA") + theme(plot.title = element_text(size=12,hjust=0.5,face="italic"))
p4 <- FeaturePlot(seurat.neg, reduction="tsne", features=c("protein_HLADR"), pt.size=1, cols=c("lightgrey","green4")) +
  theme_void() + margin + ggtitle("HLADR protein") + theme(plot.title = element_text(size=12,hjust=0.5))
fig <- grid.arrange(p1,p2,p3,p4, ncol=4)
ggsave("RNA_vs_protein.pdf", plot = fig, width=10, height = 2.4, units = "in") # Fig 1i
# write.csv(GetAssayData(seurat.neg, slot="data", assay="RNA")[c("CD3E","HLA-DRA"),], "temp.csv")
# write.csv(GetAssayData(seurat.neg, slot="data", assay="protein")[c("CD3","HLADR"),], "temp.csv")

# fig <- VlnPlot(seurat.neg, features=c("protein_CD28","protein_PD1","protein_ICAM1","protein_PDL1"), ncol=4, cols=brewer.pal(3,"Set2"))
# ggsave("4proteins.pdf", plot = fig, width=12.75, height = 3.3, units = "in")

# Cluster on protein data
seurat.neg <- RunPCA(seurat.neg, features = rownames(seurat.neg),
                     reduction.name = "proteinPCA", reduction.key = "proteinPCA_",
                     approx = FALSE, verbose = FALSE)
DimPlot(seurat.neg, reduction = "proteinPCA")
seurat.neg <- FindNeighbors(seurat.neg, reduction = "proteinPCA")
seurat.neg <- FindClusters(seurat.neg, resolution=0.1)
DimPlot(seurat.neg, reduction = "proteinPCA", cols=brewer.pal(3,"Accent"), pt.size = 1) +
  xlab("Protein PCA 1") + ylab("Protein PCA 2") + my.theme
ggsave("protein_pca.pdf", width=3.25, height = 2.8, units = "in") # Fig 1f
# write.csv(seurat.neg@reductions$proteinPCA@cell.embeddings, "temp.csv")
# write.csv(Idents(seurat.neg), "temp.csv")
DimPlot(seurat.neg, reduction = "proteinPCA", cols=brewer.pal(3,"Dark2"), group.by="rna.idents", pt.size = 1) +
  xlab("Protein PCA 1") + ylab("Protein PCA 2") + my.theme + theme(plot.title=element_blank())
ggsave("protein_pca_byRNA.pdf", width=3.25, height = 2.8, units = "in") # Fig 1h
# write.csv(seurat.neg@reductions$proteinPCA@cell.embeddings, "temp.csv")
# write.csv(Idents(seurat.neg), "temp.csv")

# RNA and protein correlation -----
margin <- theme(plot.margin = unit(c(2,2,2,2), "mm"))
p1 <- FeatureScatter(seurat.neg, feature1="rna_CD3E", feature2="protein_CD3",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("CD3E") + ylab("CD3 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p2 <- FeatureScatter(seurat.neg, feature1="rna_CD28", feature2="protein_CD28",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("CD28") + ylab("CD28 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p3 <- FeatureScatter(seurat.neg, feature1="rna_PDCD1", feature2="protein_PD1",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("PDCD1") + ylab("PD1 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p4 <- FeatureScatter(seurat.neg, feature1="rna_HLA-DRA", feature2="protein_HLADR",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("HLA-DRA") + ylab("HLADR protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p5 <- FeatureScatter(seurat.neg, feature1="rna_ICAM1", feature2="protein_ICAM1",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("ICAM1") + ylab("ICAM1 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p6 <- FeatureScatter(seurat.neg, feature1="rna_CD274", feature2="protein_PDL1",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("CD274") + ylab("PDL1 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p7 <- FeatureScatter(seurat.neg, feature1="rna_CD80", feature2="protein_B7",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("CD80") + ylab("B7 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
p8 <- FeatureScatter(seurat.neg, feature1="rna_BSG", feature2="protein_CD147",
                     group.by="rna.cell.type", cols=brewer.pal(3,"Dark2")) +
  xlab("BSG") + ylab("CD147 protein") + margin + my.theme +
  theme(axis.title.x = element_text(face="italic"))
fig <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow=3, ncol=3)
ggsave("RNA_protein_correlation.png", plot=fig, width=9, height=7, units="in", dpi=600) # Supplementary figure 2
# write.csv(GetAssayData(seurat.neg,
#                        slot = "data", assay="protein")[c("CD3","CD28","PD1","HLADR","ICAM1","PDL1","B7","CD147"),],
#           "temp.csv")
# write.csv(GetAssayData(seurat.neg,
#                        slot = "data", assay="RNA")[c("CD3E","CD28","PDCD1","HLA-DRA","ICAM1","CD274","CD80","BSG"),],
#           "temp.csv")

# Add PLA data -----
seurat.neg[["pla"]] <- CreateAssayObject(counts=df.neg.pla[,cell.barcodes])
DefaultAssay(seurat.neg) <- "pla"
seurat.neg <- NormalizeData(seurat.neg, assay="pla", normalization.method = "CLR", margin = 2)
seurat.neg <- ScaleData(seurat.neg, assay = "pla")
seurat.neg <- RunPCA(seurat.neg, features = rownames(seurat.neg), reduction.name = "plaPCA", reduction.key = "plaPCA_", verbose = FALSE)
ElbowPlot(seurat.neg, ndims=50, reduction = "plaPCA")
seurat.neg <- FindNeighbors(seurat.neg, dims=1:10, reduction = "plaPCA")
seurat.neg <- FindClusters(seurat.neg, resolution=0.3)
seurat.neg <- RunTSNE(seurat.neg, assay = "pla", reduction.key = "plaTSNE_", reduction.name = "plaTSNE", reduction = "plaPCA", dims=1:10)
DimPlot(seurat.neg, reduction = "plaTSNE", pt.size=1) +
  xlab("PLA product t-SNE 1") + ylab("PLA product t-SNE 2") + my.theme
ggsave("PLA_tsne.pdf", width=3.25, height = 2.8, units = "in") # Fig 1g
# write.csv(seurat.neg@reductions$plaTSNE@cell.embeddings, "temp.csv")
# write.csv(Idents(seurat.neg), "temp.csv")

# PLA markers
neg.pla.markers <- FindAllMarkers(seurat.neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, assay = "pla")
neg.pla.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

margin <- theme(plot.margin = unit(c(1,1,1,4), "mm"))
p1 <- FeaturePlot(seurat.neg, reduction="plaTSNE", features=c("pla_PD1:CD3"), pt.size=1, cols=c("lightgrey","coral2")) +
  theme_void() + margin + ggtitle("PD1:CD3") + theme(plot.title = element_text(size=12,hjust=0.5))
p2 <- FeaturePlot(seurat.neg, reduction="plaTSNE", features=c("pla_CD3:CD3"), pt.size=1, cols=c("lightgrey","coral2")) +
  theme_void() + margin + ggtitle("CD3:CD3") + theme(plot.title = element_text(size=12,hjust=0.5))
p3 <- FeaturePlot(seurat.neg, reduction="plaTSNE", features=c("pla_ICAM1:HLADR"), pt.size=1, cols=c("lightgrey","coral2")) +
  theme_void() + margin + ggtitle("ICAM1:HLADR") + theme(plot.title = element_text(size=12,hjust=0.5))
p4 <- FeaturePlot(seurat.neg, reduction="plaTSNE", features=c("pla_HLADR:HLADR"), pt.size=1, cols=c("lightgrey","coral2")) +
  theme_void() + margin + ggtitle("HLADR:HLADR") + theme(plot.title = element_text(size=12,hjust=0.5))
fig <- grid.arrange(p1,p2,p3,p4, ncol=4)
ggsave("4PLA_platsne.pdf", plot=fig, width=10, height=2.4, units="in") # Fig 1j
# write.csv(GetAssayData(seurat.neg, slot="data", assay="pla")[c("PD1:CD3","CD3:CD3","ICAM1:HLADR","HLADR:HLADR"),], "temp.csv")
