# Libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(cowplot)
require(gridExtra)
library(dplyr)
library(pheatmap)
library(Seurat)
library(RColorBrewer)

my.theme <- theme(text=element_text(size=12))

# Saving figure: height of a figure panel should be 300 pixels (3 inches for pdf), then scale it down

# Import RNA data -----
df <- read.table('PLA_dge.txt.gz', sep="\t",
                 header=T, row.names=1, check.names = FALSE)

# QC
hist(colSums(df), ylab = "UMI", main = "Number of UMIs")
hist(colSums(df > 0), main = "Number of PLA products")

# Keep cells with at least 10 UMIs and less than 4000 UMIs
df <- df[,(colSums(df)>=10)]

# Make seurat object from PLA data
seurat.pla <- CreateSeuratObject(df, assay = "pla")

# Chi-square test -----
probeA <- unlist(lapply(rownames(df), FUN=function(x) strsplit(x,":")[[1]][1]))
probeB <- unlist(lapply(rownames(df), FUN=function(x) strsplit(x,":")[[1]][2]))

# Chi-square test
chisq.pbmc <- as.data.frame(matrix(0, nrow=nrow(df), ncol=ncol(df)))
colnames(chisq.pbmc) <- colnames(df)
rownames(chisq.pbmc) <- rownames(df)
for (i in colnames(df))
{
  for (j in rownames(df))
  {
    AB.j <- strsplit(j, ":")[[1]]
    my.table <- matrix(c(df[j,i], sum(df[(probeA==AB.j[1])&(probeB!=AB.j[2]),i]),
                         sum(df[(probeA!=AB.j[1])&(probeB==AB.j[2]),i]), sum(df[(probeA!=AB.j[1])&(probeB!=AB.j[2]),i])), ncol=2)
    chisq.pbmc[j,i] <- chisq.test(my.table)$p.value
  }
}

############ Add protein data
# Calculate total abundance of proteins -----
probeA <- unlist(lapply(rownames(df), FUN = function(x) strsplit(x,":")[[1]][1]))
probeB <- unlist(lapply(rownames(df), FUN = function(x) strsplit(x,":")[[1]][2]))
AB <- unique(c(probeA,probeB))
df.totalprotein <- data.frame(matrix(0,nrow=length(AB),ncol=ncol(df)))
colnames(df.totalprotein) <- colnames(df)
rownames(df.totalprotein) <- AB
for (i in AB)
{
  df.totalprotein[i,] <- colSums(df[probeA==i,]) + colSums(df[probeB==i,])
}
seurat.pla[["protein"]] <- CreateAssayObject(counts = df.totalprotein)


############ Cluster cells with protein data
DefaultAssay(seurat.pla) <- "protein"
seurat.pla <- NormalizeData(seurat.pla, normalization.method = "CLR", margin = 2)
seurat.pla <- ScaleData(seurat.pla)
seurat.pla <- RunPCA(seurat.pla, features = rownames(seurat.pla), reduction.name = "protPCA", reduction.key = "protPCA_", verbose = FALSE, approx = FALSE)
ElbowPlot(seurat.pla, reduction = "protPCA", ndims=50)
seurat.pla <- RunTSNE(seurat.pla, assay = "protein", dims=1:10, reduction = "protPCA", reduction.key = "protTSNE_", reduction.name = "protTSNE")
DimPlot(seurat.pla, reduction = "protTSNE", pt.size=1.5) +
  xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
ggsave("Z:/Van/20201202 - Smart-PLAy PBMC PI/data_analysis/figures/protein_tsne.pdf", width = 3.2, height = 2.5, units="in")

margin <- theme(plot.margin=unit(c(0,4,0,4),"mm"))
p1 <- FeaturePlot(seurat.pla, features = c("CD3"), reduction = "protTSNE", cols=c("lightgrey","green4")) +
  xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("CD3 protein") + theme_void() +
  my.theme + margin + theme(plot.title=element_text(size=12, face="plain", hjust=0.5))
p2 <- FeaturePlot(seurat.pla, features = c("CD4"), reduction = "protTSNE", cols=c("lightgrey","green4")) +
  xlab("t-SNE 1") + ylab("") + ggtitle("CD4 protein") + theme_void() +
  my.theme + margin + theme(plot.title=element_text(size=12, face="plain", hjust=0.5))
p3 <- FeaturePlot(seurat.pla, features = c("CD8"), reduction = "protTSNE", cols=c("lightgrey","green4")) +
  xlab("t-SNE 1") + ylab("") + ggtitle("CD8 protein") + theme_void() +
  my.theme + margin + theme(plot.title=element_text(size=12, face="plain", hjust=0.5))
fig <- grid.arrange(p1,p2,p3, ncol=3)
ggsave("CD3CD4CD8_tsne.png",
       plot=fig, width = 7.8, height = 2.5, units="in", dpi=600) # Figure 3a
# write.csv(seurat.pla@reductions$protTSNE@cell.embeddings[,1:2], "temp.csv")
# write.csv(GetAssayData(seurat.pla, slot="data", assay="protein")[c("CD3","CD4","CD8"),], "temp.csv")