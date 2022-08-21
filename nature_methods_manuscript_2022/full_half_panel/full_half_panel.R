# Libraries
library(Matrix)
library(dplyr)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
library(cowplot)
library(Seurat)


my.theme <- theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
                  axis.text = element_text(size=12, color='black'), plot.title = element_text(size=12, face="plain"),
                  legend.position = "none")

#------ Analyze 10x data ------
tenx <- Read10X(data.dir = "filtered_feature_bc_matrix")
dim(tenx)

# Number of genes and UMIs
colSums(tenx>0) %>% median()
colSums(tenx) %>% median()

# Number of UMIs and genes
hist(colSums(tenx), breaks=100)
hist(log10(colSums(tenx)), breaks=100)
hist(colSums(tenx>0), breaks=100)

# Keep cells with UMI >= 2000 and < 100e3
tenx <- tenx[,(colSums(tenx)>=2000) & (colSums(tenx)<100e3)]
# Keep cells with number of genes < 8000
tenx <- tenx[,colSums(tenx>0)<8000]
# Keep genes detected in at least 10 cells
tenx <- tenx[rowSums(tenx>0)>=10,]
dim(tenx)

# Cluster
tenx <- CreateSeuratObject(tenx, assay = "RNA")
tenx[["percent.mt"]] <- PercentageFeatureSet(tenx, pattern = "^MT-")
VlnPlot(tenx, features=c('nCount_RNA','nFeature_RNA','percent.mt'))

# Filter cells based on percent mitochondria
tenx <- subset(tenx, percent.mt<=20)

# Normalize
tenx <- NormalizeData(tenx)
tenx <- FindVariableFeatures(tenx, selection.method = "vst", nfeatures = 2000)
tenx <- ScaleData(tenx, features=rownames(tenx))

# PCA
tenx <- RunPCA(tenx, features = VariableFeatures(object = tenx), npcs = 50)
FeaturePlot(tenx, features="nCount_RNA", reduction = "pca")
ElbowPlot(tenx, ndims=50)

# Clustering
tenx <- FindNeighbors(tenx, dims = 1:10)
tenx <- FindClusters(tenx, resolution = 0.2) # resolution = 0.2

# t-SNE: before doublet removal
tenx <- RunTSNE(tenx, dims = 1:10)
DimPlot(tenx, reduction = "tsne", label=TRUE, label.size=8, pt.size=0.2) +
  xlab("mRNA t-SNE 1") + ylab("mRNA t-SNE 2") +
  theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
        axis.text = element_text(size=12, color='black')) + NoLegend()

# Jurkat and Raji markers ====> cluster 5 is doublet
FeaturePlot(tenx, features=c("CD3D","CD3E","HLA-DRA","CD19"))

# QC metrics
FeaturePlot(tenx, features=c("nCount_RNA","nFeature_RNA","percent.mt"))

# Cluster 5: doublet -> remove
tenx <- tenx[, Idents(tenx) != 5]

# Annotate
tenx$cell.type <- "Jurkat"
tenx$cell.type[Idents(tenx) %in% c(0,4)] <- "Raji"

# Plot: after doublet removal
margin0 <- theme(plot.margin = unit(c(1,1,1,4), "mm"))
p1 <- DimPlot(tenx, reduction = "tsne", group.by = "cell.type", cols=brewer.pal(3,"Dark2"),
              label=FALSE, pt.size=0.2) +
  ggtitle("") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5))
p2 <- FeaturePlot(tenx, reduction="tsne", features=c("rna_CD3E"), pt.size=0.2) +
  ggtitle("CD3E mRNA") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5))
p3 <- FeaturePlot(tenx, reduction="tsne", features=c("rna_HLA-DRA"), pt.size=0.2) +
  ggtitle("HLA-DRA mRNA") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5))
fig <- grid.arrange(p1,p2,p3, ncol=3, nrow=1)

# Add PLA data -----
pla <- read.table("10x_PLA_count_matrix.txt.gz",
                  header=T, row.names=1, sep="\t", check.names = FALSE)
dim(pla)

# Number of PLA products and UMIs
colSums(pla>0) %>% median()
colSums(pla) %>% median()

# Number of UMIs and PLA products
hist(log10(colSums(pla)), breaks=100)
hist(colSums(pla>0), breaks=100)

# Keep cells with UMIs > 50 and < 10,000
pla <- pla[,(colSums(pla)>50) & (colSums(pla)<10e3)]
# Keep PLA products detected in at least 1 cells
pla <- pla[rowSums(pla>0)>=1,]
dim(pla)

# Filter cells that intersect with 10x RNA
rna.bc <- colnames(tenx)
pla.bc <- colnames(pla)
both.bc <- intersect(rna.bc, pla.bc)
length(both.bc)

# Keep intersect barcodes
pla <- pla[,both.bc]
tenx <- tenx[,both.bc]

# Add PLA products
tenx[["pla"]] <- CreateAssayObject(counts = pla)
DefaultAssay(tenx) <- "pla"

# Normalize
tenx <- NormalizeData(tenx, normalization.method = "CLR", margin = 2, assay = "pla")
tenx <- ScaleData(tenx, assay = "pla")
tenx <- RunPCA(tenx, features = rownames(tenx),
               reduction.name = "plaPCA", reduction.key = "plaPCA_", verbose = FALSE)
ElbowPlot(tenx, ndims=50, reduction = "plaPCA")

# Cluster on PLA data
tenx <- FindNeighbors(tenx, dims=1:10, reduction = "plaPCA")
tenx <- FindClusters(tenx, resolution = 0.1)
tenx <- RunTSNE(tenx, assay = "pla", reduction = "plaPCA", reduction.name = "plaTSNE", reduction.key = "plaTSNE_")

# Annotate by PLA
tenx$cell.type.pla <- "Jurkat_H"
tenx$cell.type.pla[Idents(tenx)==3] <- "Jurkat_F"
tenx$cell.type.pla[Idents(tenx)==2] <- "Raji_H"
tenx$cell.type.pla[Idents(tenx)==1] <- "Raji_F"
# tenx$cell.type.pla <- as.factor(tenx$cell.type.pla)

# Show cluster identity
margin0 <- theme(plot.margin = unit(c(1,1,1,4), "mm"))
p1 <- DimPlot(tenx, reduction = "plaTSNE", group.by = "cell.type.pla",
              label=FALSE, pt.size=0.2) +
  ggtitle("") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5))
p2 <- FeaturePlot(tenx, reduction="plaTSNE", features=c("rna_CD3E"), pt.size=0.2) +
  ggtitle("CD3E") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5,face="italic"))
p3 <- FeaturePlot(tenx, reduction="plaTSNE", features=c("rna_HLA-DRA"), pt.size=0.2) +
  ggtitle("HLA-DRA") + theme_void() +
  margin0 + theme(plot.title=element_text(size=10,hjust=0.5,face="italic"))
fig <- grid.arrange(p2,p3, ncol=2, nrow=1)
ggsave("figures/PLA_tsne.png", plot=p1, dpi=600, 
       width=3.3, height=2.5, units = "in") # Extended figure 5b
ggsave("figures/PLA_tsne_mRNA.png", plot=fig, dpi=600,
       width=5.5, height=2.4, units = "in") # Extended figure 5c
write.csv(tenx$cell.type.pla, "PLA_cluster.csv")
# write.csv(tenx@reductions$plaTSNE@cell.embeddings, "temp.csv")
# write.csv(GetAssayData(tenx, slot="data", assay="RNA")[c("CD3E","HLA-DRA"),], "temp.csv")

# Violin plot
margin2 <- theme(plot.margin = unit(c(0,0,4,4), "mm"))
my.theme2 <- theme(axis.title.y = element_text(size = 12), axis.text.x = element_text(angle=45),
                   axis.title.x = element_blank(),
                   axis.text = element_text(size=12, color='black'),
                   plot.title = element_text(size=12, face="plain"),
                   legend.position = "none")
p1 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_CD28:CD28", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("CD28:CD28") +
  margin2 + my.theme2
p2 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_PD1:PD1", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("PD1:PD1") +
  margin2 + my.theme2
p3 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_CD3:CD3", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("CD3:CD3") +
  margin2 + my.theme2
p4 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_HLADR:HLADR", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("HLADR:HLADR") +
  margin2 + my.theme2
p5 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_PDL1:PDL1", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("PDL1:PDL1") +
  margin2 + my.theme2
p6 <- VlnPlot(tenx, group.by="cell.type.pla", features="pla_ICAM1:ICAM1", pt.size=0.15) +
  xlab("") + ylab("Relative level") + ggtitle("ICAM1:ICAM1") +
  margin2 + my.theme2
fig <- grid.arrange(p1,p4,p2,p5,p3,p6, ncol=2, nrow=3)
ggsave("figures/PLA_violin.png", plot=fig, dpi=600,
       width=6, height=6.5, units = "in") # Extended figure 5e
# write.csv(
#   GetAssayData(tenx, slot="data", assay="pla")[c("CD28:CD28","PD1:PD1","CD3:CD3","HLADR:HLADR","PDL1:PDL1","ICAM1:ICAM1"),],
#   "temp.csv"
# )

# Cluster markers
tenx$pla.idents <- Idents(tenx)
Idents(tenx) <- tenx$cell.type.pla
pla.markers <- FindAllMarkers(tenx, test.use="wilcox",
                              assay="pla", only.pos=TRUE)

# Top 10 markers per cluster
top10 <- pla.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
fig <- DoHeatmap(tenx, features=top10$gene, lines.width=80)
ggsave("figures/PLA_markers_heatmap.png", plot=fig, dpi=600,
       width=5.5, height=7, units = "in") # Extended figure 5f
# write.csv(fig$data, "temp.csv")
