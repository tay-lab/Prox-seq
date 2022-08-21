# Libraries
library(data.table)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(reshape2)


# Import PLA data
df <- read.table('Jurkat_Raji_plate_count_matrix.txt.gz',
                 header=T, row.names=1, sep="\t", check.names = FALSE)

# Keep cells with at least 500 UMIs
df <- df[, colSums(df)>=500]

# Create factors
cell.type <- character(ncol(df))
cell.type[grepl("Tcells", colnames(df))] <- "T cells"
cell.type[grepl("Bcells", colnames(df))] <- "B cells"
coincubate <- character(ncol(df))
coincubate[grepl("coincubate", colnames(df))] <- "co-incubate"
coincubate[grepl("separate", colnames(df))] <- "separate"

# Total protein abundance
AB1 <- unlist(lapply(rownames(df), FUN=function(x) strsplit(x,":")[[1]][1]))
AB2 <- unlist(lapply(rownames(df), FUN=function(x) strsplit(x,":")[[1]][2]))
AB <- unique(c(AB1,AB2))
df.totalprotein <- data.frame(matrix(0,nrow=length(AB),ncol=ncol(df)))
colnames(df.totalprotein) <- colnames(df)
rownames(df.totalprotein) <- AB
for (i in AB)
{
  df.totalprotein[i,] <- colSums(df[AB1==i,]) + colSums(df[AB2==i,])
}

# Get seprate T cells
df.totalprotein.T <- df.totalprotein[,cell.type=="T cells" & coincubate=="separate"]

# Convert to long format
df.totalprotein.T.long <- melt(as.matrix(t(df.totalprotein.T)), varnames = c("cell.barcodes","protein"), value.name = "UMI")
fig <- ggplot(data=df.totalprotein.T.long[grep("CD3|CD28|PD1|CD147",df.totalprotein.T.long$protein),], aes(x=UMI)) +
  geom_density(aes(fill=factor(protein, levels=c("CD28","PD1","CD147","CD3"))), alpha=0.5, size=0.3) +
  # geom_density(aes(fill=protein), alpha=0.5, size=1) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(0.05,0.8)) +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.5)) +
  scale_x_log10() +
  ylab("Density") +
  guides(fill=guide_legend(title="Target"))
ggsave("Tcells_flowcomparison.pdf", plot = fig, width=4, height = 2.8, units = "in") # Extended Figure 3a
# write.csv(df.totalprotein.T[c("CD28","PD1","CD147","CD3"),], "temp.csv")

# Median
apply(df.totalprotein.T, MARGIN = 1, median)

# Get separate B cells
df.totalprotein.B <- df.totalprotein[,cell.type=="B cells" & coincubate=="separate"]
df.totalprotein.B.long <- melt(as.matrix(t(df.totalprotein.B)), varnames = c("cell.barcodes","protein"), value.name = "UMI")
fig <- ggplot(data=df.totalprotein.B.long[grep("B7|ICAM1|HLADR|PDL1|CD147",df.totalprotein.B.long$protein),], aes(x=UMI)) +
  geom_density(aes(fill=factor(protein, levels=c("B7","PDL1","HLADR","ICAM1","CD147"))), alpha=0.5, size=0.3) +
  # geom_density(aes(fill=protein), alpha=0.5, size=1) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(0.05,0.8)) +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.5)) +
  scale_x_log10() +
  ylab("Density") +
  guides(fill=guide_legend(title="Target"))
ggsave("Bcells_flowcomparison.pdf", plot = fig, width=4, height = 2.8, units = "in")  # Extended Figure 3b
# write.csv(df.totalprotein.B[c("B7","PDL1","HLADR","ICAM1","CD147"),], "temp.csv")

# Median
apply(df.totalprotein.B, MARGIN = 1, median)
