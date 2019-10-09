library(DESeq2)
library(data.table)
library(dplyr)
library(biomaRt)

directory = c("~/Documents/R/Samanta.Mariani.RNAseq/counts")
all.files <- list.files(path = directory,pattern = ".txt")


setwd(dir = "Samanta.Mariani.RNAseq/counts/")
mylist <- lapply(all.files, read.table,header=TRUE, row.names=1)
mydata <- as.data.frame(do.call('cbind',mylist))

count_data = cbind(I17.1215.01.AGM1_1_GFPpAPCm.merged.sorted.bam = mydata[,6], 
                   I17.1215.02.AGM1_2_GFPpAPCp.merged.sorted.bam = mydata[,12], 
                   I17.1215.03.AGM2_3_GFPpAPCm.merged.sorted.bam = mydata[,18], 
                   I17.1215.04.AGM2_4_GFPpAPCp.merged.sorted.bam = mydata[,24], 
                   I17.1215.05.AGM3_5_GFPpAPCm.merged.sorted.bam = mydata[,30], 
                   I17.1215.06.AGM3_6_GFPpAPCp.merged.sorted.bam = mydata[,36])


rownames(count_data)=rownames(mydata)


fileNames <- c("AGM1_GFPpAPCm", "AGM1_GFPpAPCp", "AGM2_GFPpAPCm", "AGM2_GFPpAPCp", "AGM3_GFPpAPCm", "AGM3_GFPpAPCp")
sampleCondition <- c("negative","positive", "negative","positive", "negative","positive")


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(count_data), filenames = fileNames, sampleCondition))
dds <- DESeqDataSetFromMatrix(countData=count_data, colData=coldata, design=~sampleCondition)
dds

#Filtering for rows/genes that have 0 or 1
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds$sampleCondition <- factor(dds$sampleCondition, levels=c("positive","negative"))
dds$sampleCondition <- relevel(dds$sampleCondition, ref="negative")

dds <- DESeq(dds)
res <- results(dds, contrast = c("sampleCondition", "positive", "negative"))
res
summary.DESeqResults(res)

# number of significant genes at <= 0.5
sum(res$padj <= 0.05, na.rm=TRUE)
sum(res$log2FoldChange <= -1.5 & res$padj <= 0.05, na.rm = TRUE) # down 170
sum(res$log2FoldChange >= 1.5 & res$padj <= 0.05, na.rm = TRUE) #up 145

gene.table = as.data.frame(res)
gene.table = annotation(gene.table)

sig_gene_table = gene.table[which(gene.table$padj <= 0.05),]
#sig_gene_table  = annotation(sig_gene_table)
sig_gene_table = sig_gene_table[order(sig_gene_table$padj),]
write.table(sig_gene_table, file = "../sig_gene_table.txt", sep = "\t")


upregulated_gene_table = gene.table[which(gene.table$padj <= 0.05 & gene.table$log2FoldChange >= 1.5),]
#upregulated_gene_table  = annotation(upregulated_gene_table)
upregulated_gene_table = upregulated_gene_table[order(upregulated_gene_table$padj),]
write.table(upregulated_gene_table, file = "../upregulated_gene_table.txt", sep = "\t")

downregulated_gene_table = gene.table[which(gene.table$padj <= 0.05 & gene.table$log2FoldChange <= -1.5),]
#downregulated_gene_table  = annotation(downregulated_gene_table)
downregulated_gene_table = downregulated_gene_table[order(downregulated_gene_table$padj),]
write.table(downregulated_gene_table, file = "../downregulated_gene_table.txt", sep = "\t")



counts = assay(dds)
source(file = "../remove.dots.R")
rownames(res) = remove.dots(res)

counts[which(rownames(counts) == "ENSMUSG00000064326"),]

t = data.frame(res)
t = annotation(t)
write.table(t, file = "../diff_ex_table.txt", sep = "\t")
#t$row = rownames(t)



d <- plotCounts(dds, gene= "ENSMUSG00000026712.3", intgroup="sampleCondition",returnData=TRUE)
library("ggplot2")
gg = ggplot(d, aes(x=sampleCondition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0),aes(colour = factor(sampleCondition)))
  
gg + labs(colour = "Samples") + scale_color_manual(labels = c("CD206-negative", "CD206-positive"), values = c("#f8766d","turquoise")) +
  theme_light()+
  theme(panel.background = element_blank(), text = element_text(size=12)) + xlab("") + ylab("Raw gene counts")
  


rld <- rlog(dds, blind=FALSE)
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld), 3)

rld.table = as.data.frame(assay(rld))
rld.table = annotation(rld.table)

# Heatmap of gene expression using the 50 most variable genes
library("pheatmap")
library(genefilter)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE),50)
matrix <- assay(rld)[topVarGenes,]


df <- data.frame(condition = colData(dds)[,c("sampleCondition")])
rownames(df) = colnames(assay(rld)[topVarGenes,])
names(df)[1] = "Samples"


df3 = list(Samples = c(positive = "#4DAF4A", negative="#984EA3"))


pheatmap(mat = matrix, cluster_rows=FALSE, show_rownames=FALSE, 
         show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df, annotation_colors = df3,
         clustering_distance_cols = "correlation", 
         clustering_method = "average", scale = "row",cellwidth=40,cellheight = 4)


# Heatmap of distances between samples
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rld$sampleCondition)
colnames(sampleDistMatrix) <- paste(rld$sampleCondition)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

rownames(sampleDistMatrix) = c("CD206-negative", "CD206-positive", "CD206-negative", "CD206-positive","CD206-negative", "CD206-positive")
colnames(sampleDistMatrix) = c("CD206-negative", "CD206-positive", "CD206-negative", "CD206-positive","CD206-negative", "CD206-positive")

pheatmap(sampleDistMatrix,
         clustering_distance_cols=sampleDists,
         col=colors, show_rownames = TRUE, show_colnames = TRUE)


#PCA plot
library(ggplot2)
pcaData <- plotPCA(rld, intgroup=c("sampleCondition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sampleCondition)) +
  geom_point(size=3) +
  labs(colour = "Samples") + scale_color_manual(labels = c("CD206-negative", "CD206-positive"), values = c("#984EA3","#4DAF4A")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = 2) +
  theme_light()+
  theme(panel.background = element_blank(), text = element_text(size=12))
  
# 
# plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
# {
#   rv <- rowVars(assay(object))
#   select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
#                                                      length(rv)))]
#   pca <- prcomp(t(assay(object)[select, ]))
#   percentVar <- pca$sdev^2/sum(pca$sdev^2)
#   if (!all(intgroup %in% names(colData(object)))) {
#     stop("the argument 'intgroup' should specify columns of colData(dds)")
#   }
#   intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
#   group <- if (length(intgroup) > 1) {
#     factor(apply(intgroup.df, 1, paste, collapse = " : "))
#   }
#   else {
#     colData(object)[[intgroup]]
#   }
#   d <- data.frame(PC1 = pca$x[, 2], PC2 = pca$x[, 3], group = group, 
#                   intgroup.df, name = colData(rld)[,1])
#   if (returnData) {
#     attr(d, "percentVar") <- percentVar[1:2]
#     return(d)
#   }
#   ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
#   
# }
# #view raw


rld.norm = assay(rld)
raw.counts = assay(dds)
plot(rld.norm[,1:2], cex =.1 , main = " Normalized log2 ( read counts ) ")

library(vsn)
par(mfrow=c(1,2))
meanSdPlot (raw.counts , ranks = FALSE , ylim = c(0 ,3) ,main = "sequencing depth normalized log2 ( read counts )" )
meanSdPlot (rld.norm , ranks = FALSE , ylim = c(0 ,3) ,main = "sequencing depth normalized log2 ( read counts )" )


# calculates the correlation between 2 replicates and returns the result as well as a scatterplot
# Arguments : vectors for the replicates
calc_cor = function(replicate1, replicate2){
  #caclulates correlation and prints it on the screen 
  cor =cor.test(as.matrix(replicate1), as.matrix(replicate2), method = "pearson")
  # plots the replicates on a scaterplot
  plot(as.matrix(replicate1), as.matrix(replicate2), xlab = "replicate1", ylab = "replicate2")
  return(cor)
}

#example how to run it for replicate 1 and 3 of the rld matrix
calc_cor(rld[,1], rld[,3])

#example how to run it  and save the correlation in a variable for replicate 1 and 3 of the rld matrix
corr = calc_cor(rld[,1], rld[,3])




annotation = function(table){
  source(file = "../remove.dots.R")
  rownames(table) = remove.dots(table)

  ensembl = useMart(biomart = "ensembl",dataset="mmusculus_gene_ensembl")
  all_genes= getBM(filters = "ensembl_gene_id", attributes= c("external_gene_name","ensembl_gene_id"),
                 values = rownames(table), mart = ensembl)

  test_vlookup=merge(as.data.frame(table), all_genes, by.x="row.names", by.y="ensembl_gene_id") 
  return (test_vlookup)
}

