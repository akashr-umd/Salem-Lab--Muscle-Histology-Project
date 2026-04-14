library("DEGreport")
library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library(openxlsx)


#import the data
setwd("C:/Users/orsoo/OneDrive/Documents/Aquaculture/Muscle Histology Work/RNA-seq/Filled Area")
gene_counts <- read.csv("Gene_count_filled_area.csv" , check.names = F)
row.names(gene_counts) <- gene_counts$`Gene Name`
gene_counts <- gene_counts[,-1]

metadata <- read.xlsx("Metadata_Table_Filled_Areaa.xlsx")
row.names(metadata) <- metadata$Sample_ID
metadata <- metadata[,2:5]



#Check to see if colnames are shared between gene_counts table and metadata entry
all(colnames(gene_counts[,-1]) %in% row.names(metadata))

#Check to see if colnames are ordered identically between gene_counts table and metadata entry
all(colnames(gene_counts[,-1]) == row.names(metadata))


dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = metadata,
                              design = ~ From + Phenotype)
dds



featureData <- data.frame(gene=rownames(gene_counts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)


dds$Phenotype <- relevel(dds$Phenotype, ref = "Filled_Area_Low")

dds <- DESeq(dds)

res <- results(dds)
res

summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)


plotCounts(dds, gene=which.min(res$padj), intgroup="Phenotype")


degPlot(dds = dds, res = res, n = 6, xs = "Phenotype")


write.csv(as.data.frame(res), 
          file="DESEQ2_batch_corrected_results.csv")

