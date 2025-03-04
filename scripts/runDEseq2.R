library("DESeq2")
library(ggplot2)
library(dplyr)
library(pheatmap)
library(stringr)
library("tximport")
library("readr")
library("tximportData")

directory <- "/data/Ambs_genomics/liuh/RNAseq155_KenyaStress/CS034879/Analysis_HTCYGDRX2/RSEM"

#sampleFiles <- grep(".genes.results",list.files(directory),value=TRUE)
#names(sampleFiles) <- str_remove(sampleFiles,"_rsem.genes.results")
#data <- tximport(sampleFiles, type = "rsem", txIn = FALSE, txOut = FALSE)
#names(data)
#[1] "abundance"           "counts"              "length"             
#[4] "countsFromAbundance"

#head(data$counts)
#data$counts

sampleFiles <- "RawCountFile_rsemgenes.txt"

txi.rsem <- read.table(file.path(directory,sampleFiles), header=TRUE, check.names = F)
dimnames(txi.rsem)[[1]] <- txi.rsem[,1]

txi.rsem[,1] <- str_split_fixed(txi.rsem[,1],"_",2)[,2]

data <- txi.rsem
str(data)
#sapply(data, mode)
#data2[,2:dim(data2)[2]] <- lapply(data2[,2:dim(data2)[2]], function (x) as.numeric(x))

####Average Row Duplicates
data.a <- do.call(rbind,lapply(lapply(split(data,data$gene_id),`[`,2:ncol(data)),colMeans))

#############
##############

sampleCondition <- c(rep("tumor",27),rep("normal",25))

sampleName = dimnames(data.a)[2]

sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition)


sampleTable$condition <- factor(sampleTable$condition)

write.table(data.a,'RNAseq.Genes_Counts_RSEM.txt', sep='\t', row.names = T, col.names = T, quote = F)  


#Then we build the DESeqDataSet using the following function:
library("DESeq2")

#########
forceMatrixToInteger <- function(m){
    apply (m, c (1, 2), function (x) {
         (as.integer(x))
    })
}

cnt <- forceMatrixToInteger(data.a)

dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = sampleTable,
                              design= ~  condition)


pdf("boxplot.normal_vs_tumor.pdf")
#### boxplot & density plot ----
#add argument las=2 to function boxplot() to make all labels perpendicular to ax
#par(mfrow=c(1,1))
par(cex.axis=0.5)
boxplot(log2(cnt+1),las=2, outline = T)
limma::plotDensities(log2(cnt+1), legend = T)
dev.off()

#####QC
#####
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = sampleTable,
                              design= ~  condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","tumor","normal"))
res

resOrdered <- res[order(res$pvalue),]

summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

### Exploring and exporting results

pdf("MA-plot.filtered.normal_vs_tumor.pdf")
plotMA(res, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
dev.off()

###
pdf("Volcano.filtered.normal_vs_tumor.pdf")
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05)
dev.off()

pdf("PlotCounts.filtered.normal_vs_tumor.pdf")
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

dev.off()

####
####
#####

colData(dds)

# Extracting transformed values
rld <- rlog(dds, blind=FALSE)


# this gives log2(n + 1)
ntd <- normTransform(dds)

pdf("PCA.filtered.normal_vs_tumor.pdf")

plotPCA(ntd, intgroup="condition")
dev.off()

##
pdf("boxplot_cooks.filtered.normal_vs_tumor.pdf")
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

dev.off()


###Dispersion plot and fitting alternatives

pdf("DispersionPlot.filtered.normal_vs_tumor.pdf")
plotDispEsts(dds)

dev.off()

mcols(res)$description
#############################
#############################
###############################
#############################
#############################
###############################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = sampleTable,
                              design= ~  condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","tumor","normal"))
res

## Log fold change shrinkage for visualization and ranking
#resultsNames(dds)

# Using parallelization
library("BiocParallel")
register(MulticoreParam(4))

#p-values and adjusted p-values

resOrdered <- res[order(res$pvalue),]

summary(res)

#sum(res$padj < 0.1, na.rm=TRUE)
sum(res$pvalue < 0.05, na.rm=TRUE)

# Exporting results to CSV files
write.csv(as.data.frame(resOrdered),
          file="condition_tumor_vs_control_normal_vs_tumor.csv")


resSig <- subset(resOrdered, pvalue < 0.05)
resSig

write.csv(as.data.frame(resSig),
          file="condition_tumor_vs_control_p0.05.filtered_normal_vs_tumor.csv")

###
pdf("Volcano_tumor_vs_control.filtered.normal_vs_tumor.pdf")
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05)
dev.off()

###################
######

