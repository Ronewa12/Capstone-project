library(DESeq2)
the_rawcount <- read.csv("gene_counts_final.txt",header = TRUE, row.names=1, sep = "\t")
the_rawcount  <- the_rawcount[,order(names(the_rawcount))]
the_rawcount <- as.matrix(the_rawcount)

condition <- factor(c(rep('ctrl',20),rep('case',24)))
site <- factor(c(rep('brn',9),rep('lym',11),rep('brn',12),rep('lym',12)))

coldata <- DataFrame(row.names=colnames(the_rawcount), condition, site)

dds <- DESeqDataSetFromMatrix(countData = the_rawcount,
                              colData = coldata,
                              design = ~ condition+site)
featureData <- data.frame(gene=rownames(the_rawcount))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("ctrl","case"))
dds$condition <- droplevels(dds$condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.5)
#Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_case_vs_ctrl", type="apeglm")

library("BiocParallel")
register(MulticoreParam(4))
#We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue), alpha= 0.5]
#How many adjusted p-values were less than 0.5?
sum(res$padj < 0.5, na.rm=TRUE)

#Independent hypothesis weighting
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.5, na.rm=TRUE)
metadata(resIHW)$ihwResult


library(biomaRt)
head(listMarts(), 3) 
head(listDatasets(useMart("ensembl")), 3) 

