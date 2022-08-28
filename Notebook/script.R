library(DESeq2)

the.data <- read.csv(file='gene_counts_final.txt', header=TRUE, row.names=1, sep = '\t')
the.data <- the.data[,order(names(the.data))]

the.data <- as.matrix(the.data)
condition <- factor(c(rep('ctrl',20),rep('case',24)))
site <- factor(c(rep('brn',9),rep('lym',11),rep('brn',12),rep('lym',12)))
the.coldata <- DataFrame(row.names=colnames(the.data), condition, site)

dds <- DESeqDataSetFromMatrix(countData=the.data, colData=the.coldata, design = ~ site + condition)
dds$condition <- relevel(dds$condition, ref="ctrl")
dds <- estimateSizeFactors(dds)
dds <- dds[ rowSums(counts(dds,normalized=TRUE) >= 5) >= 3, ]

DE.dds <- DESeq(dds)
res.DE.dds <- results(DE.dds, alpha=0.05, lfcThreshold=1.1, altHypothesis="greaterAbs")


