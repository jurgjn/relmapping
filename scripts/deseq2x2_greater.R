#!/usr/bin/Rscript
suppressMessages(library("DESeq2"))
options(max.print=1000000) # Default limit is aroun 16k rows, so truncates output...
options(width=10000)
countData <- read.table(file('stdin'), header=TRUE, sep="\t", row.names=1)
colData <- data.frame(row.names=colnames(countData),
                      condition = relevel(factor(c("cond1", "cond1", "cond2", "cond2")), "cond1"))
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
#suppressMessages(dds <- DESeq(dds))
dds <- DESeq(dds)
res <- results(dds, altHypothesis="greater")
as.data.frame(res)
