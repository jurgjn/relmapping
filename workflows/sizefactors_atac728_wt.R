library("DESeq2")

fp <- snakemake@input[[1]]
countData <- read.table(fp, header=TRUE, sep="\t")#, row.names=1)

cond <- relevel(factor(c( # Hard-coded(!) unlist(lapply(colnames(countData), function(s) substr(s, 1, nchar(s) - 5)))
    "wt_emb", "wt_emb",
    "wt_l1", "wt_l1",
    "wt_l2", "wt_l2",
    "wt_l3", "wt_l3",
    "wt_l4", "wt_l4",
    "wt_ya", "wt_ya")), "wt_emb")

colData <- data.frame(row.names=colnames(countData), condition=cond)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- estimateSizeFactors(dds)

# sf <- raw DESeq sizefactors
sf <- sizeFactors(dds)

# geomeansMedian <- median counts in pseudo-reference sample
loggeomeans <- rowMeans(log(countData))
geomeansMedian <- exp(stats::median(loggeomeans[is.finite(loggeomeans)]))

# final ATAC-seq normalisation factors
sfScaled <- 1 / (sf * geomeansMedian)

fp <- snakemake@output[[1]]
write.table(sfScaled, fp, sep="\t", quote=FALSE, col.names=F)
