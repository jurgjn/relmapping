library("DESeq2")

fp <- snakemake@input[[1]]
countData <- read.table(fp, header=TRUE, sep="\t", row.names=1)

cond <- relevel(factor(c( # Hard-coded(!) unlist(lapply(colnames(countData), function(s) substr(s, 1, nchar(s) - 5)))
    "lcap808_wt_emb", "lcap808_wt_emb",
    "lcap808_wt_l1", "lcap808_wt_l1",
    "lcap808_wt_l2", "lcap808_wt_l2",
    "lcap808_wt_l3", "lcap808_wt_l3",
    "lcap808_wt_l4", "lcap808_wt_l4",
    "lcap808_wt_ya", "lcap808_wt_ya",
    "lcap808_glp1_d1", "lcap808_glp1_d1",
    "lcap808_glp1_d2", "lcap808_glp1_d2",
    "lcap808_glp1_d6", "lcap808_glp1_d6",
    "lcap808_glp1_d9", "lcap808_glp1_d9",
    "lcap808_glp1_d13", "lcap808_glp1_d13")), "lcap808_wt_emb")

colData <- data.frame(row.names=colnames(countData), condition=cond)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- estimateSizeFactors(dds)

fp <- snakemake@output[[1]]
write.table(sizeFactors(dds), fp, sep="\t", quote=FALSE, col.names=F)
